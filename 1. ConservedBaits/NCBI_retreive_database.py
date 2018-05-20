#!/usr/bin/python

##############################################################
#  Script     : NCBI_retreive_database.py
#  Author     : Uma Sangumathi
#  Date       : 28/06/2013
#  Last Edited: 
#  Description: Download sequences from NCBI for query viruses
##############################################################
# Purpose: 
#  Search NCBI for '[ORGANISM] AND 1000:120000000[SLEN] NOT clone NOT vector'  # Virus and specific length
#  Create sqlite3 database 
#  Update the database  with hits
# Requirements:
#  1. NCBI blast commandline tool and database
#  2. Sqlite 

#############################################################


import os
import Bio
import pickle
from Bio import Entrez
from datetime import datetime
import collections
import operator
import time
from urllib2 import HTTPError
import httplib 
from httplib import HTTPException
import sqlite3 
import sys
import types


Entrez.email = "tmp@gmail.com"			#Log in
hm_dir = "/home/gmsusk/Project/Indentity_per_baits"
con1 = None
con1 = sqlite3.connect(hm_dir + '/db_VirusFiles/VirusCapture.db')
curi = con1.cursor()
#curi.execute("DROP TABLE IF EXISTS VirusData_4sp")
#curi.execute("CREATE TABLE VirusData_4sp(Species TEXT, Unique_ID INT Primary KEY, Primary_accession TEXT, Sequence_length INT, Other_Seq_id TEXT, Taxonomy TEXT,Source TEXT,Molecular_type TEXT,Comment TEXT, Host TEXT, Strain TEXT, Collection_Date TEXT, dbRef TEXT, Isolate TEXT, Country TEXT, Journal TEXT, Title TEXT, Sequence TEXT )")
#cur.execute("DROP TABLE IF EXISTS Hits")
#cur.execute("CREATE TABLE Hits (Virus TEXT PRIMARY KEY, hitGrt1000 INT, hitsGrt100 INT, ReferenceGenome TEXT, UniqueID TEXT)")


""" Fn: Write into file! """
def write_file(file_name, data_to_write):                       # Write dictionary or list to file
    with open(file_name, 'wb') as handle:
        pickle.dump(data_to_write, handle)


""" Fn: Read a file """
def read_file(file_name):                                       # Read Dict or list from file
    with open(file_name, 'rb') as handle:
        x = pickle.load(handle)
        return x

""" Fn: Check whether the term is in Entrez database"""
def check_entry(rec, R, s, query):
    xi = "NULL"
    if rec in R:
        for i in range(0,len(s)):
            if s[i]['GBQualifier_name']==query:
                xi = s[i]['GBQualifier_value']
    return xi


""" Fn: Sort NCBI hits by Sequence length """
def sort_by_seqlen(search_word, virus):				# Retrieve all the records and sort it by Seq length
    Ref_gen = {}
    info_dict = {}
    sear = Entrez.esearch(db="nucleotide", term=search_word, retMax=10000)   # Maximum allowed hits to be retrieved =10000
    hits = Entrez.read(sear)
    count = int(hits['Count'])
    if int(hits['Count'])<20:
        sear = Entrez.esearch(db="nucleotide", term=virus + '[ORGANISM] AND 100:120000000[SLEN]', retMax=10000)   # Maximum allowed hits to be retrieved =10000
        hits = Entrez.read(sear)

    count100 = int(hits['Count'])
    print "Count: ", hits["Count"]				# Display the num of hits
    for v_id in hits["IdList"]:
	time.sleep(1)			# Wait 2 sec to ensure not many request is send to Entrez
        try:
	    handle = Entrez.efetch(db="nuccore", id=v_id, retmode='xml', rettype="gb")      # Get information of each hit
        except HTTPError:
	    print "HTTP Error detected in Sorting...."
            time.sleep(25)			# Wait 20 sec and try http again 
	    try:
	        handle = Entrez.efetch(db="nuccore", id=v_id, retmode='xml', rettype="gb")      # Get information of each hit
	    except:
		print "Error!!! sort writing to file:", v_id 
                with open('ncbi_sort.err','a') as f:
                    f.write(str(v_id)+'\n')
		 
		f.close() 
	    continue     
        except HTTPException, e:
    	    print "Network problem: %s" % e
    	    print "Second (and final) attempt..." 
	    try:
	        handle = Entrez.efetch(db="nuccore", id=v_id, retmode='xml', rettype="gb")  
	    except:
	        with open('ncbi_sort.err','a') as f:
                    f.write(str(v_id)+'\n')

                f.close()
	    continue
	except httplib.IncompleteRead:
	    print "Incomplete Read Error!!!"
            try:
                handle = Entrez.efetch(db="nuccore", id=v_id, retmode='xml', rettype="gb")
            except:
                with open('ncbi_sort.err','a') as f:
                    f.write(str(v_id)+'\n')

                f.close()
            continue


        record = Entrez.read(handle, validate=False)
        handle.close()
        info_list = []
	info_list.append(virus )
	info_list.append(int(v_id))
       	info_list.append(record[0]["GBSeq_primary-accession"])
        info_list.append(int(record[0]["GBSeq_length"]))
	info_list.append(str(record[0]["GBSeq_other-seqids"]))	
	other_seq = record[0]["GBSeq_other-seqids"]
        info_list.append(record[0]["GBSeq_taxonomy"])
        info_list.append(record[0]["GBSeq_source"])
        info_list.append(record[0]["GBSeq_moltype"])
	if "GBSeq_comment" in record[0]:
            info_list.append(record[0]["GBSeq_comment"])
	else:
	    info_list.append("NULL")

        info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'host'))
        info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'strain'))
        info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'collection_date'))
        info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'db_xref'))
        info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'isolate'))
	info_list.append(check_entry("GBSeq_feature-table", record[0], record[0]['GBSeq_feature-table'][0]['GBFeature_quals'], 'country'))	
	if 'GBReference_journal' in record[0]['GBSeq_references'][0]:
	    info_list.append(record[0]['GBSeq_references'][0]['GBReference_journal'])
	else:
	    info_list.append("NULL")

        if 'GBReference_title' in record[0]['GBSeq_references'][0]:	
	    info_list.append(record[0]['GBSeq_references'][0]['GBReference_title'])
	else:
	    info_list.append("NULL")

	info_list.append(record[0]["GBSeq_sequence"])
	print "insert table:", v_id
	curi.execute("INSERT or REPLACE INTO NCBI_records VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", info_list)
	con1.commit()
        Ref_genome_id = [e for e in other_seq if e.startswith('ref|N') or e.startswith('ref|A')]
        if len(Ref_genome_id) > 0:
             print "Reference Sequence found!!!"
             Ref_gen['Vir'] = v_id 
	     Ref_gen['Ref'] = Ref_genome_id
             print Ref_gen['Ref']

    v = Ref_gen.get('Vir',  "NULL")	
    r = Ref_gen.get('Ref', "NULL" )
    def check(x):
	if len(x)<1:
	    x = "NULL"
	    return x
	else:
	    if (isinstance(x, types.ListType))==True:
	        return (",".join(x))
	    else:
		return (x)
    Ref_gen = {}

""" Main CODE """
os.chdir(hm_dir)
#NCBI_hit_virus('Virus_Spieces.txt')
#Retrieve_data('[ORGANISM] AND 1000:1200000[SLEN]', 'Final_Species_list.txt')

#list_name = 'Final_Species_list.txt'
#list_name = 'error_sort.bad.txt'
#virus_list = [line.rstrip('\r\n') for line in open(list_name, 'r')] # Create the virus list from file [*** \r\n for MAC and \n for windows]
add_search_term = '[ORGANISM] AND 1000:120000000[SLEN] NOT clone NOT vector'
hit_count = []
Ref_genome_id = []
genome_dict = {}
for virus in [ 'Dengue virus 1', 'Dengue virus 2', 'Dengue virus 3', 'Dengue virus 4']:
#    virus = virus.replace('\"','')
    curi.execute("Select Species from NCBI_records WHERE Species=?", (virus,))
    Vir = curi.fetchall()
#    print Vir
    if len(Vir)!=1: 
        print virus , 'is missing!!!'
        print 'Updating the database....', virus
        search_virus = virus + add_search_term
	print 'search word', search_virus 
        sort_count = sort_by_seqlen(search_virus, virus)                # Calls Fn: sort_by_seqlen
    else:
        print virus, ' present in database'
con1.commit()
