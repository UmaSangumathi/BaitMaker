#!/usr/bin/python
##############################################################
#  Script     : FindSpeciesInSample.py
#  Author     : uma sangumathi
#  Date       : 28/04/2015
#  Last Edited: 24/05/2015, uma 
#  Description: To generate blast XML file to be imported in MEGAN software to view the sample species information
# 
##############################################################
# Purpose: Generates the Blast XML file for all the reads in the fastq file. 
#	Temporary fasta files are created which can be deleted 
#       
# Requirements:
#    1. NCBI blast command line tool
#    2. nt database downloaded from NCBI
#    3. BioPython 
#
# Notes: Please delete the output xml file generated if the script is rerun 	
#		
#############################################################


import os,sys
import glob
from Bio import SeqIO

print "Usage:"
print " ./FindSpeciesInSample.py  'Fastq directory which are quality checked'  'NCBI nt database location' " 

if len(sys.argv)==3:
  fastq_dir = sys.argv[1]
  ncbi_db = sys.argv[2]
  os.chdir(fastq_dir)
  for y in glob.glob("*-qc-t15-l35.fastq"):
    print "Converting fastq file to fasta"
    fasta_fil = open(y.replace(".fastq", "-reads_nAll.fa"), 'a')
    records = SeqIO.parse(open (y, "r"), "fastq")
    c = 0
    for iii in records:
      c = c + 1
      if c <= 100000000: SeqIO.write(iii, fasta_fil , "fasta")
    print "Change the space to '='..."
    os.system("sed 's/ /=/g' " + y.replace(".fastq", "-reads_nAll.fa") + " > " + y.replace(".fastq", "-reads_fa_nAll.fa"))

    blast_cmd2 = "blastn -db " + sys.argv[2]  + "  -query " +  y.replace(".fastq", "-reads_fa_nAll.fa")  + " -outfmt 5 -num_alignments 1  -out " + y.replace(".fastq", "_blast_nAll.xml")
    print "Running: ", blast_cmd2
    if not os.path.isfile(y.replace(".fastq", "_blast_nAll.xml")) : os.system( blast_cmd2)
    else: print "Refuses to overwrite: delete the '.xml' file and rerun"
else:
  print "ERROR: Enter the following paramenters:"
  print " ./FindSpeciesInSample.py  'Fastq directory which are quality checked'  'NCBI nt database location' "


