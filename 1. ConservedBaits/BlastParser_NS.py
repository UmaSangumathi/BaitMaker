#!/usr/bin/python

##############################################################
#  Script     : BlastParser_NSpy
#  Author     : uma sangumathi
#  Date       : 28/06/2013
#  Last Edited: 
#  Description: design 120mer bait 
##############################################################
# Purpose: 
#  check if designed baits have any offtarget hits (human, mouse, bacteria)
#  Blast NCBI
#
# Requirements:

##############################################################

import os
import itertools


# ---------------------- Run Blast --------------
def runBlastn(vir_nam, blast_infile, blast_dir, refseq_fil, Ref_id):
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML

    blast_n = "/home/rongli/src/ncbi-blast-2.2.31+/bin/blastn"
    # blast_db = "/home/rongli/software/nt"
    tab_fmt = "'7 qseqid sseqid qstart qend sstart send pident qcovs evalue bitscore stitle'"
    search_org = [Ref_id + '-Refid', 'txid9606-human', 'txid10114-rat', 'txid9397-bats', 'txid9031-chicken']
    print "Changed virus:", vir_nam
    # 'txid2-bacteria', can be added back to filter for bacteria in viruses


# Check for multiple binding sites
    for txid_nam in search_org:
        [txid, nam] = txid_nam.split('-')
        if nam == "Refid":
            blast_uni = blast_dir + "/" + vir_nam + "--blast--Refid.txt"
            print "Ref:", txid + '[uid] -- Checking for multiple binding sites'
            blast_run_uni = NcbiblastnCommandline(cmd=blast_n, task='blastn', evalue=1.0, query=blast_infile, subject=refseq_fil, outfmt=tab_fmt, out=blast_uni)
            stdout, stderr = blast_run_uni()

# txid9397 bats [ORGN] OR txid9031 chicken [ORGN] OR txid9606 human [ORGN] OR txid10114 rat [ORGN])
        else:
            if nam == "human":
                query_term = '"' + 'Human[PORGN] NOT Viruses[ORGN] NOT Patient* NOT patient NOT patients NOT patient[TEXT]  NOT clone' + '"'
            else:
                query_term = '"' + txid + '[PORGN] NOT Viruses[ORGN] NOT clone NOT virus NOT virion NOT viruses NOT chimeric' + '"'
            blast_out = blast_dir + "/" + vir_nam + "--blast--" + nam + ".txt"
            blast_run = NcbiblastnCommandline(cmd=blast_n, db="nt", task='blastn', evalue=1.0, query=blast_infile, outfmt=tab_fmt, out=blast_out, entrez_query=query_term, remote=True, perc_identity=80.0)
            print "Blast against a subset database: ", query_term, nam
            stdout, stderr = blast_run()

    print "Virus :", vir_nam
    kmer_del = BlastParser(vir_nam, blast_dir, search_org)
    print "::::::::\n", kmer_del
    return kmer_del


# ----------- Parse Blast output ----------------
def BlastParser(vir, blast_dir, search_org):
  res = {}

  def parser(set_file, nam):		# Dictionary from blast op file
    print "Parsing file: ", set_file
    if nam == "Refid":
      PERC_IND = 90.0
      QUER_COV = 85.0  # 102/120 => considered as repeats....
    else:
      PERC_IND = 80.0
      QUER_COV = 40.0  # 48/120 => Failed cut-off
    blast_dict = {}
    with open(set_file, 'r') as f:
     # for key, group in itertools.groupby(f, lambda line: line.startswith('# BLASTN 2.2.28+')):
      for key, group in itertools.groupby(f, lambda line: line.startswith('# BLASTN 2.2.31+')):

        list_Q = [] 		# initiate tmp list
        if not key:
          for line in group:
            list_Q.append(line.strip('\n'))
          query = list_Q[0]
          for i in range(5, len(list_Q)):
            d = list_Q[i].split('\t')
            if len(d) != 1:
              q_coverage = (int(d[3])-int(d[2]))/1.20 		# blast site query length *100/120
              percent_ind = float(d[6])
              if q_coverage > QUER_COV and percent_ind > PERC_IND:
                if query in blast_dict:
                  blast_dict[query].append(d)
                else:
                  blast_dict[query] = ([d])
    return blast_dict

  for org in search_org:
    [txd, nam] = org.split('-')
    kmer_delete = []
    set_file = blast_dir + "/" + vir + "--blast--" + nam + ".txt"
    print set_file
    if os.path.exists(set_file):
      hits = parser(set_file, nam)
      if nam == "Refid":
        HIT_CUTOF = 1
      else:
        HIT_CUTOF = 4
      for key in hits.keys():
        if len(hits[key]) > HIT_CUTOF:		# if hits greater than 3
          q_detail = key.split('|')
          print q_detail
          loc = int(q_detail[2])
          filt = q_detail[1]
          typ = (q_detail[0].split(':')[1]).split('-')[0]
          q_det = (typ + '-' + filt).replace(' ', '')
          print "Query delete", q_det, loc
          kmer_delete.append((loc, q_det, hits[key]))
      res[nam] = kmer_delete
      # print "Resnam", res[nam]
    else:
      print "No file found ERROR !!!"
      exit
  return res, search_org


# ------------ combine the blast results with baits
def baitBlast(virus_nam, blast_infile, All_bait, blast_dir, refseq_fil, Ref_id):
  print "\n-------- Run Blastn ----------"
  Final_baits = []
  [kmer_del, search_org] = runBlastn(virus_nam, blast_infile, blast_dir, refseq_fil, Ref_id)
  print All_bait[0], len(All_bait[0])
  for bait in All_bait:
    coln = len(bait)
    for org_O in search_org:
      [txd, org] = org_O.split('-')
      coln = coln + 2
      if len(kmer_del.keys()) != 0:
        for det in kmer_del[org]:
          if int(det[0]) == int(bait[0]):
            bait = bait + ('Failed', det[2])
            print "failed bait", bait
        if len(bait) == coln - 2:
          bait = bait + ('Passed', "Null")
          print "passed bait", bait
      else:
        if len(bait) == coln - 2:
          bait = bait + ('Passed', "Null")
          print "passed bait", bait
    Final_baits.append(bait)
  return Final_baits
