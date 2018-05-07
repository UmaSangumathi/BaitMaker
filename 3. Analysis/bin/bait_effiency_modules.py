#!/usr/bin/python

import os
import sys
import glob
import datetime
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Entrez
from operator import itemgetter

Entrez.email = "umasang7@gmail.com"

def write_csv(data, file_nam):
    import csv, itertools
    with open(file_nam, 'wb') as f:
        writer = csv.writer(f, delimiter=",")
        print "Witing to :", file_nam
        writer.writerows(list(itertools.chain(t)) for t in data)

  
# Get the fragment sequence ffrom the pair end reads
def get_fragment_seqs(frag_pos_file, ref_file):
  Refs = SeqIO.parse(ref_file, "fasta")
  for r in Refs: 
    ref = str(r.seq)
  frags = []
  with open(frag_pos_file, "r") as f:
    X = f.readlines() 
  for x in X:
    y = x.strip().split(",")  
    if len(y) == 6:
      frags.append([y[0], y[3], int(y[3]) + int(y[1]), ref[int(y[3])-1 : int(y[3]) + int(y[1])-1]]) 
  return frags

# Selects the baits that might actually contribute in the enrichment process
def proccess_blasthits(Genome_Strand, hits, baits):
  header = ["S.No","Bait Id", "Bait sequence", "B.start", "B.end", "Ref.start", "Ref.end", "% Identity", "Len.covered", "Arbitrary.Ref.start", "Arbitrary.Ref.end", "Binding flag"]  
  new_map = {}
  for h in hits:
    hh = h.split("\t")
    if len(baits[hh[0]]) < 4 :
      print hh
      # For Dengue baits the mapping is plus/minus strand
      if Genome_Strand == "Negative":
        if int(hh[5]) > int(hh[4]) :          
          arb_end = int(hh[5]) + abs(int(hh[2]) - 120)
          arb_start = int(hh[4]) - abs(int(hh[3]) - 1)    
          coverage = int(hh[5]) - int(hh[4]) +1 
        else: 
          arb_end = "-"
          arb_start = "-"
          coverage = "-"
      elif Genome_Strand == "Positive":
        if int(hh[5]) < int(hh[4]):
          arb_start = int(hh[5]) - abs(int(hh[3]) - 120)
          arb_end = int(hh[4]) + abs(int(hh[2]) - 1)
          coverage = int(hh[4]) - int(hh[5]) +1
        else:
          arb_end = "-"
          arb_start = "-"
          coverage = "-"
      if coverage >= 70 or float(hh[6]) < 70:
        if coverage >= 110 and float(hh[6]) > 85:
	  tag = "V.High"
        if coverage >= 100 and float(hh[6]) > 80:
          tag = "High"
        elif coverage >= 80  and float(hh[6]) > 75:
          tag = "Medium"
        else: tag = "Low"
      else: tag = "NA"
      if coverage == "-":  tag = "Fail"
      a = [hh[0], baits[hh[0]][2], hh[2], hh[3], hh[4], hh[5], hh[6], coverage ,  arb_start , arb_end , tag]
      if hh[0] not in new_map: new_map[hh[0]] = a 
  for x in baits:
    if x not in new_map: new_map[x] = baits[x]
    else: new_map[x].insert(0, baits[x][0])
  Bait_map = [new_map[x] for x in new_map]
  Bait_map_sort = sorted(Bait_map, key=itemgetter(0))
  Bait_map_sort.insert(0, header)
  return Bait_map_sort

