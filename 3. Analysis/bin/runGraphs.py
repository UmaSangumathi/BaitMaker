#!/usr/bin/python

import os
import sys
import glob
import datetime
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Entrez
import fileinput

hmdir = "/mnt/sdc/uma-sdc/RotV"
wrdir = hmdir + "/Graphs"

# Read the summary file:
info_dict = {}
with open(hmdir + "/run_summary.csv", 'r') as csvfile:
  reader = csv.reader(csvfile)
  for xi in reader:
    print xi
    info_dict[xi[0]] = xi[1]
print info_dict.keys()

for x in info_dict:
  fil_nam = x.replace(".fq", "-") + "hybrid"
  ref_nam = info_dict[x].replace(".fa", "")
  print fil_nam, ref_nam
  if "capture" in fil_nam:
    library_type = "Capture" 
  elif "Cap" in fil_nam:
    library_type = "Capture"
  else: library_type = "Library"
  print  library_type, fil_nam
#  cmd = "/mnt/sdc/uma-sdc/Enrichment_Paper/generated_circosPlot.py " + fil_nam + ".cov  " + fil_nam + "-snps.vcf " + fil_nam + ".fa " + ref_nam + "_baits.fa  " + ref_nam + ".genes " + ref_nam + ".coords " + library_type
  cmd = hmdir + "/generated_circosPlot.py " + fil_nam + ".cov  " + fil_nam + "-snps.vcf " + fil_nam + ".fa " + "Rotv_baits.fa  " + ref_nam + ".genes " + ref_nam + ".coords " + library_type

  print cmd
  os.system(cmd)
  
