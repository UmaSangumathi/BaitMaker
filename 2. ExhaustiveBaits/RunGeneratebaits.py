#!/usr/bin/python

##############################################################
#  Script     : RunGeneatebaits.py
#  Author     : Uma sangumathi
#  Date       : 26/10/2015
#  Last Edited: 
#  Description: Wrapper to Design all possible baits for the viral species and all strains
##############################################################
# Purpose: 
# Iteratively design baits 
# Requirements:
# CD-HIT
# blast
##############################################################

import os
import sys
from Bio import SeqIO
import glob
import pickle
import difflib
from difflib import SequenceMatcher
from Bio.Blast.Applications import NcbiblastnCommandline
import time
import itertools
import csv

from DesignMinimalBaits import *


def summary(infile, renamfile, wind, reffile):
  all_baits = {}
  for i in SeqIO.parse(open(infile, "r"), "fasta"):
    all_baits[i.id] = str(i.seq)
  with open(renamfile, 'w') as f:
    n = 0
    for i in all_baits:
      n = n + 1
      if len(str(all_baits[i]).replace('N', '')) == 120:
        f.write(">bait-" + str(n) + "\r\n" + str(all_baits[i]) + "\r\n")
  [bait_loc, bait_nam, ref_seq] = blastAgainstSeq(reffile, renamfile)
  [untar_csv, map_csv] = summaryRegionPos(ref_seq, bait_loc, bait_nam, wind, renamfile)
  return untar_csv


# -------------------------------------------------------------------------
if len(sys.argv) == 4:
  print "Designing customized baits : Started on ", time.strftime("%c")
  print str(sys.argv)
  direc = str(sys.argv[1])
  initial_baits_f = str(sys.argv[2])
  orgref = str(sys.argv[3])
  os.chdir(direc)
  wind = 500  # untargetted window default 500
  delete = "TRUE"
  if delete == "TRUE":
    for filePath in glob.glob("tmp*"):
      if os.path.isfile(filePath):
        os.remove(filePath)
    for filePath in glob.glob("iter*"):
      if os.path.isfile(filePath):
        os.remove(filePath)

  # Check minimum ref sequence length as 500
  ref = orgref.replace(".fa", "-curated.fa")
  with open(ref, "w") as f:
    Ref = SeqIO.parse(open(orgref, "r"), "fasta")
    for r in Ref:
      if len(str(r.seq)) >= 500:  # default 500
        f.write(">" + r.id + "\r\n" + str(r.seq) + "\r\n")

  # Removes baits in the initial baits file if the similarity is greater than 90%
  print "I. Map the baits given to the reference genomes...."
  [bait_loc, bait_nam, ref_seq] = minimizeBait(initial_baits_f, "tmp1_minimal.fa", ref)

  print "IIa. Get Coverage Baits for the regions the baits are not present..."
  [untar_csv, map_csv] = summaryRegionPos(ref_seq, bait_loc, bait_nam, wind, "iter1.fa")
  untar_reg = genomesUntargeted(untar_csv, ref, "iter-1_untargeted_portion.fa")
  for i in ref_seq:
    if i in untar_reg:
      coverageBaits("iter1", untar_reg[i], ref_seq[i], "tmp1_coverage.fa")

  if os.path.isfile(direc + "/tmp1_coverage.fa"):
    [bait_loc, bait_nam, ref_seq] = minimizeBait("tmp1_coverage.fa", "tmp2_minimal.fa", ref)
    os.system("cat tmp1_minimal.fa  tmp2_minimal.fa >> tmp1_minimalall.fa ")
    print "IIb. Optimizing : to find minimal number of baits covering all genomes...."
    [bait_seq, bait_loc, bait_nam, bait_no, final_bait_fa, ref_seq] = optimalBaits("tmp1_minimalall.fa", ref, "tmp1_optimized.fa", 150)
    [untar_csv, map_csv] = summaryRegionPos(ref_seq, bait_loc, bait_nam, wind, "iter1.fa")
    untar_reg = genomesUntargeted(untar_csv, ref, "iter-2_untargeted_portion.fa")
  else:
    summary("tmp1_minimal.fa", initial_baits_f.replace(".fa", "-iter1.fa"), wind, ref)

  if os.path.isfile("iter-2_untargeted_portion.fa") and os.stat("iter-2_untargeted_portion.fa").st_size != 0:
# Redesign coverage baits for untargetted regions...
    print "III. Filtering & Redesigning Baits for the untargetted regions...."
    ref_seq = {}
    Ref = SeqIO.parse(open(ref, "r"), "fasta")
    for r in Ref:
      ref_seq[r.id] = str(r.seq)
    for i in ref_seq:
      if i in untar_reg:
        print untar_reg[i], ref_seq[i]
        coverageBaits("iter2", untar_reg[i], ref_seq[i], "tmp2_coverage.fa")
    [bait_loc_untar, bait_nam_untar, ref_seq_untar] = minimizeBait("tmp2_coverage.fa", 'tmp3_minimal.fa', "iter-2_untargeted_portion.fa")
    os.system("cat tmp1_optimized.fa  tmp3_minimal.fa >> tmp2_minimalall.fa ")
#    [bait_seq, bait_loc, bait_nam, bait_no, final_bait_fa, ref_seq] = optimalBaits("tmp2_minimalall.fa", ref, "tmp2_optimized.fa", 10)
#    untar_csv = summary( "tmp2_optimized.fa", initial_baits_f.replace(".fa", "-iter3.fa"), wind, ref)
    untar_csv = summary("tmp2_minimalall.fa", initial_baits_f.replace(".fa", "-iter3.fa"), wind, ref)
    untar_reg = genomesUntargeted(untar_csv, ref, "iter-4_untargeted_portion.fa")

  else:
    untar_csv = summary("tmp1_optimized.fa", initial_baits_f.replace(".fa", "-iter2.fa"), wind, ref)
    untar_reg = genomesUntargeted(untar_csv, ref, "iter-3_untargeted_portion.fa")

# Remove unnecessary files
  for todel in glob.glob('*tmp*'): os.remove(todel)
  #for todel in glob.glob('*iter*'): os.remove(todel)

  print "Completed on ", time.strftime("%c")
else:
  print "ERROR: Please enter the arguments as follows..."
  print " CustomizedBaitsPanel.py 'Directory' 'Species Conserved baits fasta file' 'Query reference genomes in fasta'"
