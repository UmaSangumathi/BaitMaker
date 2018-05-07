#!/usr/bin/python

import os
import sys
import glob
import datetime
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Entrez
from math import log10
import fileinput
from bait_effiency_modules import *
Entrez.email = "umasang7@gmail.com"
blast_n = "/usr/bin/blastn"
tab_fmt = "'6 qseqid sseqid qstart qend sstart send pident qcovs evalue bitscore stitle'"

Genome_Strand = "Positive"

def write_csv(data, file_nam):
    import csv, itertools
    with open(file_nam, 'wb') as f:
        writer = csv.writer(f, delimiter="\t")
        print "Witing to :", file_nam
        writer.writerows(list(itertools.chain(t)) for t in data)


# Input files, .cov, .vcf, .ref, .baits

if len(sys.argv)==8:
  cov_file = sys.argv[1]
  vcf_file = sys.argv[2]
  ref_file = sys.argv[3] 
  baits_file = sys.argv[4]
  genes_file = sys.argv[5]
  coords_file = sys.argv[6]
  lib_typ = sys.argv[7]
  
# 1. process the coverage file....
  with open(cov_file , "r") as f:
    cov_list = f.readlines() 
  log_cov = []
  for cc in cov_list:
    c = cc.split("\t") 
    if int(c[2].strip()) != 0: log_v = log10(int(c[2].strip()))
    else: log_v = 0 
    log_cov.append([c[0], c[1], c[1], log_v])
  write_csv(log_cov, cov_file.replace(".cov", ".logcov"))
  ref_id = c[0]

# 2. process the snp file.... 
  with open(vcf_file, "r") as f:
    vcf_list = f.readlines()
  snp_plot = []
  for vv in vcf_list:
    if not vv.startswith("#"):
      v = vv.split("\t")
      snp_plot.append([v[0], v[1], v[1], v[-1].split(";")[1].split("=")[1]])         
  write_csv(snp_plot, vcf_file.replace(".vcf", ".snps")) 

# 3. Process the baits file.... 
  if lib_typ == "Capture": 
    baits = {}
    c=0
    for b in SeqIO.parse(baits_file, "fasta"):
      c = c + 1
      baits[b.id] = [c, b.id , str(b.seq)]
  # blast baits and the consensus reference genome (blastn used to get any possible anchor positions of the baits on the genome)
    b_out = cov_file.replace(".cov", ".blastbaits")
    blast_run_uni = NcbiblastnCommandline(cmd=blast_n, task="blastn", subject=ref_file, max_target_seqs=1, query=baits_file , outfmt=tab_fmt, out=b_out)
    print "Runnin blast...", blast_run_uni
    stdout, stderr = blast_run_uni()
  # parse the blast result
    with open (cov_file.replace(".cov", ".blastbaits")) as f:
      hits = f.readlines()
    Bait_map_sort = proccess_blasthits(Genome_Strand, hits, baits)
    write_csv(Bait_map_sort,  cov_file.replace(".cov", "-baits_mapped.csv"))
    bait_pos = []
    for b in Bait_map_sort[1:len(Bait_map_sort)]:
      if len(b) > 5 and (b[-1] != "Fail" and b[-1] != "NA") : 
        if b[-1] == "High" : j = 1
        else: j = 0.75
        bait_pos.append([ref_id, b[-3], b[-2], j])
    write_csv(bait_pos, cov_file.replace(".cov", ".baitspos"))    
  else: 
    write_csv([],  cov_file.replace(".cov", "-baits_mapped.csv"))
    write_csv([], cov_file.replace(".cov", ".baitspos"))
# 4. Get the Genes names....
  os.system("cp " + genes_file + "  " + cov_file.replace(".cov", ".genes"))  
  os.system("cp " + coords_file + "  " + cov_file.replace(".cov", ".coords"))

# 5. Circos plot...
  with open("tmp.conf", "w") as f:
    for line in fileinput.input("hist.conf", inplace=False):
      f.write(line.replace("XXXXX", cov_file.replace(".cov", "")))
  circos_cmd = "circos  -conf  tmp.conf " 
  os.system(circos_cmd)        



