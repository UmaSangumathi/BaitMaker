#!/usr/bin/python
##############################################################
#  Script     : EnrichmentAnalysis_NotSegmented.py 
#  Author     : uma sangumathi
#  Date       : 14/05/2015
#  Last Edited: 23/06/2015, uma 
#  Description: Consensus genome and snps 
##############################################################
# Purpose: 
#  Create guided consensus genome from the reads and nearest reference genome for segmented viruses    
#  Map reads using bowtie2    
#  Snps detected using Lofreq2     
# Requirements:
#  1. Biopython module    
#  2. Bwa
#  3. Bowtie2
#  4. components of Vipr pipeline  
#  5. Samtools 
#  6. NCBI blast commandline tool and database
#############################################################

import os
import sys
import glob
import datetime
import csv
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Entrez


def group_nts(data):
  from operator import itemgetter
  from itertools import groupby
  ranges = []
  for k, g in groupby(enumerate(data), lambda (i,x):i-x):
    group = map(itemgetter(1), g)
    ranges.append((group[0], group[-1]))
  return ranges


# Pairend illumina 250bp sequences
if len(sys.argv)==9:
  f1 = sys.argv[1]
  Entrez.email = sys.argv[2]
  fq_dir = sys.argv[3]
  ref_dir = sys.argv[4]
  ncbi_db = sys.argv[5]
  summary_file = sys.argv[6]
  out_dir = sys.argv[7]
  cons_code = sys.argv[8]
  tab_fmt = "'6 qseqid sseqid qstart qend sstart send pident qcovs evalue bitscore stitle'"


# Read the summary file:
  info_dict = {}
  with open(summary_file, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for xi in reader:
        print xi
    info_dict[xi[0]] = xi[1]
  print info_dict.keys()

  for i in info_dict: 
    if i in f1 :
      ref_id =  info_dict[i]
  ref = ref_dir + "/" + ref_id 
  threads = str(3)
  op_prefix = f1.replace(".fq", "-" + ref_id.split(".")[0])

  f2 = f1.replace("R1", "R2").replace("val_1","val_2") 
  if not os.path.exists(out_dir): os.mkdir(out_dir)
  res_dir = out_dir + "/"+ f1.replace(".fq", "")
  if not os.path.exists(res_dir): os.mkdir(res_dir)

 # Log file
  old_stdout = sys.stdout
  log_file = open(res_dir + "/" + f1.replace(".fq", ".log"),"a")
  sys.stdout = log_file
  old_stderr = sys.stderr
  err_file = open(res_dir + "/" + f1.replace(".fq", ".err"),"a")
  sys.stderr = err_file   

  print "Processing....", f1
  print str(datetime.datetime.now())  
  #1- Consensus: Map reads to a reference genome...
  cons_ref = res_dir + "/" + op_prefix + "-1.fa"  
  cons_cmd = cons_code + "   -f " + f1 + " -g " + f2 + " --illumina  -r " + ref + " -t " + threads + " --force -o " + cons_ref    
  print cons_cmd
  if not os.path.isfile(cons_ref): os.system(cons_cmd) 
     
  # Blast the consensus formed.... 
  b_out = res_dir + "/" + f1.replace(".fq", ".tmp")
  blast_run_uni = NcbiblastnCommandline(cmd="blastn", task="megablast", db=ncbi_db, max_target_seqs=1, query=cons_ref , outfmt=tab_fmt, out=b_out)
  print "Runnin blast...", blast_run_uni
  nearest_ref1 = cons_ref.replace(".fa", "-nref1.fa")
  if not os.path.isfile(nearest_ref1): 
    stdout, stderr = blast_run_uni()
    c = 0
    for line in open(b_out, "r"):
      c = c+1
      if c ==1 : u_id=line.split("|")
    print u_id[1] 
    handle = Entrez.efetch(db="nuccore", id=u_id[1], retmode='xml', rettype="gb") 
    record = Entrez.read(handle, validate=False)
    with open( nearest_ref1, "w") as f :
      f.write(">" + record[0]["GBSeq_primary-accession"] + "\r\n" + record[0]["GBSeq_sequence"] + "\r\n")
   
  # 2-iteration Map reads to nearest reference...
  cons_ref2 = res_dir + "/" + op_prefix  +  "-cons-2.fa" 
  cons_cmd2 = cons_code + "   -f " + f1 + " -g " + f2 + " --illumina  -r " + nearest_ref1 + " -t " + threads + " --force -o " + cons_ref2    
  print cons_cmd2
  if not os.path.isfile(cons_ref2):os.system(cons_cmd2)
    
  
# get Hybrid reference genome....
  print "Getin hybrid consensus genome..."
  hybrid_reffile =  cons_ref2.replace(".fa", "-hybrid.fa")
  if not os.path.isfile(cons_ref2.replace(".fa", "-comb.fa")):os.system(" cat " + nearest_ref1 + " " + cons_ref2  + " > " +  cons_ref2.replace(".fa", "-comb.fa"))
  mafft_cmd = "mafft  --globalpair --maxiterate 16 --inputorder " + cons_ref2.replace(".fa", "-comb.fa") + "  > " + cons_ref2.replace(".fa", "-comb.align")
  os.system(mafft_cmd)
  X = SeqIO.parse(cons_ref2.replace(".fa", "-comb.align"), "fasta") 
  ref_align = []
  for rec in X: 
    ref_align.append([rec.id,  str(rec.seq)])
  hybrid_ref = []
  for i in range(0, len(ref_align[1][1])):
    if ref_align[1][1][i] == "n" or ref_align[1][1][i] == "-":
      N = ref_align[0][1][i].lower() 
    else: N = ref_align[1][1][i].upper()
    hybrid_ref.append(N) 
  from_reference = [nt + 1  for nt in range(0,len(hybrid_ref)) if hybrid_ref[nt].islower()] 
  from_consenses = [nt + 1  for nt in range(0,len(hybrid_ref)) if hybrid_ref[nt].isupper()] 
  from_reference_gp = group_nts(from_reference)
  from_consenses_gp = group_nts(from_consenses) 
  with open(hybrid_reffile.replace(".fa", ".annot"), "w") as f:
    for h in from_reference_gp:
      f.write("From reference used: " + "\t" + str(h[0]) + "-" + str(h[1]) + "\r\n")
    for h in from_consenses_gp:
      f.write("From consensus generated: " + "\t" + str(h[0]) + "-" + str(h[1]) + "\r\n")
  new_ref = "".join(hybrid_ref)  
  with open( hybrid_reffile , "w") as f: 
    f.write(">Hybrid_consensus" + "\r\n" + new_ref + "\r\n")  
  

# Map reads to consensus genome.... Bowtie2
  hyb_cons_indx = "bowtie2-build " + hybrid_reffile  + " " + res_dir + "/tmp_idx" 
  os.system(hyb_cons_indx)
  out_sam = hybrid_reffile.replace(".fa", ".sam")
  bowtie2_cmd = "bowtie2  -x " + res_dir + "/tmp_idx  -1 " + f1 + " -2 " + f2  + " -S  " + out_sam + " -p 8"
  os.system(bowtie2_cmd)
  print bowtie2_cmd
  # run lofreq2....
  bam_cmd = "samtools view -b -S " + out_sam  + " > " + out_sam.replace(".sam", ".bam")
  bam_sort = "samtools sort " + out_sam.replace(".sam", ".bam") + " " + out_sam.replace(".sam", "-sort")
  cmd = "lofreq   call  -f  " + hybrid_reffile  + " -o " + out_sam.replace(".sam", "-snps.vcf") + " " + out_sam.replace(".sam", "-sort.bam")
  print bam_cmd
  os.system(bam_cmd)
  print bam_sort
  os.system(bam_sort)
  print cmd
  os.system(cmd)
  os.remove(out_sam.replace(".sam", ".bam")) 
  print "\n\n"
  os.system("coverage_plot.py  -b " + out_sam.replace(".sam", "-sort.bam")  + " -o " + out_sam.replace(".sam", ".pdf") + " ")
 # genome coverage at particular position
  os.system("bedtools genomecov -g " + hybrid_reffile + "  -ibam  " + out_sam.replace(".sam", "-sort.bam") + " -d  > " + out_sam.replace(".sam", ".cov")) 

# close logs
  sys.stdout = old_stdout
  log_file.close() 
  sys.stderr = old_stderr
  err_file.close()

else:
  print "Error! please enter the the correct arguments"

