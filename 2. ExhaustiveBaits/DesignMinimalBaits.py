#!/usr/bin/python
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


blast_n = "/usr/bin/blastn"
tab_fmt = "'7 qseqid sseqid qstart qend sstart send pident qcovs evalue bitscore stitle'"
dist = 500  # default 500
bp_len = 120
GC_limit = "20-80"
TM_Min = 68
TM_Max = 170
utr_reg = 250
bp_cov = dist
bp_pos = (bp_cov/2) + (bp_len/2) - 1
bp_neg = (bp_cov/2) - (bp_len/2)
PERC_IND = 80  # default 90


# ------- Thermodynamic oligonucleotide calculator --------
def calGcTm(tmp_seq, loc):
  from Bio.SeqUtils import GC
  from Bio.SeqUtils.MeltingTemp import Tm_staluc
  global TM_Min, TM_Max, bp_len
  [GC_Min, GC_Max] = GC_limit.split('-')
  gc = GC(tmp_seq)
  tm = Tm_staluc(tmp_seq, rna=0)
  thermo = [round(gc, 0), round(tm, 0), tmp_seq]
  if (int(thermo[0]) >= int(GC_Min)) and (int(thermo[0]) <= int(GC_Max)):
    if (int(thermo[1]) >= TM_Min) and (int(thermo[1]) <= TM_Max):
      thermo.append('P')
      return thermo
    else:
      thermo.append('F-Tm')
      return thermo
  else:
    thermo.append('F-GC')
    return thermo


# ------------------------------------- Convert the blast output to the dictionary format--------------
def blastparser(set_file):            # Dictionary from blast op file
  global PERC_IND
  blast_dict = {}
  with open(set_file, 'r') as f:
    for key, group in itertools.groupby(f, lambda line: line.startswith('# BLASTN ')):
      list_Q = []             # initiate tmp list
      if not key:
        for line in group:
          list_Q.append(line.strip('\n'))
        query = list_Q[0]
        for i in range(4, len(list_Q)):
          d = list_Q[i].split('\t')
          if len(d) != 1:
            percent_ind = float(d[6]) * (abs(int(d[5])-int(d[4]))+1)/120         # blast site query length *100/120
            if percent_ind >= PERC_IND:
              if query in blast_dict:
                blast_dict[query].append(d)
              else:
                blast_dict[query] = ([d])
  return blast_dict


# ---------------------------------- Minimize the number of baits (redundant and below 80% similarity)------------------------
def minimizeBait(baitfile, minimalbaitfile, reffile):
  os.system("cdhit  -i " + baitfile + " -o " + minimalbaitfile + " -c 0.95 ")
  [bait_loc, bait_nam, ref_seq] = blastAgainstSeq(reffile, minimalbaitfile)
  return [bait_loc, bait_nam, ref_seq]


# ----------------- Find the optimal number of baits (remove baits overrepresenting a region(500bp)) --------------------
def optimalBaits(file_Sp, reffile, outfile, windowRemove):
  bait_seq = {}
  for i in SeqIO.parse(open(file_Sp, "r"), "fasta"):
    bait_seq[i.id] = str(i.seq)

# group the baits...
  with open(file_Sp.replace(".fa", "-combined.fa"), 'w') as f:
    for i in bait_seq:
      if len(str(bait_seq[i]).replace('N', '')) == 120:
        f.write(">" + i + "\r\n" + str(bait_seq[i]) + "\r\n")
  [bait_loc, bait_nam, ref_seq] = blastAgainstSeq(reffile, file_Sp.replace(".fa", "-combined.fa"))
  group = [(bait_nam[x], x) for x in bait_nam]
  gp = {}
  gpf = {}

  # remove baits if in 200bp window
  for g in range(0, len(group)):
    y = list(set(group[g][0]))
    y.sort()
    for i in range(0, len(y)-1):
      if y[i+1][0] - y[i][0] <= windowRemove:
        # print y[i+1], y[i]
        pair = str(y[i+1][1]) + ":" + str(y[i][1])
        if pair not in gp:
          gp[pair] = [group[g][1]]
        else:
          gp[pair].append(group[g][1])
        if pair not in gpf:
          gpf[pair] = [(y[i+1][0], y[i][0], y[i+1][0] - y[i][0])]
        else:
          gpf[pair].append((y[i+1][0], y[i][0], y[i+1][0] - y[i][0]))
  bait_no = {}
  for p in range(0, len(group)):
    for x in group[p][0]:
      if len(x) == 4:
        if x[1] not in bait_no:
          bait_no[x[1]] = [group[p][-1]]
        elif x[1] in bait_no:
          bait_no[x[1]].append(group[p][-1])
  diff = {}
  for j in gp:
    [s1, s2] = j.split(":")
    diff_s1 = [x for x in set(bait_no[s1]) if x not in set(bait_no[s2])]
    diff_s2 = [x for x in set(bait_no[s2]) if x not in set(bait_no[s1])]
    if s1 not in diff:
      diff[s1] = diff_s1
    else:
      diff[s1].extend(diff_s1)
    if s2 not in diff:
      diff[s2] = diff_s2
    else:
      diff[s2].extend(diff_s2)

# Compare and find the minimum baits
  final_baits = bait_no.keys()

#  print final_baits
  for j in gp:
    [s1, s2] = j.split(":")
    if len(diff[s1]) > len(diff[s2]):
      print "Removing: ", s2
      if s2 in final_baits:
        final_baits.remove(s2)
    elif len(diff[s2]) > len(diff[s1]):
      print "Removing: ", s1
      if s1 in final_baits:
        final_baits.remove(s1)
    if len(diff[s2]) == len(diff[s1]):
      print "Removing....", s2
      if s2 in final_baits:
        final_baits.remove(s2)
  print final_baits

# Select the final baits
  final_bait_fa = {}
  for x in final_baits:
    final_bait_fa[x] = bait_seq[x]
  with open(outfile, 'w') as f:
    for x in final_bait_fa:
      f.write(">" + x + "\r\n" + str(final_bait_fa[x]) + "\r\n")
  [bait_locF, bait_namF, ref_seq] = blastAgainstSeq(reffile, outfile)
  return bait_seq, bait_locF, bait_namF, bait_no, final_bait_fa, ref_seq


# --------------------------------------------- Blast the bait seq to all the strain respresentatives.
def blastAgainstSeq(reffile, baitfile):
  blast_res = {}
  ref_seq = {}
  print "Checking through blast"
  Ref = SeqIO.parse(open(reffile, "r"), "fasta")
  for r in Ref:
    ref_seq[r.id] = str(r.seq)
    refseq_fil = baitfile.replace(".fa", "-tmpref")
    with open(refseq_fil, "w") as f:
      f.write(">" + r.id + "\r\n" + str(r.seq) + "\r\n")
    blast_uni = baitfile.replace(".fa", "-blast")
    blast_run_uni = NcbiblastnCommandline(cmd=blast_n, task='blastn', dust='no', num_alignments=1, evalue=1.0, query=baitfile, subject=refseq_fil, outfmt=tab_fmt, out=blast_uni)
    stdout, stderr = blast_run_uni()
    blast_res[r.id] = blastparser(blast_uni)

  # print blast_res
  blast_list = []
  with open("tmp_blast.csv", "w") as f:
    for x in blast_res:
      if len(blast_res[x]) > 0:
        for i in blast_res[x]:
          for j in range(0, len(blast_res[x][i])):
            n = x + "," + i + "," + str(blast_res[x][i][j]).replace("[",  "").replace("]", "") + "\r\n"
            gg = blast_res[x][i][j]
            gg.insert(0, x)
            gg.insert(1, i)
            blast_list.append(gg)
            f.write(n)

# Get the positions and indentity score each bait against each strain
# bait_nam = {name :[position, bait_name, indentity, coverage ]
  bait_loc = {}
  bait_nam = {}
  for x in blast_list:
    if x[0] not in bait_loc:
      if int(x[6]) > int(x[7]):
        bait_loc[x[0]] = [int(x[7])]
        bait_nam[x[0]] = [(int(x[7]), x[2], float(x[8]), float(x[9]))]
    else:
      if int(x[6]) > int(x[7]):
        bait_loc[x[0]].append(int(x[7]))
        bait_nam[x[0]].append((int(x[7]), x[2], float(x[8]), float(x[9])))
  # print bait_nam
  return bait_loc, bait_nam, ref_seq


# ------------ Write csv for Graph generation -----------
def writeCsv(data, file_nam):
    import csv
    import itertools
    with open(file_nam, 'wb') as f:
        writer = csv.writer(f, delimiter=",")
        print "Writing to :", file_nam
        writer.writerows(list(itertools.chain(t)) for t in data)


# -------------- Create a fasta file containing all the untargetted region------
def genomesUntargeted(infile, reffile, untarfile):
  import csv
  import uuid
  loc_untarget = []
  ref_seq = {}
  with open(infile, 'r') as f:
    for l in csv.reader(f):
       loc_untarget.append(l)
  Ref = SeqIO.parse(open(reffile, "r"), "fasta")
  for r in Ref:
    ref_seq[r.id] = str(r.seq)
  region = {}
  for untarg in loc_untarget:
    if len(untarg) > 1:
      print untarg
      ref = untarg[0]
      for j in range(1, len(untarg)):
        if abs(int(untarg[j].split(" : ")[1].split("-")[0]) - int(untarg[j].split(" : ")[1].split("-")[1])) > 250:
          if ref not in region:
            region[ref] = [untarg[j].split(" : ")]
          else:
            region[ref].append(untarg[j].split(" : "))
  Untar_region = {}
  print "Region :",  region
  for x in region:
    for i in range(0, len(region[x])):
      if int(region[x][i][1].split("-")[0]) < int(region[x][i][1].split("-")[1]):
        Un_seq = ref_seq[x][int(region[x][i][1].split("-")[0]): int(region[x][i][1].split("-")[1])]
      else:
        Un_seq = ref_seq[x][int(region[x][i][1].split("-")[1]): int(region[x][i][1].split("-")[0])]
      if x not in Untar_region:
        Untar_region[x] = [(int(region[x][i][1].split("-")[0]), int(region[x][i][1].split("-")[1]), abs(int(region[x][i][1].split("-")[0]) - int(region[x][i][1].split("-")[1])), Un_seq)]
      else:
        Untar_region[x].append((int(region[x][i][1].split("-")[0]), int(region[x][i][1].split("-")[1]), abs(int(region[x][i][1].split("-")[0]) - int(region[x][i][1].split("-")[1])), Un_seq))
  with open(untarfile, 'w') as f:
    for i in Untar_region:
      for j in range(0, len(Untar_region[i])):
        if "N" in list(Untar_region[i][j][3]):
          print Untar_region[i][j][3].find('N')
          if len(max(Untar_region[i][j][3].split("N"), key=len)) > 250:
            f.write(">" + str(uuid.uuid4()) + i + "|" + str(Untar_region[i][j][0]) + "-" + str(Untar_region[i][j][1]) + "|" + str(Untar_region[i][j][2]) + "\r\n" + str(Untar_region[i][j][-1]) + "\r\n")
        else:
          f.write(">" + str(uuid.uuid4()) + i + "|" + str(Untar_region[i][j][0]) + "-" + str(Untar_region[i][j][1]) + "|" + str(Untar_region[i][j][2]) + "\r\n" + str(Untar_region[i][j][-1]) + "\r\n")
  # print "******** ", Untar_region
  return Untar_region


# ----------- Design coverage baits -----------------------------
def coverageBaits(name, uncovered_target_reg, Ref_Seq, cov_file):
  import uuid
  import os
  from Bio.Seq import Seq
  [GC_Min, GC_Max] = GC_limit.split('-')
  no_bait = []
  cbaits = []
  refSeq_cov = str(Ref_Seq)
  for reg in uncovered_target_reg:
    df = reg[2]/float(bp_cov)
    # print "DF: ", df
    if df <= 1.1:
      initial_start = reg[0] + int(reg[2]/2 - 60)         # select a kmer in the middle of the range
      # print "init start" , reg, initial_start
      [k, start_kmer, kmer_s] = findOptimalcoverage(reg, refSeq_cov, initial_start)
      if len(kmer_s.replace("N", '')) == 120:
        k = calGcTm(kmer_s, start_kmer)
        if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
          cbaits.append((start_kmer, name, kmer_s, k[0], k[1]))
          with open(cov_file, 'a') as f:
            f.write('>' + name + str(uuid.uuid4()) + '|' + str(len(Ref_Seq)) + '|' + str(start_kmer) + "\r\n" + str(kmer_s) + "\r\n")
        else:
          print "failed kmer:", k
          no_bait.append(reg)
      else:
        print "failed kmer:", k
        no_bait.append(reg)
    else:                               # multiple baits with a coverage of 500bp each
      reps = round(round((reg[2]/float(bp_cov)), 2))
      for i in range(0, int(reps)):
        initial_start = bp_neg + reg[0] + bp_cov * i + i
        [k, start_kmer, kmer_s] = findOptimalcoverage(reg, refSeq_cov, initial_start)
        if len(kmer_s.replace("N", '')) == 120:
          k = calGcTm(kmer_s, start_kmer)
          if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
            cbaits.append((start_kmer, name, kmer_s, k[0], k[1]))
            with open(cov_file, 'a') as f:
              f.write('>' + name + str(uuid.uuid4()) + '|' + str(len(Ref_Seq)) + '|' + str(start_kmer) + "\r\n" + str(kmer_s) + "\r\n")
          else:
            no_bait.append(reg)
        else:
          no_bait.append(reg)
  print cbaits
  print "No bait", no_bait
  return(cbaits, no_bait)


# move the window from the bait position to 60bp each side to design a bait...
def findOptimalcoverage(reg, refSeq_cov, initial_start):
  from Bio.Seq import Seq
 # initial_start = reg[0] + int(reg[2]/2 - 60)    # select a kmer in the middle of the range
  [GC_Min, GC_Max] = GC_limit.split('-')
  max_len = len(refSeq_cov)
  kmer_s = ''
  k = [0, 0, '']
  for x in range(0, 60):
    for j in ['+', '-']:
      if j == '+':
        start_kmer = initial_start + x
      if j == '-':
        start_kmer = initial_start - x
      # if start_kmer + bp_len < max_len:
      kmer_s = str((Seq(refSeq_cov[start_kmer - 1: (start_kmer - 1) + bp_len]).reverse_complement()).upper())          # 120 cDNA of target in upper case
      # else:
        # kmer_s = str((Seq(refSeq_cov[max_len - 130 : max_len - 10]).reverse_complement()).upper())
      if len(kmer_s.replace("N", '')) == 120:
        k = calGcTm(kmer_s, start_kmer)
        if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
          break
    if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
      break
  print "Kmerss....", kmer_s, reg,  start_kmer
  return k, start_kmer, kmer_s


# ------------ Find untargeted region -----------------------
def untargetRegion(index_loc, dist, Ref_Seq):
    global bp_cov, bp_len, bp_pos, bp_neg
    uncovered_target_reg = []
    len_refseq = len(Ref_Seq[0])
    if len(index_loc) > 1:
        for i in range(0, len(index_loc)):
            if i == 0:          # Start of the genome
                diff = index_loc[i] - 1        # region before the first bait
                if diff > utr_reg:
                    uncovered_target_reg.append((1, index_loc[i] - bp_neg - 1, (index_loc[i] - bp_neg - 1) - 1))
                diff = index_loc[i + 1] - index_loc[i]    # region between 1st n 2nd bait
                if diff > dist + bp_cov:
                    if i == len(index_loc)-1:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, len_refseq, (len_refseq) - (index_loc[i] + bp_pos + 1)))
                    else:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, index_loc[i + 1] - bp_neg - 1, (index_loc[i + 1] - bp_neg - 1) - (index_loc[i] + bp_pos + 1)))
            elif i == len(index_loc)-1:         # End of the Genome
                diff = len_refseq - index_loc[i]
                if diff > utr_reg:
                    if i == len(index_loc)-1:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, len_refseq, (len_refseq) - (index_loc[i] + bp_pos + 1)))
                    else:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, index_loc[i + 1] - bp_neg - 1, (index_loc[i + 1] - bp_neg - 1) - (index_loc[i] + bp_pos + 1)))
            else:
                diff = index_loc[i + 1] - index_loc[i]  # In middle
                if diff > dist + bp_cov:
                    if i == len(index_loc)-1:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, len_refseq, (len_refseq) - (index_loc[i] + bp_pos + 1)))
                    else:
                      uncovered_target_reg.append((index_loc[i] + bp_pos + 1, index_loc[i + 1] - bp_neg - 1, (index_loc[i + 1] - bp_neg - 1) - (index_loc[i] + bp_pos + 1)))
    else:
        diff = index_loc[0] - 1
        if diff > dist + bp_cov:
            uncovered_target_reg.append((1, index_loc[0] - bp_neg - 1, (index_loc[0] - bp_neg - 1) - 1))
        diff = len_refseq - index_loc[0]
        if diff > dist + bp_cov:
            uncovered_target_reg.append((index_loc[0]+bp_pos+1, len_refseq, (len_refseq) - (index_loc[0]+bp_pos+1)))
    print "UNTARGET REGIoN", uncovered_target_reg
    return uncovered_target_reg


# -------------  Summary files
def summaryRegionPos(ref_seq, bait_loc, bait_namF, untar_window, out_filnam):
  baits_mapped = []
  loc_untar = []
  for x in ref_seq:
    if x in bait_namF:
      bait_loc[x].sort()
      loc = [i[0] for i in bait_namF[x]]
      bait_pos = [i[1] + " : " + str(i[0]) for i in bait_namF[x]]
      bait_pos.insert(0, x)
      baits_mapped.append(bait_pos)
      untar = [x]
      k = bait_loc[x]
      k.sort()
      lll = untargetRegion(k, dist, ref_seq[x])
      print lll
      for l in lll:
        untar.append(str(l[1]-l[0]) + " : " + str(l[0]) + "-" + str(l[1]))
   # Untargetted region is dist + 120 = 500 + 120 =620
#      for o in range(0, len(loc)-1):
#        if loc[o+1] - loc[o] > untar_window:
#          untar.append(str(loc[o+1]-loc[o]) + " : " + str(loc[o]) + "-" + str(loc[o+1]))
      loc_untar.append(untar)  #
    else:
      loc_untar.append([x, "Full : 1-" + str(len(ref_seq[x]))])
      baits_mapped.append([x])
  writeCsv(baits_mapped, out_filnam.replace(".fa", "-mapped-position.csv"))
  writeCsv(loc_untar, out_filnam.replace(".fa", "-mapped-untargeted.csv"))

  mapp_reformat = []
  bait_name = {}
  for m in baits_mapped:
    mapp_re = [m[0], len(ref_seq[m[0]])]
    for n in range(1, len(m)):
      c = m[n].split(":")
      d = [int(c[1]) - bp_neg, int(c[1]) + bp_pos]
      if d[0] < 1:
        d[0] = 1
      if d[1] > len(ref_seq[m[0]]):
        d[1] = len(ref_seq[m[0]])
      mapp_re.append(c[0])
      mapp_re.extend(d)
      if c[0] not in bait_name:
        bait_name[c[0]] = [m[0]]
      else:
        bait_name[c[0]].append(m[0])
    mapp_reformat.append(mapp_re)
  writeCsv(mapp_reformat, out_filnam.replace(".fa", "-mapped-region.csv"))
  writeCsv([[e, bait_name[e]] for e in bait_name], out_filnam.replace(".fa", "-baits_ref.csv"))
  writeCsv([[x, len(bait_name[x])] for x in bait_name], out_filnam.replace(".fa", "-baits_ref_count.csv"))
  return out_filnam.replace(".fa", "-mapped-untargeted.csv"),  out_filnam.replace(".fa", "-mapped-position.csv")
