#!/usr/bin/python

# ----------- Global default values -------------

TM_Min = 68
TM_Max = 150
# bp_len = 120
# bp_cov = 500  # default 500
# bp_pos = 309  # (bp_cov/2)+59  # default 309
# bp_neg = 190  # (bp_cov/2)-60  # default 190


# ------- Thermodynamic oligonucleotide calculator --------
def GC_TM_calc(tmp_seq, loc, GC_comb, disRC):
    from Bio.Seq import Seq
    from Bio.SeqUtils import GC
    from Bio.SeqUtils.MeltingTemp import Tm_staluc
    global TM_Min, TM_Max
    [GC_Min, GC_Max] = GC_comb.split('-')
    if disRC == 'False':    
        gc = GC(tmp_seq)
        tm = Tm_staluc(tmp_seq, rna=0)
        thermo = [round(gc, 2), round(tm, 2), tmp_seq]
    elif disRC == 'True':
        tmp_seqd = str(Seq(tmp_seq).reverse_complement())
        gc = GC(tmp_seqd)
        tm = Tm_staluc(tmp_seqd, rna=0)
        thermo = [round(gc, 2), round(tm, 2), tmp_seq]

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


# ------------------- Covert data in query file to multi-level kmer dictionary ----------------------
def FileToDict(query_file, mismatch):
    import fileinput
    import itertools as it
    from collections import defaultdict
    print "\nMismatch allowed: ", str(mismatch)
    print "\n1. Covert data in query file to multi-level kmer dictionary"
    l = lambda: defaultdict(l)
    kmer_dict = defaultdict(lambda: defaultdict(defaultdict))

    # parse query file which is segmented in blocks : starts '>' and ends '<'
    with open(query_file, 'r') as f:
        for key, group in it.groupby(f, lambda line: line.startswith('>\n')):    # group the block with start identifier ; key = logical(start identifier or not)
            list_Q = []	 # initiate tmp list
            if not key:
                for line in group:
                    list_Q.append(line.strip('\n'))
                q_id = list_Q[0]
                seq = list_Q[2]
                num_h = list_Q[3]
                for i in range(4, len(list_Q) - 1):  # before '<' end indentifier
                    d = list_Q[i].strip(' ').split(' ')
                    vhit_id = d[0]  # virus_id
                    loc = d[2]  # location on virus_id
                    degen_info = d[3:len(d)]  # degeneration info
                    if seq in kmer_dict:  # multi-level dictionary main key: sequence; second key: kmer_id; third key: number of virus hits
                        if q_id in kmer_dict[seq]:
                          kmer_dict[seq][q_id][num_h].append((vhit_id, loc, degen_info))
                        else:
                          kmer_dict[seq][q_id][num_h] = [(vhit_id, loc, degen_info)]
                    else:
                      kmer_dict[seq][q_id][num_h] = [(vhit_id, loc, degen_info)]

    print "No. of distinct kmers in the query file - kmer_dict: ", len(kmer_dict)
    return (kmer_dict)


# ---------------- Unique kmer is used as key, removes redundant kmer position info -----------
def MergeInfo(kmer_dict):
    print "\n2. kmer-sequence is used as key, removes redundant kmer position info"
    tot_hits = {}
    for seq in kmer_dict:
        loc_tup = []
        for q in kmer_dict[seq].keys():
           for n_hit in kmer_dict[seq][q].keys():
                X = kmer_dict[seq][q][n_hit]
                for x in X:
                    if x not in loc_tup:
                      loc_tup.append(x)
        if seq in tot_hits:
          tot_hits[seq].append(loc_tup)
        else:
          tot_hits[seq] = loc_tup

    print "No. of kmer after removing duplicates - tot_hits: ", len(tot_hits)
    return (tot_hits)


# -------------- Discard the sequences with les than __% of the total no. of hits --------------
def HitsCoverage(tot_hits, filter_ratio, virus_index):
    global vir_index
    print "Hits Coverage filter ::::::: ", filter_ratio, "::::::::::"
    print "\n3. Creating NCBI hits index and filter kmer sequence if the hits coverage-"
    with open(virus_index, 'r') as f:
        vir_index = f.readlines()

    kmer_cons = {}
    for i in tot_hits.keys():
        if len(tot_hits[i]) > (len(vir_index) * (filter_ratio/100.0)):
            kmer_cons[i] = tot_hits[i]

    print "Num of Sequences representing the virus: ", len(vir_index)
    print "Filter ratio: ", filter_ratio
    print "Coverage of the hits greater than", len(vir_index) * (filter_ratio/100.0)
    if len(kmer_cons) == 0:
      print "NO KMERS clearing the filter..........TOOOOO High Set DIVERSED"
    else:
      print "Kmers passed hits coverage filter: ", len(kmer_cons)
    return kmer_cons


# --------------- Select a Reference Sequence for mapping the baits ----------------
def RefMapping(kmer_cons, cur1):
    global vir_index, virus_counter, Ref_Seq
    print "\n4. Selecting a Reference Sequence to map the kmers"
    Ref_seqloc = {}
    index_dict = {}
    for i in range(0, len(vir_index)):
        index_dict[i] = vir_index[i].strip("\r\n").split('| ')[1]

    virus_counter = {}      # get the virus_ids count
    for i in kmer_cons.keys():
        for j in range(0, len(vir_index)):
            for n in range(0, len(kmer_cons[i])):
                if int(kmer_cons[i][n][0]) == j:
                    if j in virus_counter.keys():
                        virus_counter[j] = virus_counter[j] + 1
                    else:
                        virus_counter[j] = 1

    max_v_c = 0
    for l in virus_counter:
        if virus_counter[l] > max_v_c:
            max_v_c = virus_counter[l]

    Ref_seq_list = []
    for l in virus_counter:
        if virus_counter[l] > (max_v_c - 5) and virus_counter[l] != 0:
          Ref_seq_list.append(l) 	# Reference sequence selected from eg) Maximum num mapped = 50 then sequence with 45-50 are taken into account

    print "Ref_seq:", Ref_seq_list
    ref_u_id = {}
    for r in Ref_seq_list:
        ref_u_id[index_dict[r]] = r

    target_set = []
    for u_id in ref_u_id:
        ref = cur1.execute("SELECT RefSeq_id from RefSeq_records WHERE Unique_id=?", (u_id,))
        if ref != None:
            for rf in ref:
                print rf
                if rf != "NULL":
                    target_set.append(int(u_id))

    print "Reference sequence found length: ", len(Ref_seq_list)
    if len(target_set) == 0:
      target_set.append(index_dict[Ref_seq_list[0]])

    print "target set ", target_set
    refseq = cur1.execute("SELECT Sequence from NCBI_DNA_records  WHERE Unique_id=?", (target_set[0],))
    Ref_Sequence = [t[0] for t in refseq]
    Ref_Seq = [str(a) for a in Ref_Sequence]
    virus_id = str(ref_u_id[str(target_set[0])])
    for i in kmer_cons.keys():
        loc_tup = [ti for ti in kmer_cons[i] if ti[0] == str(virus_id)]
        if len(loc_tup) != 0:
            loc = [ti[1] for ti in loc_tup]
            if loc[0] not in Ref_seqloc:
              Ref_seqloc[loc[0]] = ([[i, loc_tup]])
            else:
              Ref_seqloc[loc[0]].append([i, loc_tup])
    print "The Reference genome : ", target_set[0]
    print "Kmers mapped to reference genome: ", len(Ref_seqloc)
    return (Ref_seqloc, target_set[0], ref_u_id, virus_id, Ref_Seq)


# Map the kmers to known reference Seq
def Map_knownRef(kmer_cons, virus_id):
    Ref_seqloc = {}
    for i in kmer_cons.keys():
        loc_tup = [ti for ti in kmer_cons[i] if ti[0] == str(virus_id)]
        if len(loc_tup) != 0:
            loc = [ti[1] for ti in loc_tup]
            if loc[0] not in Ref_seqloc:
              Ref_seqloc[loc[0]] = ([[i, loc_tup]])
            else:
              Ref_seqloc[loc[0]].append([i, loc_tup])
    print "\n4. Kmers mapped to reference genome ", virus_id, ' : ', len(Ref_seqloc)
    return (Ref_seqloc)


def GetReference(vir, cur1, virus_dir):
    global Ref_Seq
    global vir_dir, virus_nam
    vir_dir = virus_dir
    virus_nam = vir
    ref_id = None
    ref = cur1.execute("SELECT Unique_ID, RefSeq_id from RefSeq_records WHERE Species=?", (vir,))
    for rf in ref:
        id = [t for t in rf]
        if id[1] != "NULL":
            refseq = cur1.execute("SELECT Sequence from NCBI_DNA_records WHERE Unique_id=?", (id[0], ))
            Ref_Sequence = [t[0] for t in refseq]
            Ref_Seq = [str(a) for a in Ref_Sequence]
            ref_id = id[0]

    if ref_id == None:
        refseq = cur1.execute("SELECT Unique_id, Sequence from NCBI_DNA_records WHERE Species=? ORDER BY Sequence_length DESC limit 1", (vir, ))
        refer = [t for t in refseq]
        ref_id = refer[0][0]
        Ref_Seq = [refer[0][1]]
    print Ref_Seq
    return ref_id, Ref_Seq


# ------ Select the least degenerate kmer on target seq ----------
def LeastDegenerate(Ref_seqloc, mismatch):
    global seq_loc
    print "\n5. Select the least degenerate kmer mapped on target seq if multiple kmers present on same location"
    kmer_seq = []
    seq_loc = {}
    for l in Ref_seqloc.keys():
        min_degen = int(mismatch)
        for n in range(0, len(Ref_seqloc[l])):
            if len(Ref_seqloc[l][n]) != 0:
                tmp = int(Ref_seqloc[l][n][1][0][2][0])
                if tmp < min_degen:
                  min_degen = tmp

        for n in range(0, len(Ref_seqloc[l])):
            if (len(Ref_seqloc[l][n]) != 0) and (Ref_seqloc[l][n][1][0][2][0] == str(min_degen)):
                kmer_seq.append((Ref_seqloc[l][n][1][0][1], Ref_seqloc[l][n][0], min_degen))
                if int(Ref_seqloc[l][n][1][0][1]) not in seq_loc:
                  seq_loc[int(Ref_seqloc[l][n][1][0][1])] = ([Ref_seqloc[l][n][0]])
                else:
                  seq_loc[int(Ref_seqloc[l][n][1][0][1])].append([Ref_seqloc[l][n][0]])

    sort_loc = sorted(kmer_seq, key=lambda tup: int(tup[0]))  # Sort based on location of kmer
    uniq_loc = set([t[0] for t in sort_loc])
    sort_u_loc = sorted(uniq_loc)
    return (sort_u_loc)


# ------------- 120mers & loc_number ----------
def seqToloc(sort_u_loc):
    global seq_loc
    comb_seq_loc = {}
    for i in sort_u_loc:
        i = int(i)
        comb_120s = seq_loc[i][0]
        comb_seq_loc[i] = comb_120s
    return(comb_seq_loc)


# ------------- GC and TM filter ------------------
def testGC_Tm(comb_seq_loc, bait_dir, virus_nam, key, GC, disRC):
    from Bio.Seq import Seq
    import os
    global max_loc
    print "\n6. Find GC and TM of the kmer"
    passed_kmer = {}
#    filtered_file = bait_dir + "/" + virus_nam + "_passed_kmers.fa"
#   if os.path.exists(filtered_file): os.remove(filtered_file)
    max_loc = 0
    for loc in comb_seq_loc:
        k = GC_TM_calc(comb_seq_loc[loc], int(loc), GC, disRC)
        if k[3] == 'P':
          passed_kmer[loc] = (comb_seq_loc[loc], k[0], k[1])
        if loc > max_loc:
          max_loc = loc
    # for l in passed_kmer:
        # with open(filtered_file, 'a') as f:
           # f.write('> ' + str(l) + ' | ' + key  + "\r\n" + str(passed_kmer[l]) + "\r\n")
    return (passed_kmer)


# ------------ Write csv for Graph generation -----------
def write_csv(data, file_nam):
    import csv
    import itertools
    with open(file_nam, 'wb') as f:
        writer = csv.writer(f, delimiter=",")
        print "Writing to :", file_nam
        writer.writerows(list(itertools.chain(t)) for t in data)


# -------------- Remove Overlapping kmers ----------------------
def RemoveOverlapBaits(passed_kmer, key, virus_dir, vir_nam, uniq_id, blast_infile, bp_cov, disRC):
    import os
    from Bio.Seq import Seq
    global max_loc, baits, Sp_cons_bait  # , bp_pos, bp_neg
    global vir_dir, virus_nam
    vir_dir = virus_dir
    virus_nam = vir_nam

    print "\n7. Discarding overlapping kmers... "
    Sp_cons_bait = {}
    prob_cov = 0
    index_loc = []
    baits = []

    Sp_cons_file = virus_dir + "/" + virus_nam + "--bait.fa"
    print "Max_loc: ", max_loc
    for i in range(1, max_loc+1):
        if i in passed_kmer:
            if i > 50 and i > prob_cov:
                print 'Overlap : passed write to file', i
                Sp_cons_bait[i] = passed_kmer[i]
                index_loc.append(i)
                prob_cov = i + bp_cov
                if disRC == 'False':
                    baits.append((i + 1, "Sp_cons-" + str(key), passed_kmer[i][0], passed_kmer[i][1], passed_kmer[i][2])) 	# location i+1 indicates the real location on the genome
                    with open(Sp_cons_file, 'a') as f:
                        f.write('>  Sp_cons-bait |' + str(key) + ' | ' + str(i+1) + ' | ' + uniq_id + "\r\n" + str(Sp_cons_bait[i][0]) + "\r\n")
                    with open(blast_infile, 'a') as f:
                        f.write('>  Sp_cons-targ |' + str(key) + ' | ' + str(i+1) + ' | ' + uniq_id + "\r\n" + str(Seq(Sp_cons_bait[i][0]).reverse_complement()) + "\r\n")
                elif disRC == 'True':
                    baits.append((i + 1, "Sp_cons-" + str(key), str(Seq(passed_kmer[i][0]).reverse_complement()), passed_kmer[i][1], passed_kmer[i][2])) 	# location i+1 indicates the real location on the genome
                    with open(Sp_cons_file, 'a') as f:
                        f.write('>  Sp_cons-bait |' + str(key) + ' | ' + str(i+1) + ' | ' + uniq_id + "\r\n" + str(Seq(Sp_cons_bait[i][0]).reverse_complement()) + "\r\n")
                    with open(blast_infile, 'a') as f:
                        f.write('>  Sp_cons-targ |' + str(key) + ' | ' + str(i+1) + ' | ' + uniq_id + "\r\n" + str(Sp_cons_bait[i][0]) + "\r\n")    

    print "Num of Non Overlapping baits: ", len(baits)
    return (index_loc, baits)


def RmvOverUntar(passed_kmer, uncovered_target_reg, key, uniq_id, dist, blast_infile, bp_len, bp_cov, disRC):
    # global bp_neg, bp_pos
    global vir_dir, virus_nam
    new_baits = {}
    for kmerloc in passed_kmer:
        for uncovReg in uncovered_target_reg:
            if kmerloc > (uncovReg[0] + dist - bp_len) and kmerloc < (uncovReg[1] - dist + bp_len):
                print "New bait found .... ", kmerloc, passed_kmer[kmerloc]
                new_baits[kmerloc] = passed_kmer[kmerloc]
    print new_baits
    [n_index_loc, n_baits] = RemoveOverlapBaits(new_baits, key, vir_dir, virus_nam, uniq_id, blast_infile, disRC)
    return (n_index_loc, n_baits)


# ------------ Find untargeted region -----------------------
def UntargetRegion(index_loc, dist, bp_len, bp_cov, bp_pos, bp_neg):
    global Ref_Seq  # bp_len, bp_pos, bp_neg, 

    print "\n8. Finding the untargeted region... start, end, difference"
    uncovered_target_reg = []
    len_refseq = len(Ref_Seq[0])
    print index_loc
    if len(index_loc) > 1:
        for i in range(0, len(index_loc)):
            if i == 0:   	# Start of the genome
                # print "location i :", i, index_loc[i]
                diff = index_loc[i] - 1  # region before the first bait
                if diff > dist + bp_cov:
                    uncovered_target_reg.append((1, index_loc[i] - bp_neg - 1, (index_loc[i] - bp_neg - 1) - 1))
                diff = index_loc[i + 1] - index_loc[i]  # region between 1st n 2nd bait
                # print "location :", index_loc[i+1], index_loc[i]
            elif i == len(index_loc) - 1:  # End of the Genome
                diff = len_refseq - index_loc[i]
                # print index_loc[i], i
            else:
                diff = index_loc[i+1] - index_loc[i]  # In middle
                # print index_loc[i], i
            if diff > dist + bp_cov:
                if i == len(index_loc)-1:
                  uncovered_target_reg.append((index_loc[i] + bp_pos + 1, len_refseq, (len_refseq) - (index_loc[i] + bp_pos + 1)))
                else:
                  uncovered_target_reg.append((index_loc[i] + bp_pos + 1, index_loc[i + 1] - bp_neg - 1, (index_loc[i + 1] - bp_neg - 1) - (index_loc[i] + bp_pos + 1)))
    else:
        diff = index_loc[0] - 1
        if diff > dist + bp_cov:
            uncovered_target_reg.append((1, index_loc[0] - bp_neg-1, (index_loc[0] - bp_neg - 1) - 1))
        diff = len_refseq - index_loc[0]
        if diff > dist + bp_cov:
            uncovered_target_reg.append((index_loc[0] + bp_pos + 1, len_refseq, (len_refseq) - (index_loc[0] + bp_pos + 1)))

    print uncovered_target_reg
    return uncovered_target_reg


# ----------- Design coverage baits -----------------------------
def CoverageBaits(vir_nam, uncovered_target_reg, bait_dir, Ref_Seq, GC, blast_infile, bp_len, bp_cov, bp_neg, disRC):
    import os
    from Bio.Seq import Seq
    global vir_dir
    # global bp_pos, bp_neg  bp_len,
    global TM_Min, TM_Max
    print "\n9. Designing Coverage baits for untargeted region"
    Cov_B_file = vir_dir + "/" + vir_nam + "--bait.fa"

    [GC_Min, GC_Max] = GC.split('-')
    no_bait = []
    cbaits = []
    refSeq_cov = Seq(Ref_Seq[0])
    for reg in uncovered_target_reg:
        print "reg :", reg
        df = reg[2]/float(bp_cov)
        if df <= 1.1:
            start_kmer = reg[0] + int(reg[2]/2 - 60)  	# select a kmer in the middle of the range
            print "start-Loc: ", start_kmer
            kmer_s = str((refSeq_cov[start_kmer - 1: start_kmer - 1 + bp_len].reverse_complement()).upper())   	# 120 cDNA of target in upper case
            k = GC_TM_calc(kmer_s, start_kmer, GC, disRC)
            if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
                if disRC == 'False':    
                    cbaits.append((start_kmer, "Cov-" + GC, kmer_s, k[0], k[1]))
                    with open(Cov_B_file, 'a') as f:
                        f.write('>  Cov-bait |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + kmer_s + "\r\n")
                    with open(blast_infile, 'a') as f:
                        f.write('>  Cov-targ |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + str(Seq(kmer_s).reverse_complement()) + "\r\n")
                elif disRC == 'True':
                    cbaits.append((start_kmer, "Cov-" + GC, str(Seq(kmer_s).reverse_complement()), k[0], k[1]))
                    with open(Cov_B_file, 'a') as f:
                        f.write('>  Cov-bait |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + str(Seq(kmer_s).reverse_complement()) + "\r\n")
                    with open(blast_infile, 'a') as f:
                        f.write('>  Cov-targ |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + kmer_s + "\r\n")
            else:
                print "failed kmer:", k
                no_bait.append(reg)
        else:                               # multiple baits with a coverage of 500bp each
            reps = round(round((reg[2]/float(bp_cov)), 2))
            print 'Reps', reps, round((reg[2]/float(bp_cov)), 2)
            for i in range(0, int(reps)):
                start_kmer = bp_neg + reg[0] + bp_cov * i + i
                kmer_s = str((refSeq_cov[start_kmer - 1: start_kmer - 1 + bp_len].reverse_complement()).upper())
                k = GC_TM_calc(kmer_s, start_kmer, GC, disRC)
                if (k[0] <= int(GC_Max) and k[0] >= int(GC_Min)) and (k[1] <= TM_Max and k[1] >= TM_Min):
                    if disRC == 'False':    
                        cbaits.append((start_kmer, "Cov-" + GC, kmer_s, k[0], k[1]))
                        with open(Cov_B_file, 'a') as f:
                            f.write('>  Cov-bait |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + kmer_s + "\r\n")
                        with open(blast_infile, 'a') as f:
                            f.write('>  Cov-targ |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + str(Seq(kmer_s).reverse_complement()) + "\r\n")
                    elif disRC == 'True':
                        cbaits.append((start_kmer, "Cov-" + GC, str(Seq(kmer_s).reverse_complement()), k[0], k[1]))
                        with open(Cov_B_file, 'a') as f:
                            f.write('>  Cov-bait |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + str(Seq(kmer_s).reverse_complement()) + "\r\n")
                        with open(blast_infile, 'a') as f:
                            f.write('>  Cov-targ |' + str(GC) + ' | ' + str(start_kmer) + "\r\n" + kmer_s + "\r\n")

                else:
                    print "failed kmer:", k
                    no_bait.append(reg)
    print "The No bait regions :", no_bait
    return(cbaits, no_bait)
#############################
