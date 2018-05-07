#!/usr/bin/python

import sys
import os
import time
import datetime

# import custom modules
from GenomeBaitProcess_NS import *  # Functions related to Baits designing from the kmers
from GetConsKmers_NS import *  # Functions related to find the conserved kmers
# from GetConsKmers_NS_top150 import *
from BlastParser_NS import *


# ------------- Set Directories -----------------
def setDirectory(virus, virus_nam, mismatch, GC, main_dir, db_nam):
    import sqlite3
    import os
    print "Setting Directory and database connection"
    global wrk_dir, blast_dir, db_dir, virus_dir, query_file, bait_dir
    global num_seq, virus_index, virus_seq_count, blast_infile
    global con1, cur1

    print "VIRUS :  ", virus
    hm_dir = main_dir + "/Data"
    wrk_dir = hm_dir + "/Process"
    if not os.path.exists(wrk_dir):
        os.makedirs(wrk_dir)
    db_dir = main_dir + "/DB"
    if not os.path.isdir(db_dir):
        os.mkdir(db_dir)
    virus_dir = wrk_dir + '/' + virus_nam
    if not os.path.isdir(virus_dir):
        os.mkdir(virus_dir)
    query_file = virus_dir + '/GC_' + GC + "/conserved_res-" + str(mismatch) + "/queries_vm_supersets"
    virus_index = virus_dir + '/GC_' + GC + "/indx_nam.table"
    bait_dir = hm_dir + "/BaitSummary"
    if not os.path.isdir(bait_dir):
        os.mkdir(bait_dir)
    blast_dir = hm_dir + "/BlastRes"
    if not os.path.isdir(blast_dir):
        os.mkdir(blast_dir)
    blast_infile = virus_dir + "/" + virus_nam + "_blastin.fa"
    graph_dir = hm_dir + "/Graphs"
    if not os.path.isdir(graph_dir):
        os.mkdir(graph_dir)

    # ------------- Connect to databases --------------
    con1 = sqlite3.connect(db_dir + '/' + db_nam)  # VirusCaptureBaits.db => includs vector ; VirusCapture_run.db -> no clon r vector
    cur1 = con1.cursor()


# ------------------ Main Function -------------
def MainFn(virus, vir_nam, main_dir, db_nam, bp_len, bp_cov, bp_pos, bp_neg, disRC):
    from operator import itemgetter
    import csv
    global cur1, con1
    global query_file, blast_infile, virus_index

    pri_cons_baits = []
    num_baits = 0
    expt_Spcons = None
    cont = 'T'
    virus_id = None
    dist = 250  # distance between two coverage region tat is allowed
    GC_sets = ["35-65", "20-80"]

    for gc in GC_sets:  # GC range
        for m_topn in ["6-3500", "12-2500", "24-1000", "30-500"]:  # Indentity : 95%, 90%, 80%, 75% ["6-5000", "12-4000", "24-3000", "30-1000"]
            [m, topN] = m_topn.split('-')
            if cont != 'F':  # continue till the num of Sp conserved baits is less than 60% of the genome length
                setDirectory(virus, vir_nam, m, gc, main_dir, db_nam)

                cmd = getConsKmers(virus, vir_nam, main_dir, wrk_dir, str(m), str(topN), gc, con1, cur1)
                print cmd

                [flag, Genomelength, num_hits] = getConsKmers(virus, vir_nam, main_dir, wrk_dir, str(m), str(topN), gc, con1, cur1)
                # Genomelength = 18794
                # flag = "Yes"
                if flag == "Yes":
                    expt_Spcons = Genomelength * 0.60/bp_cov  # bp_cov default 500
                    expt_Cov = Genomelength * 0.40/bp_cov  # bp_cov default 500
                    print "Expected Conserved baits: ", expt_Spcons
                    print "Expected Coverage baits: ", expt_Cov
                    kmer_dict = FileToDict(query_file, m)
                    print "No. of distinct kmers in the query file - kmer_dict: ", len(kmer_dict)
                    tot_hits = MergeInfo(kmer_dict)
                    print "No. of kmer after merging kmer information - tot_hits: ", len(tot_hits)

                    for i in [90.0, 80.0, 70.0, 60.0]:
                        key = str(m) + '-' + str(i) + '-' + gc
                        print "\n\n######## ", key, " ###########"
                        kmer_cons_filt = HitsCoverage(tot_hits, i, virus_index)
                        print "No. of kmers that cleared coverage filer on target - kmer_cons: ", len(kmer_cons_filt)
                        if len(kmer_cons_filt) > 0 and num_baits < expt_Spcons:
                            if virus_id is None:
                                [Ref_seqloc, RefUniqueId, ref_u_id, virus_id, Ref_Seq] = RefMapping(kmer_cons_filt, cur1)
                            # print "Reference Genome Found : ", RefUniqueId
                            else:
                                Ref_seqloc = Map_knownRef(kmer_cons_filt, virus_id)
                            sort_u_loc = LeastDegenerate(Ref_seqloc, m)
                            comb_seq_loc = seqToloc(sort_u_loc)
                            passed_kmer = testGC_Tm(comb_seq_loc, bait_dir, vir_nam, key, gc, disRC)

                            if len(pri_cons_baits) == 0:
                                pri_index_loc = []
                                pri_untar_reg = []
                                [index_loc, baits] = RemoveOverlapBaits(passed_kmer, key, virus_dir, vir_nam, str(RefUniqueId), blast_infile, bp_cov)
                                print index_loc
                                if len(index_loc) != 0:
                                    uncovered_target_reg = UntargetRegion(index_loc, dist, bp_len, bp_cov, bp_pos, bp_neg)
                                    pri_cons_baits = baits
                                    pri_index_loc = index_loc
                                    pri_untar_reg = uncovered_target_reg
                                    print "The pri untarget  region ....", pri_untar_reg
                                    num_baits = len(pri_cons_baits)
                            else:
                                [n_index_loc, baits] = RmvOverUntar(passed_kmer, pri_untar_reg, key, str(RefUniqueId), dist, blast_infile, bp_len, bp_cov, disRC)
                                print "n_index_loc:", n_index_loc
                                print "pri_index_loc:", pri_index_loc
                                pri_index_loc.extend(n_index_loc)
                                print pri_index_loc.sort()
                                print "pri_index_loc.extend ", pri_index_loc
                                print "New Baits discovered", baits
                                pri_cons_baits.extend(baits)
                                num_baits = len(pri_cons_baits)
                                pri_untar_reg = UntargetRegion(pri_index_loc, dist, bp_len, bp_cov, bp_pos, bp_neg)
                                print "Changed pri untarget....", pri_untar_reg
                                print "\n CONS BAITS.......", pri_cons_baits

                            print "\n\nExpected Conserved baits: ", expt_Spcons
                            print "Expected Coverage baits: ", expt_Cov
                            print "No. of Species conserved baits designed: ", len(pri_cons_baits)
                            if num_baits >= expt_Spcons:
                                print "60% of genome covered by conserved probes.... "
                                cont = 'F'
                else:
                    print "NO Baits designed........."
                    exit

    if flag == "Yes" and len(pri_cons_baits) != 0:
        print "Conserved Baits: \n", pri_cons_baits
        pri_cov_baits = []
        Final_baits = []
        for GC in GC_sets:
            key = GC
            [Cov_baits, no_bait] = CoverageBaits(vir_nam, pri_untar_reg, bait_dir, Ref_Seq, GC, blast_infile, bp_len, bp_cov, bp_neg, disRC)
            All_bait = pri_cons_baits + pri_cov_baits + Cov_baits
            index_loc_tm = [x[0] for x in All_bait]
            index_loc_tm.sort()
            pri_untar_reg = UntargetRegion(index_loc_tm, dist, bp_cov, bp_pos, bp_neg)
            pri_cov_baits = Cov_baits
        print "Coverage baits: ", Cov_baits
        refseq_fil = blast_dir + '/' + vir_nam + '-ref.txt'
        with open(refseq_fil, 'a') as f:
            f.write(str(RefUniqueId))
        Final_baits = baitBlast(vir_nam, blast_infile, All_bait, blast_dir, refseq_fil, str(RefUniqueId))
        ref_length = cur1.execute("SELECT Sequence_length from NCBI_DNA_records WHERE Unique_id=?", (RefUniqueId, ))
        ref_len = [i[0] for i in ref_length][0]
        Final_baits_sort = sorted(Final_baits, key=itemgetter(0))
        Final_baits_sort.insert(0, ('Start_location', 'Bait_type', 'Bait_Sequence', 'GC_content', 'Melting_temperature', 'Blast_Ref', 'Hits_Ref', 'Blast_human', 'Hits_human', 'Blast_bacteria', 'Hits_bacteria', 'Blast_rat', 'Hits_rat', 'Blast_bats', 'Hits_bats', 'Blast_chicken', 'Hits_chicken'))
        resultFile = bait_dir + '/' + vir_nam + '--' + str(RefUniqueId) + '--' + str(ref_len) + "--" + str(num_hits) + "--RES.csv"
        write_csv(Final_baits_sort, resultFile)
        # generateGraph(RefUniqueId, resultFile, virus_nam)

    else:
        print "Coverage Baits for whole genome: \n"
        pri_cov_baits = []
        Final_baits = []
        for GC in GC_sets:
            key = GC
            [RefUniqueId, RefSeq] = GetReference(virus, cur1, virus_dir)
            if len(pri_cov_baits) == 0:
                pri_untar_reg = [(1, len(RefSeq[0]), len(RefSeq[0]) - 1)]
            [Cov_baits, no_bait] = CoverageBaits(vir_nam, pri_untar_reg, bait_dir, RefSeq, GC, blast_infile, bp_len, bp_cov, bp_neg, disRC)
            All_bait = pri_cov_baits + Cov_baits
            index_loc = [x[0] for x in All_bait]
            index_loc.sort()
            print "All baits", All_bait
            if len(index_loc) != 0:
                pri_untar_reg = UntargetRegion(index_loc, dist, bp_len, bp_cov, bp_pos, bp_neg)
                pri_cov_baits = Cov_baits
        print "Coverage baits: ", Cov_baits
        ref_length = cur1.execute("SELECT Sequence_length from NCBI_DNA_records WHERE Unique_id=?", (RefUniqueId, ))
        ref_len = [i[0] for i in ref_length][0]
        refseq_fil = blast_dir + '/' + vir_nam + '-ref.txt'
        with open(refseq_fil, 'a') as f:
            f.write(str(RefUniqueId))
        # Final_baits = baitBlast(vir_nam, blast_infile, All_bait, blast_dir, refseq_fil, str(RefUniqueId))
        Final_baits = All_bait  # line is added to bypass blast filter
        Final_baits_sort = sorted(Final_baits, key=itemgetter(0))
        Final_baits_sort.insert(0, ('Start_location', 'Bait_type', 'Bait_Sequence', 'GC_content', 'Melting_temperature', 'Blast_Ref', 'Hits_Ref', 'Blast_human', 'Hits_human', 'Blast_bacteria', 'Hits_bacteria', 'Blast_rat', 'Hits_rat', 'Blast_bats', 'Hits_bats', 'Blast_chicken', 'Hits_chicken'))
        resultFile = bait_dir + '/' + vir_nam + '--' + str(RefUniqueId) + '--' + str(ref_len) + "--" + str(num_hits) + "--RESc.csv"
        write_csv(Final_baits_sort, resultFile)

        # generateGraph(RefUniqueId, resultFile, virus_nam)


def generateGraph(RefUniqueId, resultFile, virus_nam, genome_len):
    arg_list = "--args " + graph_dir + " " + RefUniqueId + " " + resultFile + " " + virus_nam + " " + genome_len
    cmd = "Rscript script.R " + arg_list
    subprocess.call(cmd, shell=True)


# ************************   main script   **************************

if len(sys.argv) == 7:
    print str(sys.argv)
    sp = str(sys.argv[1])
    print sp
    main_dir = str(sys.argv[2])
    sp_nam = sp.replace(' ', '_')
    sp_nam = sp_nam.replace('/', '_')
    sp_nam = sp_nam.replace("'", "_")
    sp_nam = sp_nam.replace("(", "_")
    sp_nam = sp_nam.replace(")", "_")

    db_nam = str(sys.argv[3])
    bp_len = int(sys.argv[4])
    bp_cov = int(sys.argv[5])
    bp_pos = (bp_cov / 2) + ((bp_len / 2) - 1)  # default 309, (500/2)+59
    bp_neg = (bp_cov / 2) - (bp_len / 2)  # default 190, (500/2)-60
    disRC = str(sys.argv[6])

    # ------------- Open Logfile -------------
    op_dir = main_dir + "/SummaryLOGS/OP_NS"
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    log_file = open(op_dir + "/OP-" + sp_nam + ".log", "w")
    old_stdout = sys.stdout
    sys.stdout = log_file
    print "\n\n------------------- Summary Log -----------------\n\n"
    print "******************************************* ", sp, " *******************************************\n"

    start_time = time.time()
    MainFn(sp, sp_nam, main_dir, db_nam, bp_len, bp_cov, bp_pos, bp_neg, disRC)
    time_sec = time.time() - start_time
    time_hr = str(datetime.timedelta(seconds=time_sec))
    print "The time taken for Prog execution: ", time_hr
    # close logfile....
    sys.stdout = old_stdout
    log_file.close()
else:
    print "Error: Please enter the correct arguments."
    print str(sys.argv)
