#!/usr/bin/python

##############################################################
#  Script     : GetConsKmers.py
#  Author     : uma sangumathi
#  Date       : 28/06/2013
#  Last Edited: 
#  Description: Wrapper script to find conserved viral regions and design 120mer bait 
##############################################################
# Purpose: 
# Select diverse set of representative sequences from the database and make a fasta file
# Reverse complement of the fasta sequences
# Find conserved kmers using Primux software for different mismatches
#  
# Requirements:
#  1. Primux   
#  2. Vmatch
#  3. NCBI blast commandline tool and database
#  4. Sqlite 

#############################################################


# ------------- Set Default values and directory ------------
def set_dir(wrk_dir, proj_dir, mismatch, vir_nam):
    import os
    global exec_dir, primux_dir
    primux_dir = wrk_dir + '/src/primux'
    #primux_di = "/home/rongli/src/primux_uma"
    exec_dir = proj_dir + '/' + vir_nam
    if not os.path.isdir(exec_dir):
        os.mkdir(exec_dir)
        print "Created directory: ", exec_dir


# ---------- Create fasta files (1) Target sequence fasta (2) Complementary to Target sequence i.e. for bait sequence (primux processing) --------------
def write_to_fasta(u_id, numkmer, primux_flag, mismatch, vir_nam, GC):
    import sqlite3
    import os
    import os.path
    import shutil
    from Bio.Seq import Seq
    global exec_dir, primux_dir
    global cur1, con1
    if os.path.exists(exec_dir + '/' + vir_nam + "_cDNA.txt"):
        if primux_flag == "Yes":
            file_name = exec_dir + '/' + vir_nam + "_cDNA.txt"
            run_Primux(file_name, vir_nam, numkmer, mismatch, GC)
    else:
        for id in u_id:
            for row in cur1.execute("SELECT * from NCBI_DNA_records WHERE Unique_id=?", (id,)):
                Seqq = [str(a) for a in row]
                vir_nam = str(Seqq[0])
                vir_nam = vir_nam.replace(" ", "_")
                vir_nam = vir_nam.replace("/", "_")
                vir_nam = vir_nam.replace("'", "_")
                vir_nam = vir_nam.replace("(", "_")
                vir_nam = vir_nam.replace(")", "_")
                # print "vir_nam =" ,vir_nam
                Seq_fasta = ">" + str(Seqq[0]) + "| " + str(Seqq[1]) + '| ' + str(Seqq[3]) + "| " + str(Seqq[9]) + '| ' + str(Seqq[10]) + "| " + str(Seqq[11]) + '| ' + str(Seqq[12]) + "| " + str(Seqq[14]) + "\r\n" + str(Seqq[17]) + "\r\n"
                fast = Seq(str(Seqq[17]))
                compl_target = fast.reverse_complement()
                Seq_compl_fasta = ">" + str(Seqq[0]) + "| " + str(Seqq[1]) + '| ' + str(Seqq[3]) + "| " + str(Seqq[9]) + '| ' + str(Seqq[10]) + "| " + str(Seqq[11]) + '| ' + str(Seqq[12]) + " | " + str(Seqq[14]) + "\r\n" + str(compl_target) + "\r\n"		# write the sequence into fasta file with the header:: Species, unique_id, Seq_length, host, strain, collection date, tax_id, country, sequence
                fi = open(exec_dir + '/' + vir_nam + '_Target.txt', 'a')  # Store the original target squence in Fasta folder
                fi.write(Seq_fasta)
                fi.close()
                f = open(exec_dir + '/' + vir_nam + "_cDNA.txt", 'a')  # store complementary of the target sequence in fasta file for Primux Conserved kmers
                f.write(Seq_compl_fasta)
                f.close()
                # print ">> fasta file: update.. "

        if primux_flag == "Yes":
            file_name = exec_dir + '/' + vir_nam + "_cDNA.txt"
            run_Primux(file_name, vir_nam, numkmer, mismatch, GC)
        else:
            print "Hits less than 4: No PRIMUX ", vir_nam


# -------------- Run Primux for conserved kmers within the species on the complementary Target sequence ----------------
# PRIMUX:: Process fasta files, compute suffix array, filter kmers by UNAfold, get the final set of conserved and unique probes with target coverage
def run_Primux(file_name, vir_nam, numkmer, mismatch, GC):
    import shutil
    import os
    import subprocess
    global exec_dir, primux_dir
    wrk_dir = exec_dir + '/GC_' + GC
    if not os.path.isdir(wrk_dir):
        os.mkdir(wrk_dir)
    os.chdir(wrk_dir)
    [GC_min, GC_max] = GC.split('-')
    shutil.copy2(primux_dir + '/options', wrk_dir + '/options_sp_C')
    shutil.copy2(primux_dir + '/run_Primux_tmp.sh', wrk_dir + '/run_Primux_sub.sh')
    lines = open(wrk_dir + '/options_sp_C', 'r').readlines()
    lines[0] = lines[0].replace('-dir=results', '-dir=' + wrk_dir + '/conserved_res-' + mismatch)
    lines[1] = lines[1].replace('-input_sequences=tmp.fa', '-input_sequences=' + wrk_dir + '/indx_nam')
    lines[2] = lines[2].replace('-topN_min=0', '-topN_min=' + numkmer)
    lines[3] = lines[3].replace('-topN=1', '-topN=1')
    lines[4] = lines[4].replace('-bottomN=0', '-bottomN=0')
    lines[5] = lines[5].replace('-topN_pp=0', '-topN_pp=' + numkmer)
    lines[6] = lines[6].replace('-topN_ss=0', '-topN_ss=' + numkmer)
    lines[7] = lines[7].replace('-max_mm=0', '-max_mm=' + mismatch)
    lines[8] = lines[8].replace('-min_percent_gc=0', '-min_percent_gc=' + GC_min)
    lines[9] = lines[9].replace('-max_percent_gc=0', '-max_percent_gc=' + GC_max)
    open('options_sp_C', 'w').write(''.join(lines))
    print "Runing PrimuX Software:"
    subprocess.call(['./run_Primux_sub.sh', file_name])


# -------------------- Similarity Score -------------------
def find_Similarity(sp, hit_details, SIM_FILTER, MAX_ncbihits):
    global cur1, con1
    HIT_FILTER = 4
    X = []
    comSet = []
    len_hits = len(hit_details)
    for i in range(len_hits - 1):   # Pairwise comparison between list of tuples to dicard similar sequences increase diversity
        for j in range((i + 1), (len_hits)):
            if j < (len_hits):
                comSet.append(str(hit_details[i][0]) + " : " + str(hit_details[j][0]))             # Store the Unique id's of of pairwise comparison
                intsec = len(set(hit_details[i]).intersection(set(hit_details[j])))         # Indentity or Similarity Score based on the common elements in tuple
                X.append(intsec)
                if intsec >= SIM_FILTER:            # discard a hit if indentity is 6 or 7 in the list of tuple
                    check_ref = cur1.execute("SELECT Other_Seq_id from NCBI_DNA_records WHERE Unique_ID=?", (hit_details[j][0],))
                    Ref = [str(a) for a in check_ref]
                 #   print "Ref::: ", Ref
                    Ref_flag = [e for e in Ref if e.startswith('ref|N') or e.startswith('ref|A')]           # if the hit is a reference genome then discard the other pair being compared
                 #   print Ref_flag
                    if len(Ref_flag) == 0:
                 #       print "Deleting", hit_details[j][0]
                        hit_details.remove(hit_details[j])
                        len_hits = len_hits - 1
                    else:
                 #       print "Ref_seq found... Deleting", (hit_details[i])
                        hit_details.remove(hit_details[i])
                        len_hits = len_hits - 1

    print 'Re hit_details: simfilter', len(hit_details), SIM_FILTER
    return [hit_details, SIM_FILTER, X, comSet]


# ---------------------------------------------------
def filter_NCBIhits(sp, sp_nam, SIM_FILTER, mismatch, Numkmer, hm_dir, GC):
    import sqlite3
    import pickle
    global exec_dir
    global Genome_length, con1, cur1
    HIT_FILTER = 4          # Numbers of hits per species required to run Primux
    MAX_ncbihits = 125
    hit_details_org = []

    print 'NCBI filter on ...', sp
    max_len_hit = cur1.execute("SELECT * from NCBI_DNA_records  WHERE Species=? ORDER BY Sequence_length DESC limit 1", (sp,))
    for dd in max_len_hit:
        max_F = [str(a) for a in dd]
        Genome_length = int(max_F[3])
        min_len = (Genome_length * 0.65) 			# min_len of 3/4  maximum length of the genome; hits with seq length below this is discarded
        print 'The range of squence length selected is from :', min_len, '--', max_F[3]
    print cur1.execute("SELECT * from NCBI_DNA_records  WHERE Species=? ORDER BY Sequence_length DESC limit 1", (sp,))
    sub_hit = cur1.execute("SELECT * from NCBI_DNA_records  WHERE Species=? AND Sequence_length>=? ORDER BY Sequence_length DESC",  (sp, min_len))
    ai = 0
    for rec in sub_hit:
        itm = [str(a) for a in rec]
        ai = ai + 1
        set_tup = (itm[1], itm[3], itm[9], itm[10], itm[11], itm[13], itm[13], itm[14], itm[12])
        # Store the hit details tot:8 values (unique_id, Sequence length, host, strain, collection date, Isolate, country, dbRef(taxon:id)) in tuple for comparison
        hit_details_org.append(set_tup)

    if len(hit_details_org) > HIT_FILTER:
        num_hits = len(hit_details_org)
        [hit_details, Sim_filter, X, comSet] = find_Similarity(sp, hit_details_org, SIM_FILTER, MAX_ncbihits)

        if len(hit_details) > MAX_ncbihits and Sim_filter > 3:
            num_hits = len(hit_details)
            print 'hit_details: simfilter', len(hit_details), Sim_filter
            status = "repeat"
            while status == "repeat":
                [hit_details, Sim_filter, X, comSet] = find_Similarity(sp, hit_details, Sim_filter - 1, MAX_ncbihits)
                if len(hit_details) <= MAX_ncbihits or Sim_filter == 3:
                    status = "stop"
                else:
                    status = "repeat"

        if len(hit_details) <= MAX_ncbihits or Sim_filter == 3:
            if len(hit_details) < HIT_FILTER:
                primux_flag = "No"
            else:
                primux_flag = "Yes"
            num_hits = len(hit_details)
            ID_list = [x[0] for x in hit_details]
            write_to_fasta(ID_list, Numkmer, primux_flag, mismatch, sp_nam, GC)			# Calls function to write the selected sequence to file and proceed to primux
            print "Storing the selected hit ids and similarity score...", Sim_filter
            list_ins = (sp_nam, ai, str(comSet), str(X), str(ID_list), str(len(ID_list)))
            with open(exec_dir + '/' + sp_nam + '-hit_dict.txt', 'wb') as handle:
                pickle.dump(list_ins, handle)
        else:
            print "U r not possible.... LOGIC ERROR"
            exit
    else:
        primux_flag = "No"
        num_hits = len(hit_details_org)
        ID_list = [x[0] for x in hit_details_org]
        write_to_fasta(ID_list, Numkmer, primux_flag, mismatch, sp_nam, GC)
    return (primux_flag, Genome_length, num_hits)


# ---------------------------------------------------
# Main code::
def getConsKmers(sp, sp_nam, hm_dir, proj_dir, mismatch, topNumkmer, GC, con, cur):
    global con1, cur1
    con1 = con
    cur1 = cur

    set_dir(hm_dir, proj_dir, mismatch, sp_nam)
    print "Calling NCBI filter........."
    SIM_FILTERG = 6          # Selecting the hits with diverse set of features (6 out of 8 total features indicate 6 similar features **unique_id included)
    [primux_flag, Genome_length, num_hits] = filter_NCBIhits(sp, sp_nam, SIM_FILTERG, mismatch, topNumkmer, hm_dir, GC)
    return (primux_flag, Genome_length, num_hits)
