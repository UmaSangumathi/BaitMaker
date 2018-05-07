#!/usr/bin/python
##############################################################
#  Script     : EnrichmentPipeline.py
#  Author     : uma sangumathi
#  Date       : 23/06/2015
#  Last Edited: 23/06/2015, uma 
#  Description: Wrapper script to run Enrichment downstream analysis for Segmented or not segmented virus
#     		Optinal to run in parallel [screen mode]
#		Please edit the path of the softwares 
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
import csv
from optparse import OptionParser


def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "usage: %prog [options]\nExample: ./EnrichmentPipeline.py -r `pwd`/Reference -q `pwd`/Fastq -e umasang7@gmail.com -d /share/ncbi/nt -o `pwd`/Test -i `pwd`/run_summary-tmp.csv -p F -s F -b `pwd` -t stringent"
    parser = OptionParser(usage=usage)
    parser.add_option("-r", "--refdir",
                      dest="ref_dir",
                      help="Directory of Reference fasta file")
    parser.add_option("-q", "--fastq_dir",
                      dest="fq_dir",
                      help="Directory of Fastq files")
    parser.add_option("-e", "--email_address",
                      dest="em",
                      help="Email address to access online NCBI")
    parser.add_option("-d", "--database",
                      dest="db",
                      help="NCBI nt database")
    parser.add_option("-o", "--outdir",
                      dest="out_dir",
                      help="Output directory")
    parser.add_option("-i", "--info",
                      dest="csv_file",
                      help="csv file containing the sample and reference information")
    parser.add_option("-p", "--parallel",
                      dest="par",
                      help="Run each fastq file in the info file in Parallel or not. Takes T or F")
    parser.add_option("-s", "--segment_type",
                      dest="seg",
                      help="Virus in the info file in Segmented or not. Takes T or F")
    parser.add_option("-b", "--code_dir",
                      dest="code_dir",
                      help="The directory were the scripts are present")
    parser.add_option("-t", "--stringency",
                      dest="param",
                      help="Stringency for the consensus genome generation: Takes value 'stringent' or 'lenient'. stringent = MIN_COV 10 and Mismatches bwa defualt; lenient = minimum coverage 1 and 20bp mismatch ")

    return parser

def main():
  """The main function
  """
  parser = cmdline_parser()
  try: 
    (opts, args) = parser.parse_args()
    print opts,args
    if len(args):
      parser.error("Unrecognized arguments found: %s." %(
                ' '.join(args)))
      sys.exit(1)
    else:
      threads = str(3)
      print opts.csv_file, args
    # Call the enrichment script according to the parameters
      os.chdir(opts.fq_dir)
    # Read the summary file:
      for xi in csv.reader(open(opts.csv_file, 'r')):
        x = xi[0]
        print "Processing : ", x
        os.environ["PICARDDIR"] = opts.code_dir + "/src/picard-tools-1.105"
        sys.path.append(opts.code_dir + "/bin")
        sys.path.append(opts.code_dir + "/src")
        print "PATH: ....",  sys.path
        print os.system("which EnrichmentAnalysis_Segmented.py")
        if opts.param == "stringent": cons_code = "bam2cons_iter-stringent.sh"
        elif opts.param == "lenient": cons_code = "bam2cons_iter-lenient.sh"
        else: 
	  print "Error : 'stringent' or 'lenient' "
          parser.print_help()
	  sys.exit(0)
        add_args = " " + opts.em + " " + opts.fq_dir + " " + opts.ref_dir + " " + opts.db + " " + opts.csv_file + " " + opts.out_dir + " " + cons_code + " " 
        print  add_args
        if opts.seg == "T": script = "EnrichmentAnalysis_Segmented.py"
        else: script = "EnrichmentAnalysis_NotSegmented.py"
        if opts.par == "T": cmd = "screen -d -m -S " +  x  + script + " " +  x  + add_args
        else: cmd = script + " " +  x  + add_args
        print "Running...."
        print cmd
        os.system(cmd)
  except:
    parser.print_help()
    sys.exit(0)

if __name__ == "__main__":
    main()

