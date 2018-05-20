# PathNrich
updated : Dec 2015
 
## Design baits targetting the genome for enrichment (120bp)  
### A) Conserved baits: 

runVirusDNA_NS.py = Main wrapper script that calls different functions to design conserved baits
 
### B) Exhaustive baits: 
 
Scripts:   
RunGeneratebaits.py  = Main wrapper script that calls different functions to design optimal number of baits for the given set of sequences.  
DesignMinimalBaits.py = Contains the functions to design baits   
 
To Run:   
RunGeneratebaits.py  'Directory'  'Species Conserved baits fasta file'  'Query reference genomes in fasta'   
 
Notes:   
'Directory': the initial files required for the designing and the output files are written  
'Species Conserved baits fasta file': fasta file with designed conserved baits designed by option A. This is to make sure that the Species   consevered baits already designed are included.  
'Query reference genomes in fasta': This file is a fasta file with all the genome sequences for which the baits are designed (.fa)   
Each step generated a number of intermediate files and these are deleted in the final step   

DesignMinimalBaits.py. Hence if you have any other file with the following terms they will be deleted.  json , finalbaits, Optimized, Blast, blast, Redesign, txt , Untar   
Output files:  
*-All_BAITS.fa: Fasta file with all the baits sequences   
*-All_BAITS-mapped.csv: File contains the information on where the baits map to each of the genome sequence used.   
*-All_BAITS-untargeted.csv: Untargeted region of the genome where the adjacent baits are in a distance of more than 500bp apart from each other.  
 

-------------------------------------------------------------------------------------------------------- 
###  C) Analysis Pipeline    
 
Procedure to run Enrichment Analysis Pipeline:   
 
Requirements:    
bin folder = Contains all the scripts   
src folder = Contains all the softwares   
EnrichmentAnalysisPipeline.py  
run_QC.sh  
 
Before you run the scripts export the folders    
export PATH=pathtofolder/bin:$PATH   
export PATH=pathtofolder/src:$PATH   
 
 
I) Run QC   
Have all the raw fastq files in a single folder and run the following script. The script uses fastqc and trim_galore software to do the quality filter of 20 and minimum bp as 35. Make sure the software is imported before you run the script.   
 
./run_QC.sh  'Raw fastq files directory'   
 
Find Unknown species:   
./FindSpeciesInSample.py  'Fastq directory which are quality checked'  'NCBI nt database location'    
Export the XML file in Megan software   
 
II) Find consensus genome, mapping and snps   
 
The "EnrichmentAnalysisPipeline.py" is a wrapper script which depending upon the parameters calls either   
1) EnrichmentAnalysis_NotSegmented.py   
2) EnrichmentAnalysis_Segmented.py   
 
The sample information and the fasta file to be used as a reference is given in a csv file:   
CD-Lib-0004_S4_L001_R1_001_val_1.fq,Rotv-CU938.fa   
CD-Lib-0005_S5_L001_R1_001_val_1.fq,Rotv-YR062.fa   
 
where, CD-Lib-0004_S4_L001_R1_001_val_1.fq is the forward read fastq file present in Fastq directory and Rotv-CU938.fa is the reference file with 11 segments in Reference directory. Please note that only the forward read fastq file should be given in the file, the code will find its pair by searching "R2" instead of R1.   
 
For parallel run, for every fastq file in the csv file a ‘screen’ is created and the pipeline is run. You can check the jobs running by screen -ls command.   
 
Usage: EnrichmentPipeline.py [options]   
Example: ./EnrichmentPipeline.py -r `pwd`/Reference -q `pwd`/Fastq -e gmail.com -d /share/ncbi/nt -o `pwd`/Test -i `pwd`/run_summary-tmp.csv -p F -s F -b `pwd` -t stringent   
 
Options: 
  -h, --help            show this help message and exit 
  -r REF_DIR, --refdir=REF_DIR 
                        Directory of Reference fasta file 
  -q FQ_DIR, --fastq_dir=FQ_DIR 
                        Directory of Fastq files 
  -e EM, --email_address=EM 
                        Email address to access online NCBI 
  -d DB, --database=DB  NCBI nt database 
  -o OUT_DIR, --outdir=OUT_DIR 
                        Output directory 
  -i CSV_FILE, --info=CSV_FILE 
                        csv file containing the sample and reference information 
 -p PAR, --parallel=PAR 
                        Run each fastq file in the info file in Parallel or not. Takes T or F 
-s SEG, --segment_type=SEG 
                        Virus in the info file in Segmented or not. Takes T or F 
-b CODE_DIR, --code_dir=CODE_DIR 
                        The directory were the scripts are present 
 -t PARAM, --stringency=PARAM 
                        Stringency for the consensus genome generation: Takes 
                        value 'stringent' or 'lenient'. stringent = MIN_COV 10 
                        and Mismatches bwa defualt; lenient = minimum coverage 
                        1 and 20bp mismatch 
 
Output: 
The pipeline outputs, 
i) Consensus genome generated: *-hybrid.fa and *-hybrid.annot 
ii) SNPs found: *.vcf   
iii) Mapped bam file: *-sort.bam  
iv) Per position coverage: *.cov  
 
 
III) Generate graphs  
runGraphs.py is the wrapper script that calls generated_circosPlot.py (to generate the graph) and bait_effiency_modules.py (For plotting the baits position)  
 
Files required for Graph..   
Sample related:   
Create a new directory and copy all the *vcf , *-hybrid.fa , *.cov from the previous output.  
Genome related:    
Rotv-CU938.coords (Gives the genome structure [represent as single chromosome or segmented])  
Rotv-CU938.genes (Gives the names of the segments or if the genome is split into genes then gene names)  
Rotv-CU938.baits: If it’s a capture and baits are present, this file should contain the baits in fasta format. The baits are then mapped  to the *hybrid.fa and positions are obtained. If it’s the Library, please make an empty file with name is Rotv-CU938.baits.   
Software related:   
Copy all the *.conf files from the ./bin folder to this directory. The “hist.conf” is the configuration file for the layers in the plot (e.g. Coverage layer, baits layer, snps layer etc..), change the display parameters in this if needed. ideogram.conf gives the config for the outer genome display.  
Copy the log_axis.txt from the ./bin folder and edit it accordingly. By default the log axis is marked from 10^1 to 10^5. The changes for this is also in the hist.conf   
 
consensus genome fasta file headers should be used as the genome or segment name in the graph related files (e.g. Rotv-CU938.coords).   
 
* *These scripts are not packaged. Hence please edit the parameters and values within the script.  
