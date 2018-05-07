#!/bin/bash

# export the required paths to run PRIMUX : vmatch, mkvtree, unafold
export PRIMUX_DIR=/home/rongli/src/primux_uma/bin
export LD_LIBRARY_PATH=/home/rongli/src/primux_uma/bin
export UNAFOLDDAT=/home/rongli/src/primux_uma/src/unafold_data
export PATH=$PATH:/home/rongli/src/vmatch-2.2.4-Linux_x86_64-64bit

primux_dir=/home/rongli/src/primux_uma

if [ $# -ne 1 ]
then
        echo Usage: $0 file_process
        exit 1
else
        fasta_file=$1                    #file name
fi


$primux_dir/bin/prep_fasta_file.pl $fasta_file indx_nam 
$primux_dir/bin/conseredKmers_us.py  options_sp_C
