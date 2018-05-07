#!/bin/bash

# Andreas Wilm <wilma@gis.a-star.edu.sg>
# Copyright: 2012, 2013 Genome Institute of Singapore
# License: GPL2


# keep in sync with arg parsing below
keep_last_bam=0
illumina=0
threads=4
force=0
usage() {
cat <<EOF
$(basename $0): Derive a new consensus seq through iterative mapping

Simple bash script that calls bam2cons.py repeatedly until the maximum
number of iterations is reached or not more maps than before map.
Similar in design to iCORN (Otto et al., 2010) PMID 20562415, but faster.

  Options:
    -f | --fastq1   : Input fastq[.gz] file
    -g | --fastq2   : Fastq[.gz], second in pair (optional)
    -h | --help     : Display this help
       | --illumina : Phred qualities are ASCII64, i.e. Illumina 1.3-1.7.
                      Check with e.g. FastQC if unsure.
                      BWA will assume Sanger/Illumina 1.8+ otherwise,
                      which is default for Casava>=1.82
    -k | --keep     : Keep last BAM (mainly for debugging)
    -r | --ref      : Initial reference fasta file
    -t | --threads  : Number of threads to use (default=$threads)
    -o | --outref   : The output reference
         --force    : Force overwriting of files
EOF
}


while [ "$1" != "" ]; do
    case $1 in
        -f | --fastq1 )
            shift
            fastq1=$1
            ;;
        -g | --fastq2 )
            shift
            fastq2=$1
            ;;
            --illumina )
            illumina=1
            ;;
        -r | --ref )    
            shift
            init_ref=$1
            ;;
        -o | --outref )
            shift
            out_ref=$1
            ;;
        -k | --keep )
            keep_last_bam=1
            ;;
        -t | --threads )
            shift
            threads=$1
            ;;
        -h | --help )
            usage
            exit
            ;;
             --force )
            force=1
            ;;
        * ) 
            echo "FATAL: unknown argument \"$1\""
            usage
            exit 1
    esac
    shift
done


if [ -z $init_ref ]; then
    echo "FATAL: missing initial reference fasta argument" 1>&2   
    echo
    usage
    exit 1
fi
if [ -z $out_ref ]; then
    echo "FATAL: missing output file argument" 1>&2
    echo
    usage
    exit 1
fi
if [ -z $fastq1 ]; then
    echo "FATAL: missing fastq argument" 1>&2
    echo
    usage
    exit 1
fi

if [ -s $out_ref ] && [ $force -ne 1 ]; then
    echo "FATAL: refusing to overwrite $out_ref" 1>&2
    exit 1
fi
if [ ! -e $fastq1 ]; then
    echo "FATAL: non-existing file $fastq1" 1>&2
    exit 1
fi
if [ ! -z $fastq2 ]; then
    if [ ! -e $fastq2 ]; then
        echo "FATAL: non-existing file $fastq2" 1>&2
        exit 1
    fi
fi
if [ ! -e $init_ref ]; then
    echo "FATAL: non-existing file $init_ref" 1>&2
    exit 1
fi


# ----------

# check for existance of external tools we rely on.
# 
which bwa_unique-lenient.sh >/dev/null || exit 1
which samtools >/dev/null || exit 1
which bam2cons-lenient.py >/dev/null || exit 1
bam2cons-lenient.py -h >/dev/null
if [ $? -ne 0 ]; then
    echo "FATAL: dry run of bam2cons.py failed. Is pysam installed?" 1>&2
    exit 1
fi


# ----------

# prevent accidental overwriting of files
set -o noclobber

MAX_ITER=10
DEBUG=0

tmpdir=$(mktemp -d "$(basename $0).XXXXXX")

# ----------

function fa_to_md5 {  
    grep -v '^>' | tr -d '\n' | md5sum  | cut -f 1 -d ' '
}


if [ $(grep -c '^>' $init_ref) -ne 1 ]; then
    echo "FATAL: can only process fasta files with one sequence/chromosome" 1>&2
    exit 1
fi

echo "INFO: saving temporary results to $tmpdir"

# copy refseq and change id which is used for reconstruction
old_ref=$tmpdir/$(basename $init_ref)
newid=$(basename $out_ref .fa | sed -e 's,^\.*,,')
cat $init_ref | sed -e "/^>/s/>.*/>${newid}/" > $old_ref


old_num_mapped=0
iter=1
while [ $iter -le $MAX_ITER ]; do
    outprefix=$tmpdir/newref.iter-$iter
    bam=$outprefix.bam

    echo "INFO: iteration $iter: aligning to $old_ref"

    bwa_unique_args="-f $fastq1"
    if [ ! -z $fastq2 ]; then
        bwa_unique_args="$bwa_unique_args -g $fastq2"
    fi
    if [ $illumina -eq 1 ]; then
        bwa_unique_args="$bwa_unique_args --illumina"
    fi
    bwa_unique_args="$bwa_unique_args -r $old_ref"
    bwa_unique_args="$bwa_unique_args -t $threads"
    bwa_unique_args="$bwa_unique_args -o $bam"
    echo "INFO: calling bwa_unique-lenient.sh with the following args: $bwa_unique_args"
    #test $DEBUG -eq 1 && echo "DEBUG iter $iter: bwa_unique_args=$bwa_unique_args" 1>&2
    bwa_unique-lenient.sh $bwa_unique_args || exit 1
    samtools index $bam || exit 1
    # extract new ref
    #
    echo "INFO: iteration $iter: extracting new consensus"
    new_ref=${outprefix}.fa
    seq=$(samtools view -H $bam | grep '^@SQ' | cut -f 2 | cut -f 2 -d :)
    test $DEBUG -eq 1 && echo "DEBUG iter $iter: calling bam2cons-lenient.py -b $bam -r $seq -o $new_ref --no-cov-char N" 1>&2
    bam2cons-lenient.py -b $bam -r $seq -o $new_ref --no-cov-char "N" -v || exit 1

    # exit if new and old ref are the same
    #
    old_md5=$(cat $old_ref | fa_to_md5)
    new_md5=$(cat $new_ref | fa_to_md5)
    if [ $old_md5 == $new_md5 ]; then
        echo "INFO iteration $iter: no more changes in consensus sequence"
        break
    fi

    new_num_mapped=$(samtools idxstats $bam | awk '{print $3; exit}')
    if [ $(($new_num_mapped-$old_num_mapped)) -lt 1000 ]; then
	    echo "INFO: iteration $iter: Not enough improvement after re-aligning"
	    break
    elif [ -z $new_num_mapped ] || [ $new_num_mapped -eq 0 ]; then
        echo "FATAL: numbers of mapped reads = $new_num_mapped" 1>&2
        exit 1
    else
	    echo "INFO: iteration $iter: Got $new_num_mapped mapped reads (previously $old_num_mapped)"
    fi

    old_num_mapped=$new_num_mapped
    old_ref=$new_ref
    let iter=iter+1
done
if [ $iter == $MAX_ITER ]; then
    echo "INFO: maximum number of iterations ($MAX_ITER) reached."
fi


echo "INFO: cleaning up and copying final reference to $out_ref"

# id was already updated. now just copy
#cat $new_ref | sed -e "/^>/s/>.*/>${newid}/" > $out_ref || exit 1
cp $new_ref $out_ref || exit 1
if [ $keep_last_bam -eq 1 ]; then
    #echo "WARN: keeping old bam. remember that the used refname might be different" 1>&2
    cp $bam $(dirname $out_ref)/
fi
rm -rf $tmpdir
