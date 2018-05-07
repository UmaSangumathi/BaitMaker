#!/usr/bin/env python
"""Extract a consensus sequence from a read mapping (BAM) based on
simple majority rules. Makes no ploidy assumptions (initially designed
for viral genomes). Able to handle indels. Columns with no coverage
can be specially marked (see option --no-cov-char).
"""

#--- standard library imports
#
from __future__ import division
import os
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
from collections import Counter

#--- third-party imports
#
import pysam
# FIXME SeqIO only used for id parsing. unnecesssary!
from Bio import SeqIO

#--- project specific imports
#
# /



__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013 Genome Institute of Singapore"
__license__ = "GPL2"
__credits__ = [""]
__status__ = "eternal alpha"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')

# Don't make this too small! I've seen high coverage cases, where the
# first 1000 nucleotides were N's the other X thousands were A. Too
# high numbers will make this veeeery slow in ultra-high coverage
# data-sets
MAX_DEPTH = 100000

MIN_COV = 10

def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog:\n" \
        + __doc__ + "\n" \
        "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", 
                      dest="debug",
                      help="debugging")
    parser.add_option("-b", "--bam",
                      dest="fbam", # type="string|int|float"
                      help="BAM input file (- for stdin)")
    parser.add_option("-o", "--out",
                      dest="fout", # type="string|int|float"
                      default="-",
                      help="Output fasta file (- for stdout)")
    parser.add_option("-r", "--seq",
                      dest="plpref",
                      help="Optional: reconstruct only this seqeunce"
                      " chromosome/sequence")
    parser.add_option("", "--seqs-from-fa",
                      action="store_true", 
                      dest="refs_from_fa",
                      help="Optional: reconstruct all sequences mentioned in fasta file."
                      " Will otherwise use all seqeunces mentioned in BAM header")
    parser.add_option("-f", "--ref-fasta",
                      dest="ref_fa",
                      help="Optional: If given, any ambiguities (Ns and"
                      " no coverage!) in reconstructed consensus will"
                      " be resolved by using the reference seq from this"
                      " file (ignoring terminal ends)")
    default = "~"
    parser.add_option("", "--no-cov-char",
                      dest="no_cov_char",
                      default=default,
                      help="Optional: Character denoting zero coverage"
                      " (default %c). Not used with -r" % default)
    parser.add_option("-s", "--start",
                      dest="plpstart",
                      type="int",
                      help="Optional: Start reconstrunction at this position")
    parser.add_option("-e", "--end",
                      dest="plpend",
                      type="int",
                      help="Optional: End reconstruction at this position")

    return parser



def main():
    """The main function
    """

    pysam_version = tuple(int(x) for x in pysam.__version__.split('.'))
    if pysam_version < (0 , 6):
        LOG.critical("Using untested pysam version (%s)." % (
            pysam_version))

    if sys.version_info < (2 , 6) or sys.version_info >= (2 , 8):
        LOG.critical("Using untested Python version (%s)." % (
            sys.version_info))

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if len(args):
        parser.error("Unrecognized arguments found: %s." % (
            ' '.join(args)))
        sys.exit(1)

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)

    for (fin, descr) in [(opts.fbam, 'BAM input')]:
        if not fin:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(fin) and fin != "-":
            LOG.fatal("%s does not exist.\n" % fin)
            sys.exit(1)

    for (fout, descr) in [(opts.fout, 'Fasta output')]:
        if not fout:
            parser.error("%s argument missing." % descr)
            sys.exit(1)
        if os.path.exists(fout):
            LOG.fatal("%s already exists.\n" % fout)
            sys.exit(1)

    if (opts.plpstart and not opts.plpend) or (not opts.plpstart and opts.plpend):
        LOG.fatal("Need end argument if start is set, and vice versa")
        sys.exit(1)

    if opts.ref_fa:
        fastafile = pysam.Fastafile(opts.ref_fa)
    else:
        fastafile = None

    if opts.plpref and opts.refs_from_fa:
        LOG.fatal("Can't reconstruct given sequence %s and use the one from fasta at the same time" % (opts.plpref))
        sys.exit(1)
        
    if opts.plpref:
        references = [opts.plpref]
        
    elif opts.refs_from_fa:
        if not opts.ref_fa:
            LOG.fatal("If I'm supposed to reconstruct sequences from a fasta file, then please provide me with that fasta file")
            sys.exit(1)
        # FIXME SeqIO only used for id parsing. unnecesssary!
        references = [s.id for s in SeqIO.parse(open(opts.ref_fa, 'r'), 'fasta')]
                
    else:
        samfile = pysam.Samfile(opts.fbam, "rb" )
        references = samfile.references
        
    LOG.info("Will determine %d consensus seqs" % len(references))
    #LOG.debug("consensus seqs: %s" % (', '.join(references)))

    if opts.fout == '-':
        fhout = sys.stdout
    else:
        fhout = open(opts.fout, 'w')

    for plpref in references:
        samfile = pysam.Samfile(opts.fbam, "rb" )
       
        LOG.debug("starting pileup with: ref=%s, opts.plpstart=%s, opts.plpend=%s max_depth=%d" % (
            plpref, opts.plpstart, opts.plpend, MAX_DEPTH))
   
        skip_del_cols = []
        consseq = []
        colctr = -1
       
        # http://www.cgat.org/~andreas/documentation/pysam/api.html#pysam.Samfile.pileup
        # Samfile.pileup(self, reference=None, start=None, end=None, region=None, callback=None, **kwargs)
        # WARN: set bitmask as well?
        # note: pileup ignores start and end if no ref is given
        for plpcol in samfile.pileup(reference=plpref,
                                     start=opts.plpstart, end=opts.plpend,
                                     max_depth=MAX_DEPTH):
            colctr += 1
            if plpcol.pos+1 % 1000 == 0:
                LOG.debug("Working on pileup of column %d" % (plpcol.pos+1))
               
            if plpcol.pos in skip_del_cols:
                LOG.debug("Skipping col %d" % (plpcol.pos+1))
                skip_del_cols.remove(plpcol.pos)
                continue
           
            # reads overlapping given positions will also be part of
            # pileup. skip if outside given pileup positions.
            if opts.plpstart and opts.plpend:
                if plpcol.pos < opts.plpstart-1 or plpcol.pos > opts.plpend-1:
                    continue
   
            if plpcol.pos != colctr:
                diff = plpcol.pos - colctr
               
                # don't add 5' gap. because then we would have to add the
                # 3' gap as well which is tricky.
                if colctr == 0:
                    pass
               
                elif fastafile:
                    LOG.info("Missing coverage between cols %d-%d..."
                             " filling in from original ref %s" % (
                        colctr+1, plpcol.pos+1, plpref))
                    region = fastafile.fetch(plpref, colctr, plpcol.pos)
                    if len(region) == 0:
                        LOG.debug("End of seq reached at col %d", colctr+1)
                        break
                    if len(region) != diff:
                        LOG.fatal("Couldn't fetch region %s:%d-%d from %s" % (
                            plpref, colctr+1, plpcol.pos+1, opts.ref_fa))
                        sys.exit(1)
                    consseq.append(region)
                    
                else:
                    LOG.info("Missing coverage between cols %d-%d..."
                             " adding '%c' to denote coverage gap" % (
                                 colctr+1, plpcol.pos+1, opts.no_cov_char))
                    for i in range(diff):
                        consseq.append(opts.no_cov_char)
                   
                colctr += diff
                continue
   
            # DEBUG
            #for plpread in plpcol.pileups:
            #    print '\tbase in read %s = %s' % (plpread.alignment.qname, plpread.alignment.seq[plpread.qpos])
            #print 'coverage at base %s = %s' % (plpcol.pos , plpcol.n)
   
            basecounts = Counter([plpread.alignment.seq[plpread.qpos]
                                  for plpread in plpcol.pileups])
                                  
            assert plpcol.n == sum(basecounts.values()) # paranoia                                  

            if plpcol.n < MIN_COV:
                consbase = 'N'
            else:
                consbase = basecounts.most_common(1)[0][0]
            LOG.debug("pos %d: %s (%d out of %d)" % (
                plpcol.pos+1, consbase, basecounts[consbase], plpcol.n))
            #
            # NOTE: elements with equal counts are ordered arbitrarily in most_common()
            # Alternative, if we also want to handle ties:
            ### dictionary sorted by value, converted into a list
            ##basecounts = sorted(basecounts, key=basecounts.get, reverse=True)
            ##if len(basecounts)>1 and basecounts[0]==basecounts[1]:
            ##    consbase = 'N'
            ##else:
            ##    consbase = basecounts[0]
            if consbase == 'N' and fastafile:
                # FIXME never using refbase? consbase as arg to fetch
                refbase = fastafile.fetch(plpref, plpcol.pos, plpcol.pos+1, consbase)
                if len(refbase) == 0:
                    LOG.warn("Couldn't fetch region %s:%d-%d from %s." 
                             " Will keep ambigious cons base" % (
                                 plpref, plpcol.pos, plpcol.pos+1, opts.ref_fa))
   
            if fastafile:
                LOG.debug("Pos %d Ref %s Cons %s" % (
                    plpcol.pos+1, fastafile.fetch(
                        plpref, plpcol.pos, plpcol.pos+1), consbase))
                
            consseq.append(consbase)
   
   
            # length of indels and their counts
            # indel: indel length; 0 for no indel, positive for ins and negative for del
            indelcounts = Counter([plpread.indel
                                   for plpread in plpcol.pileups])
            # FIXME if we have indels with length 2 they should also count as indels of length 1
            consindel = indelcounts.most_common(1)[0]
            # consindel: first element: length; second is count
            #print 'indelcounts = %s' % (indelcounts)
            #print "consindel = %s" % str(consindel)
            
            if consindel[0] > 0:
                #print "\nreads here:"
                #for  plpread in plpcol.pileups:
                #    print "read %s  @qpos %d:%c  @pos+1 %c  indel? %s" % (
                #        plpread.alignment.seq,
                #        plpread.qpos,
                #        plpread.alignment.seq[plpread.qpos],
                #        plpread.alignment.seq[plpread.qpos+1] if plpread.indel else '/',
                #        'Y' if plpread.indel else 'N')
                #[plpread.alignment.seq[plpread.qpos+i+1 - plpread.alignment.pos]
                for i in range(consindel[0]):
                    basecounts = Counter(
                        [plpread.alignment.seq[plpread.qpos+i+1]
                         for plpread in plpcol.pileups
                         if plpread.indel > i])
                    consbase = basecounts.most_common(1)[0][0]
                    consseq.append(consbase)
                    
            elif consindel[0] < 0:
                #print "%d %s" % (plpcol.pos+1, ''.join(['-' for i in range(abs(consindel[0]))]))
                for i in range(abs(consindel[0])):
                    LOG.debug("Appending skip/del col %d" % (plpcol.pos+i+1))
                    skip_del_cols.append(plpcol.pos+i+1)
   
        # if we're not ignoring terminal ends then we would have to check
        # here is consseq has the same length as the ref seq, for which we
        # have to load it into memory?!
        
        samfile.close()
       
        if len(consseq)==0:
            LOG.warn("consensus sequence for %s is empty (no coverage)."
                     " Skipping" % (plpref))
            continue
           
        fhout.write(">%s\n" % plpref)
        fhout.write("%s\n" % ''.join(consseq))

    if fhout != sys.stdout:
        fhout.close()


if __name__ == "__main__":

    main()
    LOG.info("Successful exit")

