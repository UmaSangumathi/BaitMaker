#!/usr/bin/env python
"""Create a coverage plot from a BAM (or RazerS) file
"""

# FIXME this needs a proper overhaul


#--- standard library imports
#
import sys
import logging
# optparse deprecated from Python 2.7 on
from optparse import OptionParser
import gzip
import os
import subprocess
import tempfile


#--- third-party imports
#
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as pyplot
import numpy

#--- project specific imports
#
# /


# max sequence len
MAX_LEN = 5e6


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



def genome_from_bam(fbam, fh_genome, samtools="samtools"):
    """
    Create A bedtools 'genome' file

    Arguments:
    - `fbam`: bam file to create genome file from
    - `fh_genome`: bedtools genome file handle
    - `samtool`: path to samtools binary
    """

    cmd = "%s view -H %s" % (samtools, fbam)
    process = subprocess.Popen(cmd.split(),
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    (stdoutdata, stderrdata) =  process.communicate()

    retcode = process.returncode
    if retcode != 0:
        raise OSError("Called command exited with error code '%d'." \
                      " Command was '%s'. stderr was: '%s'" % (
                          retcode, cmd, stderrdata))
    for line in str.splitlines(stderrdata):
        if len(line.strip()) == 0:
            continue
        LOG.warn("Got following stderr message from samtools: %s" % line)


    # local requirement would be == 1
    #if stdoutdata.count("@SQ") != 1:
    #     raise ValueError(
    #         "Did not find exactly one @SQ line in BAM file '%s'"  % (
    #             fbam))

    for line in str.splitlines(stdoutdata):
        if line.startswith("@SQ"):
            line_split = line.split('\t')
            chrname = line_split[1]
            assert chrname.startswith("SN:")
            chrname = chrname[3:]

            chrlen = line_split[2]
            assert chrlen.startswith("LN:")
            chrlen = int(chrlen[3:])

            fh_genome.write("%s\t%d\n" % (chrname, chrlen))



def parse_coverage_bam(fbam, genomecoveragebed="genomeCoverageBed"):
    """

    Arguments:
    - `fbam`: file to read coverage from. Needs samtools and bedtools (genomeCoverageBed)
    - `genomecoveragebed`: path to genomeCoverageBed binary
    """

    try:
        cmd = [genomecoveragebed]
        process = subprocess.Popen(cmd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        (stdoutdata, stderrdata) =  process.communicate()
    except OSError:
        LOG.fatal("Can't run %s. Please make sure it's in your path" % cmd[0])
        sys.exit(1)

    # forward and reverse coverage
    # int is essential; see note below
    fw_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)
    rv_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)


    file_genome = tempfile.NamedTemporaryFile(delete=False)
    genome_from_bam(fbam, file_genome)
    file_genome.close()

    basic_cmd = "%s -ibam %s -g %s -d" % (
        genomecoveragebed, fbam, file_genome.name)

    for strand in ["+", "-"]:
        cmd = "%s -strand %s" % (basic_cmd, strand)
        process = subprocess.Popen(cmd.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        (stdoutdata, stderrdata) =  process.communicate()
        retcode = process.returncode
        if retcode != 0:
            raise OSError("Called command exited with error code '%d'." \
                          " Command was '%s'. stderr was: '%s'" % (
                              retcode, cmd, stderrdata))
        for line in str.splitlines(stderrdata):
            if len(line.strip()) == 0:
                continue
            LOG.warn("Got following stderr message from genomeCoverageBed: %s" % (line))
        for line in str.splitlines(stdoutdata):
            if len(line) == 0:
                continue
            (chrom, pos, cov) = line.split('\t')
            pos = int(float(pos))-1 # we use zero offset
            cov = int(float(cov))

            assert pos < MAX_LEN

            if strand == '+':
                fw_cov[pos] = cov
            elif strand == '-':
                rv_cov[pos] = cov

    os.unlink(file_genome.name)

    try:
        max_fw_pos = max(numpy.where(fw_cov>0)[0])
    except ValueError:# no coverage!
        max_fw_pos = 0
    try:
        max_rv_pos = max(numpy.where(rv_cov>0)[0])
    except ValueError:# no coverage!
        max_rv_pos = 0

    max_pos = max(max_fw_pos, max_rv_pos)

    return (fw_cov[:max_pos+1],
            rv_cov[:max_pos+1])


def parse_coverage_plp(fhandle):
    """
    """

    cov_at_pos = dict()
    for line in fhandle:
        (pos, cov) = line.rstrip().split('\t')
        (pos, cov) = (int(pos)-1, int(cov))
        cov_at_pos[pos] = cov

    cov_arr = numpy.zeros(max(cov_at_pos.keys())+1, dtype=numpy.int32)

    for (k, v) in cov_at_pos.items():
        cov_arr[k] = v
    return cov_arr[:MAX_LEN+1]



def parse_coverage_razers(fhandle):
    """
    Returns forward and reverse coverage as numpy array.

    Arguments:
    - `fhandle`: file handle to read razers mapping from
    """

    # forward and reverse coverage
    # int is essential; see note below
    fw_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)
    rv_cov = numpy.zeros((MAX_LEN), dtype=numpy.int32)

    readlen = None
    max_pos_seen = 0

    for line in fhandle:
        if line.startswith('#'):
            continue
        #print "DEBUG %s" % line,

        # taken from razers.py
        line_split = line.split('\t')
        assert len(line_split) == 8 or len(line_split) == 11, (
            "Was expecting 8 or 11 elements but got %d in line: %s" % (
                len(line_split), line_split))
        #if len(line_split) > 8:
        #    raise ValueError, (
        #        "Not tested with paired end reads")

        (rname,
         rbegin,
         rend,
         gstrand,
         gname,
         gbegin,
         gend,
         percid) = line_split[:8]
        # don't need pairing information

        # razers position format: zero-based & half open. as python
        # slices/ranges
        gbegin = int(gbegin)
        gend = int(gend)

        rbegin = int(rbegin)
        rend = int(rend)
        if readlen:
            assert rend-rbegin == readlen, (
                "Can't handle reads of different lengths")
        else:
            readlen = rend-rbegin

        if gend > max_pos_seen:
            max_pos_seen = gend

        # this will work even with indels. however, the indel position
        # will also be treated as covered! to prevent this, we'd have
        # to parse the sequences as well which is not implemented at
        # the moment
        assert gend-gbegin == readlen, (
            "Can't handle indels")

        if gstrand == 'F':
            fw_cov[gbegin:gend] += 1
        elif gstrand == 'R':
            rv_cov[gbegin:gend] += 1
        else:
            raise ValueError, "Unknown strand orientation '%s'" % gstrand

    #return (numpy.trim_zeros(fw_cov, trim='b'),
    #        numpy.trim_zeros(rv_cov, trim='b'),
    #        readlen)
    return (fw_cov[0:max_pos_seen], rv_cov[0:max_pos_seen])



def cmdline_parser():
    """
    creates an OptionParser instance
    """

    # http://docs.python.org/library/optparse.html
    usage = "%prog:" + __doc__ + "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="be verbose")
    parser.add_option("", "--debug",
                      action="store_true", dest="debug",
                      help="debugging")

    parser.add_option("-b", "--bam",
                      dest="fbam", # type="string|int|float"
                      help="BAM input file (needs bedtools installed!; conflicts with --razers and --plp)")
    parser.add_option("-r", "--razers",
                      dest="frazers", # type="string|int|float"
                      help="razers input file (gzip supported; use '-' for stdin; conflicts with --bam amd --plp)")
    parser.add_option("-s", "--plp",
                      dest="fplp", # type="string|int|float"
                      help="samtools pileup | cut -f 2,4 (gzip supported; use '-' for stdin; conflicts with --bam and --razers)")
    #parser.add_option("-o", "--output",
    #                  dest="fout", # type="string|int|float"
    #                   help="output file")
    parser.add_option("-o", "--plot",
                      dest="fplot", # type="string|int|float"
                      help="plot output file")
    parser.add_option("-l", "--log",
                      dest="flog", # type="string|int|float"
                       help="print stats to this file (- for stdout)")
    parser.add_option("", "--force",
                      dest="force_overwrite",
                      action="store_true",
                      help="Force overwriting of files")
    parser.add_option("", "--start",
                      dest="startpos", type="int",
                      help="optional: position to start counting at (unit-offset; inclusive)")
    parser.add_option("", "--end",
                      dest="endpos", type="int",
                      help="optional: position to end counting at (unit-offset; inclusive)")
    parser.add_option("", "--plot-title",
                      dest="plottitle", # type="string|int|float"
                       help="optional: plot title (only useful with --plot)")
    parser.add_option("", "--no-legend",
                      dest="nolegend", action="store_true",
                      help="optional: don't print legend (only useful with --plot)")
    parser.add_option("", "--no-log",
                      dest="nolog", action="store_true",
                      help="optional: don't use log for coverage (only useful with --plot)")
    parser.add_option("", "--plot-xmin",
                      dest="plotxmin", type="int",
                       help="x-min (only useful with --plot)")
    parser.add_option("", "--plot-xmax",
                      dest="plotxmax", type="int",
                       help="optional: x-max (only useful with --plot)")
    parser.add_option("", "--plot-ymax",
                      dest="plotymax", type="int",
                       help="optional: y-max (only useful with --plot)")

    return parser



def main():
    """
    The main function
    """

    parser = cmdline_parser()
    (opts, args) = parser.parse_args()

    if opts.verbose:
        LOG.setLevel(logging.INFO)
    if opts.debug:
        LOG.setLevel(logging.DEBUG)


    if not any([opts.frazers, opts.fbam, opts.fplp]):
        parser.error("Need either BAM, pileup or RazerS input.")
        sys.exit(1)
    if sum([1 for f in [opts.frazers, opts.fbam, opts.fplp] if f]) != 1:
        parser.error("Can only use one: BAM, pileup or RazerS as input.")
        sys.exit(1)

    if not opts.fplot:
        LOG.fatal(
            "Missing plot output argument")
        sys.exit(1)
    if os.path.exists(opts.fplot) and not opts.force_overwrite:
        LOG.fatal(
            "Refusing to overwrite already existing file '%s'\n" % opts.fplot)
        sys.exit(1)

    if opts.flog and opts.flog != "-":
        if os.path.exists(opts.flog) and not opts.force_overwrite:
            LOG.fatal(
                "Refusing to overwrite already existing file '%s'\n" % (opts.flog))
            sys.exit(1)
    if not opts.flog or opts.flog == '-':
        fh_log = sys.stdout
    else:
        fh_log = open(opts.flog, 'w')

    
    combined_cov = fw_cov = rv_cov = None
    if opts.frazers:
        if opts.frazers == "-":
            fhandle = sys.stdin
        elif opts.frazers[-3:] == ".gz":
            fhandle = gzip.open(opts.frazers, 'r')
        else:
            fhandle = open(opts.frazers, 'r')
        (fw_cov, rv_cov) = parse_coverage_razers(fhandle)
    
    elif opts.fbam:
        (fw_cov, rv_cov) = parse_coverage_bam(opts.fbam)
        fhandle = None
        
    elif opts.fplp:
        if opts.fplp == "-":
            fhandle = sys.stdin
        elif opts.fplp[-3:] == ".gz":
            fhandle = gzip.open(opts.fplp, 'r')
        else:
            fhandle = open(opts.fplp, 'r')

        # only get combined_cov from plp
        combined_cov = parse_coverage_plp(fhandle)
        fw_cov = rv_cov = None
        if opts.startpos:
            combined_range_start = opts.startpos-1
        else:
            try:
                combined_range_start = min(numpy.where(combined_cov>0)[0])
            except ValueError: # no coverage
                combined_range_start = 0
        if opts.endpos:
            combined_range_end = opts.endpos-1
        else:
            try:
                combined_range_end = max(numpy.where(combined_cov>0)[0])
            except ValueError: # no coverage
                combined_range_end = 0


    else:
        LOG.fatal("Internal error: got neither BAM nor RazerS input")
        sys.exit(1)

    if fhandle:
        if fhandle != sys.stdout:
            fhandle.close()

    # print coverage report

    fh_log.write("Coverage (ignoring leading and trailing zeros): Min Max Avg\n")

    if not opts.startpos:
        try:
            fw_range_start = min(numpy.where(fw_cov>0)[0])
        except ValueError: # no coverage
            fw_range_start = 0
        #fw_range_start = 0
        #while fw_cov[fw_range_start] == 0:
        #    fw_range_start += 1

        try:
            rv_range_start = min(numpy.where(rv_cov>0)[0])
        except ValueError: # no coverage
            rv_range_start = 0
        #rv_range_start = 0
        #while rv_cov[rv_range_start] == 0:
        #    rv_range_start += 1
    else:
        fw_range_start = opts.startpos-1
        rv_range_start = opts.startpos-1


    if not opts.endpos:
        #fw_range_end = len(fw_cov)
        #while fw_cov[fw_range_end-1] == 0:
        #    fw_range_end -= 1
        try:
            fw_range_end = max(numpy.where(fw_cov>0)[0])
        except ValueError: # no coverage
            fw_range_end = 0

        #rv_range_end = len(rv_cov)-1
        #while rv_cov[rv_range_end] == 0:
        #    rv_range_end -= 1
        try:
            rv_range_end = max(numpy.where(rv_cov>0)[0])
        except ValueError: # no coverage
            rv_range_end = 0

    else:
        fw_range_end = opts.endpos
        rv_range_end = opts.endpos


    try:
        fw_min = numpy.min(fw_cov[fw_range_start:fw_range_end])
    except (ValueError, TypeError):
        fw_min = -1

    try:
        fw_max = numpy.max(fw_cov[fw_range_start:fw_range_end])
    except (ValueError, TypeError):
        fw_max = -1

    try:
        fw_mean = numpy.mean(fw_cov[fw_range_start:fw_range_end])
    except (ValueError, TypeError):
        fw_mean = -1

    try:
        rv_min = numpy.min(rv_cov[rv_range_start:rv_range_end])
    except (ValueError, TypeError):
        rv_min = -1

    try:
        rv_max = numpy.max(rv_cov[rv_range_start:rv_range_end])
    except (ValueError, TypeError):
        rv_max = -1

    try:
        rv_mean = numpy.mean(rv_cov[rv_range_start:rv_range_end])
    except (ValueError, TypeError):
        rv_mean = -1


    fh_log.write("Forward[%d:%d]: %d %d %f\n" % (
        fw_range_start, fw_range_end,
        fw_min, fw_max, fw_mean))


    fh_log.write("Reverse[%d:%d]: %d %d %f\n" % (
        rv_range_start, rv_range_end,
        rv_min, rv_max, rv_mean))

    if combined_cov == None:
        combined_cov = fw_cov+rv_cov
        combined_range_start = numpy.min([fw_range_start, rv_range_start])
        combined_range_end = numpy.max([fw_range_end, rv_range_end])
    #import pdb; pdb.set_trace()
    # FIXME will fail if empty:
    #  e.g.
    # /mnt/software/lib/python2.7/site-packages/numpy/core/fromnumeric.py:2374: RuntimeWarning: invalid value encountered in double_scalars
    #     return mean(axis, dtype, out)
    # Traceback (most recent call last):
    # File "/mnt/software/bin/coverage_plot.py", line 536, in <module>
    # main()
    # File "/mnt/software/bin/coverage_plot.py", line 480, in main
    # numpy.min(combined_cov[combined_range_start:combined_range_end]),
    # File "/mnt/software/lib/python2.7/site-packages/numpy/core/fromnumeric.py", line 1895, in amin
    #     return amin(axis, out)
    # ValueError: zero-size array to minimum.reduce without identity
    # make: *** [Sample_PHH123_s1_bwa-uniq-mapping-success.txt] Error 1

    fh_log.write("Combined[%d:%d]: %d %d %f\n" % (
        combined_range_start, combined_range_end,
        numpy.min(combined_cov[combined_range_start:combined_range_end]),
        numpy.max(combined_cov[combined_range_start:combined_range_end]),
        numpy.mean(combined_cov[combined_range_start:combined_range_end])))

    if fh_log != sys.stdout:
        fh_log.close()
    #import pdb; pdb.set_trace()


    # plotting
    #
    fig = pyplot.figure()
    ax1 = pyplot.subplot(111)
    if not opts.nolog:
        ax1.set_yscale('log')
    if fw_cov != None:
        fw_handle = ax1.plot(fw_cov, color='green', label = 'Forward')
    if rv_cov != None:
        rv_handle = ax1.plot(rv_cov, color='red', label = 'Reverse')
    combined_handle = ax1.plot(combined_cov, color='black', label = 'Combined')
    #ax1.axis([1000, 5000, 0, 5000])
    (xmin, xmax, ymin, ymax) = ax1.axis()
    if opts.plotxmin:
        xmin = opts.plotxmin
    if opts.plotxmax:
        xmax = opts.plotxmax
    if opts.plotymax:
        ymax = opts.plotymax
    ax1.axis([xmin, xmax, ymin, ymax])
    if not opts.nolegend:
        # handles, labels = ax1.get_legend_handles_labels()
        # ax1.legend(handles[::-1], labels[::-1])
        #
        # or:
        #
        # put legend outside of plot/axes.
        # see http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
        # First shink current axis by 20% and then place the legend outside axis.
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        leg = ax1.legend(loc='upper left', ncol=1, bbox_to_anchor=(1, 1))
        # set legend fontsize (needs labels)
        for t in leg.get_texts():
            t.set_fontsize('small')
    ax1.set_ylabel("Coverage")
    ax1.set_xlabel("Position")
    ax1.grid()
    if opts.plottitle:
        ax1.set_title(opts.plottitle)
    # write to file
    #
    fileext = os.path.splitext(opts.fplot)[1]
    pyplot.savefig(opts.fplot, format=fileext[1:])
    #sys.stderr.write("Pseudointeractivity...\n")
    #import pdb; pdb.set_trace();


if __name__ == "__main__":
    main()
    LOG.info("successful exit")
