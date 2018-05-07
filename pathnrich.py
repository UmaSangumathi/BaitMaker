#!/usr/bin/python

import argparse
import sys
import os.path

def check_file(parser, arg):
    if not ".db" in arg:
        parser.error("db: The file '%s' is not a db file" % arg)
    #elif not os.path.exists(arg):
    #    parser.error("db: The file '%s' does not exist" % arg)
    else:
        return arg


def check_len(value):
    ivalue = int(value)
    if 80 <= ivalue <= 150:
        return ivalue
    else:
        raise argparse.ArgumentTypeError("--len %s is out of range, 80-150" % value)
    

def check_cov(value):
    ivalue = int(value)
    if 80 <= ivalue <= 60000:
        return ivalue
    else:
        raise argparse.ArgumentTypeError("--cov %s is out of range, 500-60000" % value)
    

class pathnrich(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='pathnrich to generate baits for NGS Library enrichment preparation',
            usage='''pathnrich.py <command> [<args>]

Commands for pathnrich are:
   genbaits     Generate baits from sequence database
   minbaits     Minimise the number of baits that were generated using genbaits
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print 'Unrecognized command'
            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()
        

    def genbaits(self):
        parser = argparse.ArgumentParser(
            description='Generate baits from sequence database',
            usage='''pathnrich.py genbaits [sp] [op] [db] [<args>]
            ''')
        # prefixing the argument with -- means it's optional

        parser.add_argument('sp',  # '--species',
            type = str, # required = True,
            help = "specify the species name for bait generation that is also specified in the dbase file")
        parser.add_argument('op',  # '--opath',
            type = str, # required = True,
            help = "specify output path for bait generation")
        parser.add_argument('db',  # '--dbase', 
            type = lambda x: check_file(parser,x), # required = True,
            help = "specify database file '.db' for bait generation")
        parser.add_argument('--len', '-l',
            default = 120, type = check_len,
            help = "specify length of bait, 80-150, default: 120")
        parser.add_argument('--cov', '-c',
            default = 500, type = check_cov,
            help = "specify the distance between the start of a bait to the next bait, 500-60000, default: 500")
        parser.add_argument('--disRC', '-d',
            action='store_true',
            help = "disable reverse complement function for '-'gRNA viruses")
        parser.add_argument('--verbose', '-v',
            action='store_true',
            help = "enable verbosity")

        # now that we're inside a subcommand, ignore the first
        # TWO argvs, ie the command (pathnrich.py) and the subcommand (genbaits)
        args = parser.parse_args(sys.argv[2:])
        dbname = args.db
        spname = args.sp
        outpath = args.op
        len = str(args.len)
        cov = str(args.cov)
        disRC = str(args.disRC)

        print ''
        print '*** Running pathnrich genbaits ***'
        print 'species =', spname
        print 'output folder =', outpath
        print 'database =', dbname
        print 'len =', len
        print 'cov =', cov
        print 'disRevCom =', disRC
        print ''

        command = "python ./src/genbaits/runVirusDNA_NS.py " + spname + " " + outpath + " " + dbname + " " + len + " " + cov + " " + disRC
        print command
        print ''
        os.system(command)
        print ''
        

    def minbaits(self):
        parser = argparse.ArgumentParser(
            description='Minimise the number of baits that were generated using genbaits',
            usage='''pathnrich.py minbaits [<args>]
            ''')
        # NOT prefixing the argument with -- means it's not optional
        parser.add_argument('repository')
        args = parser.parse_args(sys.argv[2:])
        print 'Running pathnrich minbaits'
        print args

  
if __name__ == '__main__':
    pathnrich()
