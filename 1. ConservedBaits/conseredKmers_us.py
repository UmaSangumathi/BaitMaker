#!/usr/bin/python

#==============================================================================
import sys

usage = """
usage: %s Options_fn [-help | -h | --help]

   print help message, and all pipeline steps as they would be executed,
   then exit

usage: %s Options_fn [additional options]

  run the pipeline (which consists of a number of steps (each step is
  a separate program. Some steps are in C++, some in perl, some in python.
  Their inputs and outputs are tied together via the file system.  I know
  this is a bit hoaky and archaic, but it works ;)

additional options:
  -s j         (start at step number j.  Use -help to show step numbers)
  -x h i j     (omit steps number h, i, and j)
  -e h i j     (run only steps number h, i, and j)

  Notes: - most of the "additional Options_fn" were implemented for debugging.
         - you can mix and match -s, -x, -e as desired

""" % (sys.argv[0], sys.argv[0])

#==============================================================================
if len(sys.argv) < 2:
  print usage
  sys.exit(-1)

#==============================================================================
import sys
import os
import time
import Options
from common import *

#==============================================================================
# main starts here
#==============================================================================
tt1 = time.time()

parseCmdLine(sys.argv)

sys.stderr.write("findProbes.py is starting with options file: " + Options_fn + '\n')

fp = None

# -------  delete results directory  -------
# if len(sys.argv) > 2 :
if True:
  cmd = ['rm', '-rf', Results_dir]
  msg = """
  Delete everything in the result directory (if it exists).
  This means you'll be starting with a blank slate, which is the safest
  way to proceeed.  However, if you're dealing with a large dataset
  (such as all_bacteria), you'll be computing the suffix array from
  scratch, which could take several hours.

  Note: this step will not appear in the logfile
  """
  runCmd(cmd, msg, fp)

# ------- create result and vtest directories -------
os.system('pwd')

cmd = ['mkdir', '-p', Vtest_dir]
msg = """
  Create result and vtest directories

  Note: this step will not appear in the logfile
"""
runCmd(cmd, msg, fp)
print "made dir: ", Vtest_dir

# ------- open logfile -------
# Write_cmd = True
logfile = None
if Write_logfile == "1":
  fp = open(Results_dir + '/logfile', 'w')
  print 'opened logfile: ' + Results_dir + '/logfile'
else:
  print 'logfile disabled'

# ------- copy Options_fn file -------
copyOptionsFile()
cmd = ['cp', '-f', Options_fn, Results_dir]

# ------- compute suffix array -------
cmd = [Exec_dir + '/sa', Options_fn]
msg = 'Compute suffix array'
t = runCmd(cmd, msg, fp)

# ------- read the sequence count -------
cmd = ['ls', Results_dir + '/seq_count']
msg = 'Test that seq_count file exists'
runCmd(cmd, msg, fp)

a = open(Results_dir + '/seq_count').readlines()
Seq_count = a[0][:-1]

# ------- copy query file -------
cmd = ['cp', Results_dir + '/queries', Results_dir + '/queries_org']
msg = 'Copy queries to queries_org'
runCmd(cmd, msg, fp)

# ------- run mkvtree -------
cmd = [Exec_dir + '/run_mkvtree.pl', Seq_fn, Vtest_dir]
msg = 'Run mkvtree'
runCmd(cmd, msg, fp)

# ------- run vmatch for query file -------
cmd = [Exec_dir + '/run_vmatch.pl', Options_fn, Vtest_dir + '/vtest', Results_dir + '/queries']
msg = 'Run vmatch'
runCmd(cmd, msg, fp)

# ------- parse vmatch output -------
cmd = [Exec_dir + '/pvmo.py', Results_dir + '/queries_vm']
msg = 'Parse vmatch output'
runCmd(cmd, msg, fp)

if not Help:
  ss_count_1 = getCount(Results_dir + '/queries_vm_supersets')

# ------- copy superset file -------
cmd = ['cp', Results_dir + '/queries_vm_supersets', Results_dir + '/queries_vm_supersets_org']
msg = 'Copy queries_vm_supersets to queries_vm_supersets_org'
runCmd(cmd, msg, fp)

# ------- merge the query superset file -------
# this reads queries_vm_supersets_fwd/rc, filters, and writes
# results to the same filename
cmd = [Exec_dir + '/mss', Options_fn, 'queries_vm_supersets', 'x']
msg = 'Merge the queries superset file'
runCmd(cmd, msg, fp)

if not Help:
  ss_count_3 = getCount(Results_dir + '/queries_vm_supersets')


if Help == False:
  print "query superset count after parsing vmatch:    ", ss_count_1
  print "query superset count after merging:           ", ss_count_3
