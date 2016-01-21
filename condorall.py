description = """an all-in-one pointy analysis tool"""
usage = """doitall.py [--options] gps gps gps..."""

import os
import subprocess as sp

from optparse import OptionParser

#=================================================

sub_header = """universe = %s
executable = %s/doitall.py
arguments = " %s "
getenv = True
log = /usr1/%s/pointy-Poisson_doitall-$(Cluster)-$(Process).log
error = %s/$(JobID)/logs/doitall.err
output = %s/$(JobID)/logs/doitall.out
notification = never
accounting_group = %s
accounting_group_user = %s
queue 1"""

#=================================================

parser = OptionParser(description=description, usage=usage)

#====

parser.add_option("-v", "--verbose", default=False, action="store_true")

#==== options about data ranges

parser.add_option("-o", "--observatory", default=False, type="string")

parser.add_option("-t", "--frame-type", default="C", type="string", help="either \"C\" (default) or \"R\"")

parser.add_option("", "--whitelist", default=False, type="string", help="a file containing a list of channel names defining which channels should be analyzed. Only the channels in this list will be analyzed and we require an exact match. Channel names should be as they come out of the frames (less IFO prefix). Default is False, so all channels present will be allowed.")
parser.add_option("", "--error-if-absent", default=False, action="store_true", help="raise an error if one of the channels in --whitelist is not found in the frames")

parser.add_option("", "--blacklist", default=False, type="string", help="a file containing a list of channel names defining which channels should not be analyzed. Channel names should be as they come out of the frames (less IFO prefix). Default is False, so all channels present will be allowed.")

#==== options about KW processing and frame width

parser.add_option("-s", "--kwstride", default=32, type="int", help="desired stride for KW analysis. Should be a power of 2")
parser.add_option("-S", "--framestride", default=64, type="int", help="duration of frames. Should be a power of 2")

#==== options about pointy analysis

parser.add_option("", "--pointybin", default="/home/idq/fishyStatistics/pointy-Poisson/", type="string", help="the directory which contains the pointy executables and sub files" )

parser.add_option("-w", "--coinc-window", default=1.0, type="float")
parser.add_option("-W", "--signif-window", default=5000.0, type="float")
parser.add_option("-e", "--exclude", default=0.0, type="float", help="passed on to pointy commands")

parser.add_option("", "--signifmin", default=[], action="append", type="float")
parser.add_option("", "--signifmax", default=1e6, type="float")
parser.add_option("", "--fmin", default=0, type="float")
parser.add_option("", "--fmax", default=16384, type="float")
parser.add_option("", "--durmin", default=0.0, type="float")
parser.add_option("", "--durmax", default=100.0, type="float")

#==== options about output and post-processing

parser.add_option("-p", "--pvalue", default=1.0, type="float")
parser.add_option("-u", "--unsafe", default=None, type="string")

parser.add_option("-O", "--output-dir", default=".", type="string")

#==== options specifically about condor processing

parser.add_option("", "--condor_submit_dag", default=False, action="store_true", help="automatically submit the dag")

parser.add_option("", "--universe", default="vanilla", type="string" )
parser.add_option("", "--accounting_group", default=False, type="string" )
parser.add_option("", "--accounting_group_user", default=False, type="string" )

parser.add_option("", "--retry", default=0, type="int")

#====

opts, args = parser.parse_args()

if not len(args):
    raise ValueError("please supply at least one gps time as an argument")
args = [float(arg) for arg in args]

if not opts.signifmin:
    opts.signifmin.append( 0.0 )

if not opts.observatory:
    opts.obervatory = raw_input("observatory = ")

if not opts.accounting_group:
    opts.accounting_group = raw_input("accounting_group = ")
if not opts.accounting_group_user:
    opts.accounting_group_user = raw_input("accounting_group_user = ")

cwd = os.getcwd()
if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

#=================================================

### figure out whoami
whoami = sp.Popen(["whoami"], stdout=sp.PIPE).communicate()[0].strip()

#=================================================

### we need a sub file for doitall.py
subfile = "%s/pointy-Poisson_%s.sub"%(opts.output_dir, opts.frame_type)

subfile_args = " $(gpsTarget) -o %s -t %s -s %d -S %d --pointybin %s -w %.6f -W %.6f -e %.6f --fmin %.6f --fmax %.6f --durmin %.6f --durmax %.6f -p %.6f -O %s/$(JobID)/ --signifmax %.3e --signifmin %s"%(opts.observatory, opts.frame_type, opts.kwstride, opts.framestride, opts.pointybin, opts.coinc_window, opts.signif_window, opts.exclude, opts.fmin, opts.fmax, opts.durmin, opts.durmax, opts.pvalue, opts.output_dir, opts.signifmax, " --signifmin ".join("%.6f"%signifmin for signifmin in opts.signifmin))
if opts.verbose:
    subfile_args += " -v"
if opts.unsafe:
    subfile_args += " -u %s"%(opts.unsafe)
if opts.whitelist:
    subfile_args += " --whitelist %s"%(opts.whitelist)
if opts.error_if_absent:
    subfile_args += " --error-if-absent"
if opts.blacklist:
    subfile_args += " --blacklist %s"%(opts.blacklist)

if opts.verbose:
    print "writing %s"%(subfile)
subfile_obj = open(subfile, "w")
subfile_obj.write( sub_header%(opts.universe, opts.pointybin, subfile_args, whoami, opts.output_dir, opts.output_dir, opts.accounting_group, opts.accounting_group_user) )
subfile_obj.close()

#=================================================

dagfile = "%s/pointy-Poisson_%s.dag"%(opts.output_dir, opts.frame_type)
if opts.verbose:
    print "writing %s"%(dagfile)
dagfile_obj = open(dagfile, "w")

N = len(args)
for ind, gps in enumerate(args):

    jobID = "%s-%d-%d"%(opts.observatory, int(gps), ind)

    if opts.verbose: 
        print "%d / %d : jobID -> %s"%(ind+1, N, jobID)

    output_dir = "%s/%s/logs/"%(opts.output_dir, jobID)
    if not os.path.exists(output_dir):
        os.makedirs( output_dir )

    print >> dagfile_obj, "JOB   %s %s"%(jobID, subfile)
    print >> dagfile_obj, "VARS  %s JobID=\"%s\" gpsTarget=\"%.9f\""%(jobID, jobID, gps)
    print >> dagfile_obj, "RETRY %s %d"%(jobID, opts.retry)

dagfile_obj.close()

if opts.verbose:
    print "%s now ready for submission"%(dagfile)

if opts.condor_submit_dag:
    cmd = "condor_submit_dag %s"%(dagfile)
    if opts.verbose:
        print "    "+cmd

    sp.Popen(cmd.split()).wait()
