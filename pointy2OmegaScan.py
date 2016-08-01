#!/usr/bin/python
usage = "pointy2OmegaScan.py [--options] gps,pointy.out gps,pointy.out ..."
description = "builds an OmegaScan config file based off the channels with pvalues below --pvalueThr"
author = "reed.essick@ligo.org"

import numpy as np

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option('-p', '--pvalueThr', default=1e-3, type='float')

parser.add_option("", "--frame-type", default="H1_R", type="string")
parser.add_option("", "--timeRange", default=64, type="int")

parser.add_option("", "--freq-map", default=None, type="string", help="the output of FrChannels, used to map channel names to sample frequencies")

parser.add_option("", "--gwchan", default="H1:CAL-DELTAL_EXTERNAL_DQ", type="string")

parser.add_option("", "--output-dir", default=".", type="string")

parser.add_option("", "--condor", default=False, action="store_true", help="write a condor_sub file instead of a shell script")
parser.add_option("", "--accounting-group", default="ligo.dev.o2.detchar.explore.test", type="string")
parser.add_option("", "--accounting-group-user", default="reed.essick", type="string")
parser.add_option("", "--request-memory", default=2000000, type="int", help="measured in kB")

opts, args = parser.parse_args()

if len(args)<1:
    raise ValueError("Please supply exactly at least one input argument\n%s"%(usage))
pointys = []
gps = []
for arg in args:
    g, p = arg.split(',')
    point.append( p )
    gps.append( float(g) )

if opts.freq_map==None:
    opts.freq_map = raw_input("--freq-map=")

if not os.path.exists(opts.output_dir):
    os.makedirs( opts.output_dir )

ifo = opts.frame_type[0]

#-------------------------------------------------

### set up commands and such

gwdf_cmd = "gw_data_find -o %s --type %s"%(ifo, opts.frame_type) + " -s %d -e %d -u file"
os_cmd = "/home/omega/opt/omega/bin/wpipeline scan %.9f -r -c %s -o %s -f %s"

header = """# Q Scan configuration file
# Automatically generated with wconfigure.sh
# by user bhughey on 2009-07-09 10:33:18 PDT
# from sample frame files:
#   /archive/frames/S6/L1/LHO/H-H1_RDS_R_L1-9311/H-H1_RDS_R_L1-931194752-64.gwf

[Context,Context]

[Parameters,Parameter Estimation]

[Notes,Notes]

[Aux Channels,Identified interesting Aux channels]
"""

template = """{
  channelName:                 '%s'
  frameType:                   '%s'
  sampleFrequency:             %s
  searchTimeRange:             %d
  searchFrequencyRange:        [0 Inf]
  searchQRange:                [4 64]
  searchMaximumEnergyLoss:     0.2
  whiteNoiseFalseRate:         1e-3
  searchWindowDuration:        0.5
  plotTimeRanges:              [0.1 1 4 16]
  plotFrequencyRange:          []
  plotNormalizedEnergyRange:   [0 25.5]
  alwaysPlotFlag:              1
}"""%('%s', opts.frame_type, '%s', opts.timeRange)

if opts.condor:
    cmd_file = "%s/run_Qscan.sub"%(opts.output_dir)
    cmd_obj = open(cmd_file, "w")
    print >> cmd_obj, """universe = vanilla
executable = /home/omega/opt/omega/bin/wpipeline
getenv = True
accounting_group = %s
accounting_group_user = %s
log = %s/Qscan.log
error = %s/Qscan-$(cluster)-$(process).err
output = %s/Qscan-$(cluster)-$(process).out
request_memory = %d KB
notification = never"""%(opts.accounting_group, opts.accounting_group_user, opts.output_dir, opts.output_dir, opts.output_dir, opts.request_memory)

else:
    cmd_file = "%s/run_Qscan.sh"%(opts.output_dir)
    cmd_obj = open(cmd_file, "w")

#-------------------------------------------------

if opts.verbose:
    print "reading in sample frequencies from :"+opts.freq_map

file_obj = open(opts.freq_map, "r")
freq_map = dict( [l.strip().split() for l in file_obj] )
file_obj.close()

### extract relevant channels from pointy.out
for t, pointy in zip(gps, pointys):

    if opts.verbose:
        print "writing Qscan config files for: %d -> %s"%(t, pointy)

    outdir = "%s/%.6f"%(opts.output_dir, gps)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    chans = []

    if opts.verbose:
        print "processing : %s"%(pointy)

    file_obj = open(pointy, "r")
    lines = file_obj.readlines()
    file_obj.close()

    nlines = len(lines)
    ind = 0
    while ind < nlines:
        if "channel" in lines[ind]:
            chan = lines[ind].strip().split("=")[-1]
            nind = ind + 4
            if "pvalue" not in lines[nind]:
                nind += 4
            pvalue = float(lines[nind].strip().split("=")[-1])

            if pvalue <= opts.pvalueThr: ### only include if pvalue is small enough
                chans.append( chan )

        ind += 1

    if opts.verbose:
        print "    found %d relevant KW channels"%(len(chas))

    ### convert to raw channel names
    channels = set()
    for chan in chans:
        chan = chan.split("_")
        chan = "%s:%s"%(chan[0], "_".join(chan[1:-2]))
        channels.add( chan )

    if opts.verbose:
        print "    corresponds to %d raw channels"%(len(channels))

    ### set up channel frequency stuff
    channels = dict( (chan, freq_map[chan]) for chan in channels )

    raise NotImplementedError("WRITE ACTUAL CONFIG FILE")

    ### write config file
    conf_file = "%s/Qscan.cnf"%(outdir)
    if opts.verbose:
        print "        "+conf_file
    conf_obj = open(conf_file, "w")

    print >> conf_obj, header
    print >> conf_obj, template%(opts.gwchan, freq_map[opts.gwchan])

    for chan in sorted(channels.keys()): ### assumes KW naming conventions
        print >> conf_obj, template%(chan, channels[chan])

    conf_obj.close()

    ### set up command
    this_cmd = gwdf_cmd%(int(gps), int(gps)+1)
    if opts.verbose:
        print "        "+this_cmd
    frame = sp.Popen( this_cmd.split(), stdout=sp.PIPE).communicate()[0].split()[0]
    directory = os.path.dirname( frame.replace("file://localhost","") )

    that_cmd = os_cmd%(gps, conf_file, outdir, directory)
    if opts.verbose:
        print "        "+that_cmd
    if opts.condor:
        print >> cmd_obj, "arguments = \" %s \"\nqueue 1"%(" ".join(that_cmd.split()[1:]))
    else:
        print >> cmd_obj, that_cmd

cmd_obj.close()
if opts.verbose:
    if opts.condor:
        print "now run :\ncondor_submit %s"%(cmd_file)
    else:
        print "now run :\n%s"%(cmd_file)
