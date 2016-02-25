#!/usr/bin/python
usage = "multiPopVectors2OmegaScan.py [--options] vectors.txt"
description = "writes OmegaScan config files based on the multi-population vectors supplied. Assumes KW channel naming conventions"
author = "reed.essick@ligo.org"

from optparse import OptionParser

#-------------------------------------------------

parser=OptionParser(usage=usage, description=description)

#parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--frame-type", default="H1_R", type="string")
parser.add_option("", "--timeRange", default=64, type="int")

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError("Please supply exactly one input argument\n%s"%(usage))
vectors = args[0]

#-------------------------------------------------

file_obj = open(vectors, "r")
chans = file_obj.readline().strip().split()[1:] ### skip first column because that is the filename
file_obj.close()
Nchans = len(chans)
#if opts.verbose:
#    print "    found %d channels"%(Nchans)

### assume KW channel naming convention
channels = {}
for chan in chans:
    chan = chan.split("_")
    hF = int(chan[-1])
    chan = "%s:%s"%(chan[0], "_".join(chan[1:-2]))
    if channels.has_key(chan):
        channels[chan] = max(2*hF, channels[chan])
    else:
        channels[chan] = 2*hF

#-------------------------------------------------

print """# Q Scan configuration file
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
}"""%('%s', opts.frame_type, '%d', opts.timeRange)

for chan, sampleFreq in channels.items(): ### assumes KW naming conventions
    print template%(chan, sampleFreq)

