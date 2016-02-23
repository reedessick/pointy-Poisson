#!/usr/bin/python
usage = "stripRawUnsafe.py unsafe.txt interesting.txt"
description = "reads in a list of unsafe channels from unsafe.txt. If these are not \"raw\" channel names, it converts them to that form. I then reads in a channel list from interesting.txt and performs a filter based on the unsafe channels. Channels not flagged as unsafe are printed to stdout while channels flagged as unsafe are printed to stderr"
author = "reed.essick@ligo.org"

import sys
from collections import defaultdict
from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

opts, args = parser.parse_args()

if len(args)!=2:
    raise ValueError("Please supply exactly 2 input arguments\n%s"%(usage))
unsafe, interesting = args

#-------------------------------------------------

### read in unsafe channel list
file_obj = open(unsafe, "r")
unsafe_chans = defaultdict( set() )
for chan in file_obj:
    chan = chan.strip()
    if chan[2] == "-": ### interpret at KW channel name -> convert!
        chan = chan.split("_")
        ifo, chan = chan[0], "%s"%("_".join(chan[1:-2]))
    else:
        ifo, chan = chan.split(":")
    unsafe_chans[ifo].add( chan )
file_obj.close()

#-------------------------------------------------

### read in interesting channel list and parse
file_obj = open(interesting, "r")
for channel in file_obj:
    channel = channel.strip()
    chan = channel
    if chan[2] == "-": ### interpret at KW channel name -> convert!
        chan = chan.split("_")
        ifo, chan = chan[0], "%s"%("_".join(chan[1:-2]))
    else:
        ifo, chan = chan.split(":")

    if chan in unsafe_chans[ifo]:
        print >> sys.stderr, channel
    else:
        print >> sys.stdout, channel
file_obj.close()
