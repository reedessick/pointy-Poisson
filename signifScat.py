#!/usr/bin/python

usage = """signifScat.py [--options] gps gps gps..."""
description="""a script that generates scatter plots of KW signif amplitudes"""

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams['text.usetex'] = True

import numpy as np
import math

from laldetchar.idq import idq
from laldetchar.idq import event

from collections import defaultdict

from optparse import OptionParser

from ConfigParser import SafeConfigParser

#=================================================

parser = OptionParser(usage="", description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("", "--kwverbose", default=False, action="store_true", help="make the retrieve_kwtrigs() call verbose")

parser.add_option("-c", "--config", default="config.ini", type="string")
parser.add_option("-w", "--window", default=0.500, type="float")

parser.add_option("-o", "--output-dir", default="./", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-f", "--force", default=False, action="store_true", help="continue on past times with no coverage")

parser.add_option("", "--pairs", default=[], action="append", type="string", help="eg: H1:CAL-DELTAL_EXTERNAL_OUT_DQ,H1:ASC-AS_B_RF36_Q_YAW_OUT_DQ")

opts, args = parser.parse_args()

if not len(args):
    if opts.verbose:
        print "no gps times specified"
    import sys
    sys.exit(0)
else:
    args = [float(arg) for arg in args]

if opts.tag:
    opts.tag = "_"+opts.tag

opts,pairs = [pair.split(",") for pair in opts.pairs]

#=================================================

config = SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')

#===========

kwgdsdirs = config.get("kleinewelle", "gdsdirs").split()
kwbasename = config.get("kleinewelle", "basename")
kwstride = config.getint("kleinewelle", "stride")

kwsignifmin = config.getfloat("kleinewelle", "signifmin")
kwsignifmax = config.getfloat("kleinewelle", "signifmax")

fmin = config.getfloat("kleinewelle", "fmin")
fmax = config.getfloat("kleinewelle", "fmax")

durmin = config.getfloat("kleinewelle", "durmin")
durmax = config.getfloat("kleinewelle", "durmax")

#===========

### iterate over gps and load timing differences
trgdata = {}
for gps in args:
    print "gps : %.9f"%(gps)

    minwin = opts.window

    ### go find triggers
    if opts.verbose:
        print "\tdiscoverying KW triggers within [%.9f, %.9f]"%(gps-opts.window, gps+opts.window)

    ### figure out which files you want
    filenames = []
    coverage = []
    for gdsdir in kwgdsdirs:
        for filename in idq.get_all_files_in_range(gdsdir, gps-opts.window, gps+opts.window, pad=0, suffix=".trg"):
            seg = idq.extract_start_stop(filename, suffix=".trg")
            if not event.livetime(event.andsegments([coverage, [seg]])):
                coverage = event.fixsegments( coverage + [seg] )
                filenames.append( filename )

    ### figure out the extent of the coverage
    if len(event.include([[gps]], coverage, tcent=0)) == 0:
        if opts.force:
            if opts.verbose:
                print "no triggers found for gps : %.3f"%(gps)
            continue
        else:
            raise ValueError("no triggers found for gps : %.3f"%(gps))
    for s, e in coverage:
        if s < gps:
            if gps-s < minwin:
                minwin = gps-s
        elif e > gps:
            if e-gps < minwin:
                minwin = e-gps
    else:
        print "\t\tminwin = %.3f"%(minwin)
 
    ### load triggers
    trgdict = event.trigdict()
    for filename in filenames:
        if opts.verbose:
            print "\t\t%s"%(filename)
        trgdict.add( event.loadkwm( filename ) )
    trgdict.include([[gps-opts.window, gps+opts.window]], tcent=event.col_kw['tcent']) ### keep only those within big window
    for key, values in trgdict.items(): ### add shifted times so they are measured relative to gps
        trgdict[key] = [trg+[trg[event.col_kw['tcent']]-gps] for trg in values]
    trgdata[gps] = (minwin, trgdict)

#-----------------------------------

if opts.verbose:
    print "plotting signif scatters for :"

for chanA, chanB in opts.pairs:
    if opts.verbose:
        print "    %s\n    %s"%(chanA, chanB)


