#!/usr/bin/python

usage = """pointed.py [--options] gps gps gps..."""
description="""a script that generates a pointed follow-up of possible auxiliary couplings"""

import numpy as np

from laldetchar.idq import idq
from laldetchar.idq import event

from optparse import OptionParser

from ConfigParser import SafeConfigParser

#=================================================

parser = OptionParser(usage="", description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--config", default="config.ini", type="string")
parser.add_option("-w", "--window", default=10, type="float")

parser.add_option("-o", "--output-dir", default="./", type="string")

opts, args = parser.parse_args()

if not len(args):
    if opts.verbose:
        print "no gps times specified"
    import sys
    sys.exit(0)
else:
    args = [float(arg) for arg in args]

#=================================================

config = SafeConfigParser()
config.read(opts.config)

gdsdir = config.get("general", "gdsdir")
kwbasename = config.get("general", "kwbasename")
kwstride = config.getint("general", "kwstride")

channels = config.get("general", "channels").split()

for chan in channels:
    if not config.has_section(chan):
        raise ValueError("no section for channel=%s found in %s"%(chan, opts.config))

#=================================================

for gps in args:

    ### go find triggers
    if opts.verbose:
        print "gps : %.9f"%gps
        print "\tdiscoverying triggers within [%.9f, %.9f]"%(gps -opts.window, gps+opts.window)
    trgdict = idq.retrieve_kwtrigs(gdsdir, kwbasename, int(np.floor(gps-opts.window)), 2*opts.window+1, kwstride)
    
    ### keep only the relevant channels
    if opts.verbose:
        print "\tdownselecting to only the following channels:"
        for chan in channels:
            print "\t\t%s"%chan
    trgdict.keep_channels(channels)

    ### trim the edges
    trgdict.include([[gps-opts.window, gps+opts.window]])
   

    print "\n\n\t\tWARNING: you may want to cluster triggers!\n\n"

 
    ### ensure we have all the channels requested
    if opts.verbose:
        print "\tdownselecting triggers:"
    for chan in channels:
        if trgdict.has_key(chan): ### apply windows, thresholds
            if opts.verbose:
                print "\t\tchannel=%s, found %d triggers"%(chan, len(trgdict[chan]))

            signifmin = config.getfloat(chan, "signifmin")
            signifmax = config.getfloat(chan, "signifmax")
            trgdict[chan] = [trg for trg in trgdict[chan] if (trg[event.col_kw['signif']] >= signifmin) and (trg[event.col_kw['signif']] <= signifmax) ]

            fmin = config.getfloat(chan, "fmin")
            fmax = config.getfloat(chan, "fmax")
            trgdict[chan] = [trg for trg in trgdict[chan] if (trg[event.col_kw['fcent']] >= fmin) and (trg[event.col_kw['fcent']] <= fmax) ]

            durmin = config.getfloat(chan, "durmin")
            durmax = config.getfloat(chan, "durmax")
            trgdict[chan] = [trg for trg in trgdict[chan] if (trg[event.col_kw['tstop']]-trg[event.col_kw['tstart']] >= durmin) and (trg[event.col_kw['tstop']]-trg[event.col_kw['tstart']] <= durmax) ]

        else:
            if opts.verbose:
                print "\t\tWARNING: channel=%s not found, inserting an empty list"%chan
            trgdict[chan] = []


        if opts.verbose:
            print "\t\t\tchannel=%s -> %d triggers"%(chan, len(trgdict[chan]))

    if opts.verbose:
        print "\tcomputing statistics, generating plots"
    for chan in channels:
        n = len(trgdict[chan]) ### number of triggers
        r = 0.5*n/opts.window ### point estimate of the rate

        print "\n\n\t\tWARNING: need to compute CI for poisson rates and propogate this through to pvalues!\n\n"
        ### closed form posterior based on a complimentary prior (Gamma Distribution)?


        if n:
            min_dt = min(np.abs(np.array([trg[event.col_kw['tcent']] for trg in trgdict[chan]]) - gps))
        else:
            min_dt = opts.window

        if r > 0:
            pvalue = 1 - np.exp(-r*min_dt) ### cumulative probability of observing min_dt <= observed min_dt | estimated rate
        else:
            pvalue = 1 ### limit is not great here...need to add CI
        ### SHOULD BE ABLE TO TEST THAT THIS STATISTIC DOES WHAT WE THINK IT SHOULD BY CHECKING RANDOM TIMES AND THE DISTRIBUTION DERIVED THEREFROM

        print "channel=%s\n\t-> Ntrg=%d\n\t-> rate=%.9e Hz\n\t-> min_dt=%.9e sec\n\t-> pvalue=%.9e"%(chan, n, r, min_dt, pvalue)

    ### DO SOMETHING
    '''
plot histograms
    TIME-DOMAIN ONLY
	zero scale so t=0 corresponds to the gps time
    FREQ-DOMAIN ONLY
	
    TIME-FREQ DOMAIN ONLY

    ALL
	PLOT FOR EACH CHANNEL SEPARATELY, in some sort of stacked way
	PLOT for each gps separately?

	markers for each trigger separately
	stacke multiple gps times into a global histogram
	compute and report basic statistics:
		aux rate (need to estimate this over a large window... rely on the usre to do this correctly?)
			set CI on this rate
		poisson probability of finding a trigger within "x" of the requested time. set "x" to be the closest time?
			set CI on this statistic using CI on the rate

scatter aux SNR vs DARM SNR
	--> possibly inform safety criteria that we can eventually add back into the downselection of events when defining a veto

histogram the dt between triggers in each channel -> confirm "poisson-ianity" of the triggers and the validity of our pvalue estimate

    '''





