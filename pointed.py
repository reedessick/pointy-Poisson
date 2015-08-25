#!/usr/bin/python

usage = """pointed.py [--options] gps gps gps..."""
description="""a script that generates a pointed follow-up of possible auxiliary couplings"""

import numpy as np

from laldetchar.idq import idq
from laldetchar.idq import event

import greedyCI as gci

from optparse import OptionParser

from ConfigParser import SafeConfigParser

#=================================================

parser = OptionParser(usage="", description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("", "--kwverbose", default=False, action="store_true", help="make the retrieve_kwtrigs() call verbose")
parser.add_option("", "--Omicronverbose", dest="overbose", default=False, action="store_true", help="make the retrieve_OmicronTrigs() call verbose")
parser.add_option("", "--OfflineOmicronverbose", dest="ooverbose", default=False, action="store_true", help="make the retrieve_OfflineOmicronTrigs() call verbose")

parser.add_option("-c", "--config", default="config.ini", type="string")
parser.add_option("-w", "--window", default=10, type="float")

parser.add_option("", "--pvalue-print-thr", default=1.0, type="float", help="only print channels/pvalues if they are smaller than this")
parser.add_option("-C", "--confidence-intervals", default=False, action="store_true", help="compute confidence intervals for rates and pvalues")

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

conf = config.getfloat('general', 'conf')
ifo = config.get('general', 'ifo')

#===========

kwgdsdir = config.get("kleinewelle", "gdsdir")
kwbasename = config.get("kleinewelle", "basename")
kwstride = config.getint("kleinewelle", "stride")
kwchannels = config.get("kleinewelle", "channels").split()
for chan in kwchannels:
    if not config.has_section(chan):
        raise ValueError("no section for channel=%s found in %s"%(chan, opts.config))

#===========

oogdsdir = config.get('OfflineOmicron', 'gdsdir')
oochannels = config.get('OfflineOmicron', 'channels').split()
for chan in oochannels:
    if not config.has_section(chan):
        raise ValueError('no section for channel=%s in %s'%(chan, opts.config))

#===========

ogdsdir = config.get('Omicron', 'gdsdir')
ostride = config.getint('Omicron', 'stride')
ochannels = config.get('Omicron', 'channels').split()
for chan in ochannels:
    if not config.has_section(chan):
        raise ValueError('no section for channel=%s in %s'%(chan, opts.config))

if oochannels and ochannels:
    print "WARNING: you've specified both Omicron and OfflineOmicron channels. In the event of a conflict, the OfflineOmicron data will be preferred!"

#=================================================

for gps in args:
    print "gps : %.9f"%(gps)

    #=============================================
    # KW triggers
    #=============================================
    if kwchannels:
        ### go find triggers
        if opts.verbose:
            print "\tdiscoverying KW triggers within [%.9f, %.9f]"%(gps-opts.window, gps+opts.window)
        kwtrgdict = idq.retrieve_kwtrigs(kwgdsdir, kwbasename, int(np.floor(gps-opts.window)), 2*opts.window+1, kwstride, verbose=opts.kwverbose)
    
        ### keep only the relevant channels
        if opts.verbose:
            print "\tdownselecting to only the following channels:"
            for chan in kwchannels:
                print "\t\t%s"%chan
        kwtrgdict.keep_channels(kwchannels)

        ### trim the edges
        kwtrgdict.include([[gps-opts.window, gps+opts.window]], tcent=event.col_kw['tcent'])
   
        ### ensure we have entries for all requested channels and downselect as needed
        if opts.verbose:
            print "\tdownselecting triggers:"
        for chan in kwchannels:
            if kwtrgdict.has_key(chan): ### apply windows, thresholds
                if opts.verbose:
                    print "\t\tchannel=%s, found %d triggers"%(chan, len(kwtrgdict[chan]))

                signifmin = config.getfloat(chan, "signifmin")
                signifmax = config.getfloat(chan, "signifmax")
                kwtrgdict[chan] = [trg for trg in kwtrgdict[chan] if (trg[event.col_kw['signif']] >= signifmin) and (trg[event.col_kw['signif']] <= signifmax) ]

                fmin = config.getfloat(chan, "fmin")
                fmax = config.getfloat(chan, "fmax")
                kwtrgdict[chan] = [trg for trg in kwtrgdict[chan] if (trg[event.col_kw['fcent']] >= fmin) and (trg[event.col_kw['fcent']] <= fmax) ]

                durmin = config.getfloat(chan, "durmin")
                durmax = config.getfloat(chan, "durmax")
                kwtrgdict[chan] = [trg for trg in kwtrgdict[chan] if (trg[event.col_kw['tstop']]-trg[event.col_kw['tstart']] >= durmin) and (trg[event.col_kw['tstop']]-trg[event.col_kw['tstart']] <= durmax) ]

            else:
                if opts.verbose:
                    print "\t\tWARNING: channel=%s not found, inserting an empty list"%chan
                kwtrgdict[chan] = []

            if opts.verbose:
                print "\t\t\tchannel=%s -> %d triggers"%(chan, len(kwtrgdict[chan]))
    else:
        kwtrgdict = event.trigdict()

    #=============================================
    # Omicron triggers
    #=============================================
    if ochannels:
        if opts.verbose:
            print "\tdiscovering Omicron triggers within [%.9f, %.9f]"%(gps-opts.window, gps+opts.window)
        otrgdict = idq.retrieve_OmicronTrigs(ogdsdir, ifo, int(np.floor(gps-opts.window)), 2*opts.window+1, ostride, ochannels, verbose=opts.overbose)

        ### trim edges
        otrgdict.include([[gps-opts.window, gps+opts.window]], tcent=event.col_snglBurst['tcent'])

        ### downselect as needed
        if opts.verbose:
            print "\tdownselecting triggers"
        for chan in ochannels:
            if otrgdict.has_key(chan):
                if opts.verbose:
                    print "\t\tchannel=%s, found %d triggers"%(chan, len(otrgdict[chan]))

                snrmin = config.getfloat(chan, "snrmin")
                snrmax = config.getfloat(chan, "snrmax")
                otrgdict[chan] = [trg for trg in otrgdict[chan] if (trg[event.col_snglBurst['snr']] >= snrmin) and (trg[event.col_snglBurst['snr']] <= snrmax) ]

                fmin = config.getfloat(chan, "fmin")
                fmax = config.getfloat(chan, "fmax")
                otrgdict[chan] = [trg for trg in otrgdict[chan] if (trg[event.col_snglBurst['fcent']] >= fmin) and (trg[event.col_snglBurst['fcent']] <= fmax) ]

                durmin = config.getfloat(chan, "durmin")
                durmax = config.getfloat(chan, "durmax")
                otrgdict[chan] = [trg for trg in otrgdict[chan] if (trg[event.col_snglBurst['duration']] >= durmin) and (trg[event.col_snglBurst['duration']] <= durmax) ]

            else:
                if opts.verbose:
                    print "\t\tWARNING: channel=%s not found, inserting an empty list"%chan
                otrgdict[chan] = []

            if opts.verbose:
                print "\t\t\tchannel=%s -> %d triggers"%(chan, len(otrgdict[chan]))
    else:
        otrgdict = event.trigdict()

    #=============================================
    # OfflineOmicron triggers
    #=============================================
    if oochannels:
        if opts.verbose:
            print "\tdiscovering OfflineOmicron triggers within [%.9f, %.9f]"%(gps-opts.window, gps+opts.window)
        ootrgdict = idq.retrieve_OfflineOmicronTrigs(oogdsdir, ifo, int(np.floor(gps-opts.window)), 2*opts.window+1, channels=oochannels, verbose=opts.ooverbose)

        ### trim edges
        ootrgdict.include([[gps-opts.window, gps+opts.window]], tcent=event.col_snglBurst['tcent'])

        ### downselect as needed
        if opts.verbose:
            print "\tdownselecting triggers"
        for chan in oochannels:
            if ootrgdict.has_key(chan):
                if opts.verbose:
                    print "\t\tchannel=%s, found %d triggers"%(chan, len(ootrgdict[chan]))

                snrmin = config.getfloat(chan, "snrmin")
                snrmax = config.getfloat(chan, "snrmax")
                ootrgdict[chan] = [trg for trg in ootrgdict[chan] if (trg[event.col_snglBurst['snr']] >= snrmin) and (trg[event.col_snglBurst['snr']] <= snrmax) ]

                fmin = config.getfloat(chan, "fmin")
                fmax = config.getfloat(chan, "fmax")
                ootrgdict[chan] = [trg for trg in ootrgdict[chan] if (trg[event.col_snglBurst['fcent']] >= fmin) and (trg[event.col_snglBurst['fcent']] <= fmax) ]

                durmin = config.getfloat(chan, "durmin")
                durmax = config.getfloat(chan, "durmax")
                ootrgdict[chan] = [trg for trg in ootrgdict[chan] if (trg[event.col_snglBurst['duration']] >= durmin) and (trg[event.col_snglBurst['duration']] <= durmax) ]

            else:
                if opts.verbose:
                    print "\t\tWARNING: channel=%s not found, inserting an empty list"%chan
                ootrgdict[chan] = []

            if opts.verbose:
                print "\t\t\tchannel=%s -> %d triggers"%(chan, len(ootrgdict[chan]))

    else:
        ootrgdict = event.trigdict()

    #=============================================
    # combine all trgdicts
    #=============================================
    trgdict = event.trigdict()
    channels = []
    
    ### add kw triggers
    trgdict.add( kwtrgdict )
    channels += kwchannels

    ### add Omicron triggers
    trgdict.add( otrgdict )
    channels += ochannels

    ### add OfflineOmicron triggers
    trgdict.add( ootrgdict )
    channels += oochannels

    ### clean up channel list to get a unique set
    channels = list(set(channels))

    #=============================================
    # cluster triggers?
    #=============================================

#    print "\n\n\tWARNING: you may want to cluster triggers!\n\n"

    #=============================================
    # generate statistics, plots
    #=============================================
    if opts.verbose:
        print "\tcomputing statistics, generating plots"

    for chan in channels:
        n = len(trgdict[chan]) ### number of triggers
        r = 0.5*n/opts.window ### point estimate of the rate

        if n:
            dt = np.array([trg[event.col_kw['tcent']] for trg in trgdict[chan]]) - gps
            arg = np.argmin(np.abs(dt))
            min_dt = dt[arg]
        else:
            min_dt = opts.window
        absmin_dt = abs(min_dt)

        if r > 0:
            pvalue = 1 - np.exp(-r*2*absmin_dt) ### cumulative probability of observing min_dt <= observed min_dt | estimated rate, the factor of 2 comes from looking on either side of the specified gps time
        else:
            pvalue = 1 ### limit is not great here...need to add CI

        if (pvalue <= opts.pvalue_print_thr):
            print "\n\tchannel=%s\n\t-> Ntrg=%d\n\t-> rate=%.9e Hz\n\t-> min_dt=%.9e sec\n\t-> pvalue=%.9e"%(chan, n, r, min_dt, pvalue)

        if opts.confidence_intervals:
            r_l, r_h = np.array( gci.poisson_bs(conf, n) ) * 0.5 / opts.window
            pvalue_l = 1 - np.exp(-r_l*2*absmin_dt)
            pvalue_h = 1 - np.exp(-r_h*2*absmin_dt)

            print "\t-> %.5e confidence:\n\t\tlow rate =%.9e Hz\n\t\thigh rate=%.9e Hz\n\t\tlow pvalue =%.9e\n\t\thigh pvalue=%.9e"%(conf, r_l, r_h, pvalue_l, pvalue_h)







    ### SHOULD BE ABLE TO TEST THAT THIS STATISTIC DOES WHAT WE THINK IT SHOULD BY CHECKING RANDOM TIMES AND THE DISTRIBUTION DERIVED THEREFROM
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





