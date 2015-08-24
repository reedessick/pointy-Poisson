#!/usr/bin/python

usage = """grand_tour.py [--options] gps gps gps..."""
description="""a script that generates a pointed follow-up of possible auxiliary couplings, sweeping through a variety of windows and snr thresholds"""

import numpy as np

from laldetchar.idq import idq
from laldetchar.idq import event

import greedyCI as gci

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from optparse import OptionParser

from ConfigParser import SafeConfigParser

#=================================================

def grand_tour(gps, trgs, trgtype, windows, conf=0.68, plot=False):
    """
    a function that compute the grand tour statistics 
    """
    if trgtype == "kw":
        col = event.col_kw
        snrkey = 'signif'
    elif trgtype == "Omicron":
        col = event.col_snglBurst
        snrkey = 'snr'
    else:
        raise ValueError("do not understand trgtype=%s"%(trgtype))

    figax = []
    for win in windows:
        print "window=%.6f"%win

        ### downelect to only this window
        trgs = [ trg for trg in trgs if (trg[col['tcent']] >= gps-win) and (trg[col['tcent']] <= gps+win) ]

        Pvalue = np.infty
        Snr = None
        Dt = None
        Frq = None
        for snrThr in sorted(list(set([ trg[col[snrkey]] for trg in trgs ]))): ### iterate over all snrs present
            print "snrThr=%.6f"%snrThr

            ctrgs =  [ trg for trg in trgs if trg[snrkey] >= snrThr ]
            n = len( ctrgs ) ### number of triggers
            r = 0.5*n/opts.window ### point estimate of the rate

            if n:
                dt = np.array([ctrgs[col['tcent']] for trg in ctrgs]) - gps
                arg = np.argmin(np.abs(dt))

                snr = ctrgs[arg][col[snrkey]]
                frq = ctrgs[arg][col['fcent']]
                min_dt = dt[arg]
                absmin_dt = abs(min_dt)
            else:
                min_dt = win
                snr = None
                frq = None

            if r > 0:
                pvalue = 1 - np.exp(-r*2*absmin_dt) ### cumulative probability of observing min_dt <= observed min_dt | estimated rate, the factor of 2 comes from looking on either side of the specified gps time
            else:
                pvalue = 1 ### limit is not great here...need to add CI

            print "\n\tchannel=%s\n\t-> Ntrg=%d\n\t-> rate=%.9e Hz\n\t-> min_dt=%.9e sec\n\t-> pvalue=%.9e"%(chan, n, r, min_dt, pvalue)

            r_l, r_h = np.array( gci.poisson_bs(conf, n) ) * 0.5 / win
            pvalue_l = 1 - np.exp(-r_l*2*absmin_dt)
            pvalue_h = 1 - np.exp(-r_h*2*absmin_dt)

            print "\t-> %.5e confidence:\n\t\tlow rate =%.9e Hz\n\t\thigh rate=%.9e Hz\n\t\tlow pvalue =%.9e\n\t\thigh pvalue=%.9e"%(conf, r_l, r_h, pvalue_l, pvalue_h)

            if pvalue_h < Pvalue:
                pvalue = pvalue_h
                Snr = snr
                Dt = dt
                Frq = frq

        ### separate plot for each window
        ### xaxis : time (relative to gps)
        ### yaxis SNR
        ### color = pvalue (upper limit?)
        if plot:                 
            fig = plt.figure()
            ax = plt.subplot(1,1,1)
            figax.append( (fig, ax) )

            ### plot all triggers
#            dts = np.array([ trg[col['tcent'] for trg in trgs ]) - gps
#            snrs = [ trg[col[snrkey]] for trg in trgs ]
#            frqs = [ trg[col['fcent']] for trg in trgs ]

            for trg in trgs:
                dt  = trg[col['tcent']] - gps
                snr = trg[col[snrkey]]
                frq = trg[col['fcent']]

                color = snr_map( snr )
                ax.plot( dt, frq, markerfacecolor=color, markeredgecolor='none', marker='o', linestyle='none', alpha=0.50, markersize=2 )

            ax.set_xlabel( 'time relative to %.6f [sec]'%gps )
            ax.set_ylabel( 'frequency [Hz]' )

            if Dt:
                color = snr_map( Snr )
                ax.plot( Dt, Frq, markerfacecolor='none', markeredgecolor=color, marker='o', linestyle='none', alpha=1.00, markersize=5 )
                ax.text( Dt, Frq, '%.3e'%Pvalue, ha='left', va='center' )

            ax.set_xlim(xmin=-win, xmax=win)
            ymax = 1
            maxfrq = max( [ trg[col['fcent']] for trg in trgs ] )
            while maxfrq > ymax:
                ymax *= 2
            ax.set_ylim(ymin=0, ymax=ymax)

    return figax
 

def snr_map( snr ):
    """
    map snr into color
    """
    return 'k' 
 
#=================================================

parser = OptionParser(usage="", description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("", "--kwverbose", default=False, action="store_true", help="make the retrieve_kwtrigs() call verbose")
parser.add_option("", "--Omicronverbose", dest="overbose", default=False, action="store_true", help="make the retrieve_OmicronTrigs() call verbose")
parser.add_option("", "--OfflineOmicronverbose", dest="ooverbose", default=False, action="store_true", help="make the retrieve_OfflineOmicronTrigs() call verbose")

parser.add_option("-c", "--config", default="config.ini", type="string")

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
windows = sorted([float(l) for l in config.get('general','windows').split()], reverse=True)

maxwindow = max(windows)

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
    if opts.verbose:
        print "gps : %.9f"%(gps)

    #=============================================
    # KW triggers
    #=============================================
    if kwchannels:
        ### go find triggers
        if opts.verbose:
            print "\tdiscoverying KW triggers within [%.9f, %.9f]"%(gps-maxwindow, gps+maxwindow)
        kwtrgdict = idq.retrieve_kwtrigs(kwgdsdir, kwbasename, int(np.floor(gps-maxwindow)), 2*maxwindow+1, kwstride, verbose=opts.kwverbose)
    
        ### keep only the relevant channels
        if opts.verbose:
            print "\tdownselecting to only the following channels:"
            for chan in kwchannels:
                print "\t\t%s"%chan
        kwtrgdict.keep_channels(kwchannels)

        ### trim the edges
        kwtrgdict.include([[gps-maxwindow, gps+maxwindow]], tcent=event.col_kw['tcent'])
   
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
            print "\tdiscovering Omicron triggers within [%.9f, %.9f]"%(gps-maxwindow, gps+maxwindow)
        otrgdict = idq.retrieve_OmicronTrigs(ogdsdir, ifo, int(np.floor(gps-maxwindow)), 2*maxwindow+1, ostride, ochannels, verbose=opts.overbose)

        ### trim edges
        otrgdict.include([[gps-maxwindow, gps+maxwindow]], tcent=event.col_snglBurst['tcent'])

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
            print "\tdiscovering OfflineOmicron triggers within [%.9f, %.9f]"%(gps-maxwindow, gps+maxwindow)
        ootrgdict = idq.retrieve_OfflineOmicronTrigs(oogdsdir, ifo, int(np.floor(gps-opts.window)), 2*opts.window+1, channels=oochannels, verbose=opts.ooverbose)

        ### trim edges
        ootrgdict.include([[gps-maxwindow, gps+maxwindow]], tcent=event.col_snglBurst['tcent'])

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
    
    ### add kw triggers
    trgdict.add( kwtrgdict )

    ### add Omicron triggers
    trgdict.add( otrgdict )

    ### add OfflineOmicron triggers
    trgdict.add( ootrgdict )

    #=============================================
    # cluster triggers?
    #=============================================

    print "\n\n\tWARNING: you may want to cluster triggers!\n\n"

    #=============================================
    # generate statistics, plots
    #=============================================
    if opts.verbose:
        print "\tcomputing statistics, generating plots"

    for chan in kwchannels:
        grand_tour(gps, trgdict[chan][:], "kw", windows, conf=conf)

    for chan in list(set(ochannels+oochannels)):
        grand_tour(gps, trgdict[chan][:], "Omicron", windows, conf=conf)


