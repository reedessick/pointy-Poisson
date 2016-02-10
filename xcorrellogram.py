#!/usr/bin/python

usage = """xcorrellogram.py [--options] gps gps gps..."""
description="""a script that generates cross-correlograms to examine possible auxiliary couplings"""

import numpy as np
import math

from laldetchar.idq import idq
from laldetchar.idq import event

from collections import defaultdict

from optparse import OptionParser

from ConfigParser import SafeConfigParser

#=================================================

def recursive_binning_logp( dts, start, end, Nt, prior="data", maxleafsize=1 ):
    """
    recursively divides the space until there is either zero or 1 elements in each bin
    computes pvalues for each level of binnning/decomposition and returns the set
    """
    N = len(dts)
    w = end-start
    if prior=="uniform":
        raise StandardError("I don't trust this derivation")
        logp = -(1+N)*np.log(w) + fast_logfactorial( N )

    elif prior=="jeffreys":
        raise StandardError("I don't trust this derivation")
        logp = -(0.5+N)*np.log(w) + fast_loggamma( 0.5 + N )

    elif prior=="data":
        logp = (1+N)*np.log(Nt) - (2+N)*np.log(Nt+1) + np.log(1+N) - N*np.log(w)

    else:
        raise ValueError("prior=%s not understood"%(prior))

    dlogp = np.array( recursive_binning_logp_helper( dts, start, end, maxleafsize=maxleafsize ) )
    return logp + dlogp

def fast_logfactorial( N, Nmax=100 ):
    """
    compute factorials of big numbers
    """
    if N < Nmax:
        return math.log(math.factorial(N))
    else:
        return 0.5*np.log(2*np.pi*N) + N*(np.log(N)-1)

def fast_loggamma( N, Nmax=100 ):
    """
    compute gamma function of big numbers
    """
    if N < Nmax:
        return math.log(math.gamma(N))
    else:
        return 0.5*np.log(2*np.pi*(N-1)) + (N-1)*(np.log(N-1)-1) ### stirling's approximation for factorial

def recursive_binning_logp_helper( dts, start, end, maxleafsize=1 ):
    """
    recurse through range and return 
    """
    N = len(dts)
    if N < maxleafsize+1: ### termination condition met!
        dlogp = N*np.log(end-start) - fast_logfactorial(N)
        return [dlogp]
    else: ### compute this level and bisect!
        dlogp = N*np.log(end-start) - fast_logfactorial(N)
        mid = 0.5*(start+end)
        left = []
        right = []
        for _ in dts:
            if _ < mid:
                left.append( _ )
            else:
                right.append( _ )
        leftlogp = recursive_binning_logp_helper( left, start, mid, maxleafsize=maxleafsize )
        rightlogp = recursive_binning_logp_helper( right, mid, end, maxleafsize=maxleafsize )
        ans = [dlogp]

        ans.append( max(leftlogp) + max(rightlogp) )

        #### this is bad... blows up
#        for l in leftlogp:
#            for r in rightlogp:
#                ans.append( l + r )

        return ans

#=================================================

parser = OptionParser(usage="", description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("", "--kwverbose", default=False, action="store_true", help="make the retrieve_kwtrigs() call verbose")

parser.add_option("-c", "--config", default="config.ini", type="string")
parser.add_option("-w", "--window", default=10, type="float")

parser.add_option("-e", "--exclude", default=0.0, type="float", help="exclude triggers falling closer than exclude when computing the rate. This is NOT applied when finding the closest trigger and computing a pvalue, only for the rate estimation.")

parser.add_option("-o", "--output-dir", default="./", type="string")
parser.add_option("-t", "--tag", default="", type="string")
parser.add_option("-p", "--plot", default=False, action="store_true")
parser.add_option("-s", "--stat", default=False, action="store_true")

parser.add_option("-f", "--force", default=False, action="store_true", help="continue on past times with no coverage")

opts, args = parser.parse_args()

if not len(args):
    if opts.verbose:
        print "no gps times specified"
    import sys
    sys.exit(0)
else:
    args = [float(arg) for arg in args]

if opts.exclude >= opts.window:
    raise ValueError("--exclude is larger than --window. That doesn't make sense")

if opts.plot:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams['text.usetex'] = True

if opts.tag:
    opts.tag = "_"+opts.tag

#=================================================

config = SafeConfigParser()
config.read(opts.config)

ifo = config.get('general', 'ifo')

prior = config.get('general', 'prior')
if prior not in ['uniform', 'jeffreys', 'data']:
    raise ValueError('prior=%s not understood'%(prior))

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
    if minwin < opts.exclude:
        raise ValueError("minwin < opts.exclude caused by gps : %3f"%(gps))
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

### find all channels and map them to gps times and minwin
absMinWin = opts.window
channels = defaultdict( list ) 
for gps, (minwin, trgdict) in trgdata.items():
    for chan in trgdict.keys():
        channels[chan].append( (gps, minwin) )
    if minwin < absMinWin:
        absMinWin = minwin
for chan in channels.keys():
    channels[chan] = np.array(channels[chan])

#-------------------------------------

if opts.verbose:
    print "computing correlograms and statistics"

## iterate over channels and compute statistics, plots
chans = sorted(channels.keys())
ans = {}
for chan in sorted(channels.keys()):
    if opts.verbose:
        print "    %s"%(chan)

    dts = []
    for gps, minwin in channels[chan]:
        dts += [ trg[-1] for trg in trgdata[gps][1][chan] if abs(trg[-1]) < absMinWin ]
    if opts.stat:
        N = len(dts)
        Nt = len(args)
        if opts.verbose:
            print "        N = %d"%(N)
        ### compute statistics here!
        logps = recursive_binning_logp( dts, -absMinWin, absMinWin, Nt, prior=prior, maxleafsize=max(1, N/20) )
        if opts.verbose:
            print "        logp = %.6e"%(np.max(logps))
        ans[chan] = (logps, N)

    ### plot!
    if opts.plot:
        fig = plt.figure()
        ax = fig.gca()

        if len(dts):
            ax.hist( dts, bins=min(100, max(10, len(dts)/10)) )
        else:
            ax.text( 0, 0.5, 'no triggers found', ha='center', va='center' )

        ax.set_xlabel('$\Delta t\ [\mathrm{sec}]$')
        ax.set_ylabel('$\mathrm{count}$')

        ax.set_title(chan.replace("_","\_"))

        ax.set_xlim(xmin=-absMinWin, xmax=absMinWin)

        figname = "%s/xcorrel_%s%s.png"%(opts.output_dir, chan, opts.tag)
        if opts.verbose:
            print "\t%s"%(figname)
        fig.savefig( figname )
        plt.close( fig )

if opts.stat:
    fig = plt.figure()
    ax = fig.gca()

    ax.plot( [np.max(ans[chan][0]) for chan in chans], [ans[chan][1] for chan in chans], marker='.', linestyle='none' )

    ax.set_xlabel('log(p)')
    ax.set_ylabel('No. glitches (total)')

    figname = "%s/xcorrel_scaling%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "%s"%(figname)
    fig.savefig( figname )
    plt.close( fig )
