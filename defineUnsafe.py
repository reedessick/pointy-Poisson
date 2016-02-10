usage = "python defineUnsafe.py [--options] pointy.out pointy.out pointy.out ..."
description = "stacks pvalues to determine safety"
author = "Reed Essick (reed.essick@ligo.org)"

import glob
import numpy as np
from collections import defaultdict

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-p", "--pvalue", default=1, type="float", help="the pvalue threshold used to define participation")

parser.add_option("-P", "--plot", default=False, action="store_true")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-n", "--nperbin", default=8, type="int")

parser.add_option("", "--xmin", default=None, type="float")

parser.add_option("", "--snglchan-histograms", default=False, action="store_true")

parser.add_option("", "--present-histogram", default=False, action="store_true")

parser.add_option("", "--expected-unsafe", default=None, type="string", help="only used to determine colors on histograms")

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_%s"%(opts.tag)

nargs = len(args)
if not nargs:
    raise ValueError("please supply at least one pointy.out file as an input argument")

if opts.plot or opts.snglchan_histograms:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

if opts.expected_unsafe != None:
    file_obj = open(opts.expected_unsafe, "r")
    expected_unsafe = [line.strip() for line in file_obj.readlines() if line.strip()]
    file_obj.close()
else:
    expected_unsafe = []

#=================================================
### load in the lists
if opts.verbose:
    print "loading observed pvalues from :"

chans = defaultdict( list )
for filename in args:
    if opts.verbose:
        print "\t%s"%(filename)
    file_obj = open(filename, "r")
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
            p = float(lines[nind].strip().split("=")[-1])
            chans[chan].append( (p, lines[ind:nind+1], filename) )
            ind = ind
        ind += 1
Ntrials = len(args) ### number of trials

if opts.snglchan_histograms:
    if opts.verbose:
        print "plotting single-channel histograms:"
    for chan, pvalues in chans.items():
        pvalues = np.log10([p[0] for p in pvalues])

        fig = plt.figure()
        ax = fig.gca()

        n = len(pvalues)
        ax.hist( pvalues )

        ax.plot( [np.exp(np.sum(np.log(pvalues)/Ntrials))]*2, ax.get_ylim(), 'k-', linewidth=2, alpha=0.5 )

        ax.set_xlabel('log10(pvalue)')
        ax.set_ylabel('count')

        ax.set_title(chan.replace("_","\_"))

        fig.text( 0.15, 0.95, 'N = %d'%(n), ha='left', va='top' )

        figname = "%s/%04d%s%s.png"%(opts.output_dir, n, chan, opts.tag)
        if opts.verbose:
            print "    %s\tN=%d"%(figname, n)
        fig.savefig( figname )
        plt.close( fig )

if opts.present_histogram:
    if opts.verbose:
        print "plotting histogram of number of appearences"

    fig = plt.figure()
    ax = fig.gca()

    ax.hist( [len(pvalues) for pvalues in chans.values()], histtype='step' )

    ax.set_xlabel('No. of appearences')
    ax.set_ylabel('count')

    ax.set_yscale('log')
    ax.set_ylim(ymin=0.8)

    figname = "%s/appearences%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "    %s"%(figname)
    fig.savefig( figname )
    plt.close( fig )

chans = dict( [ (key, np.exp(np.sum([np.log(p) for p, _, _ in chans[key]])/Ntrials) ) for key in chans.keys() ] )

#=================================================
### write the summary file
filename = "%s/present%s.out"%(opts.output_dir, opts.tag)
fileNAME = "%s/unsafe%s.out"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "writing :\n\t%s\n\t%s"%(filename, fileNAME)
file_obj = open(filename, "w")
file_OBJ = open(fileNAME, "w")

### order the lists
if opts.plot:
    pvalues = []
    unsafes = []
items = chans.items()
items.sort(key=lambda l:l[1])
for key,pvalue in items:

    print >> file_obj, "%.6e\t%s"%(pvalue, key)
    if pvalue <= opts.pvalue:
        print >> file_OBJ, key

    if opts.plot:
        if key in expected_unsafe:
            unsafes.append( pvalue )
        else:
            pvalues.append( pvalue )

#=================================================
### make the plot
if opts.plot:
    if opts.verbose:
        print "plotting"
    plt.rcParams['text.usetex'] = True

    pvalues = [ pvalue for pvalue in pvalues if pvalue > 0]

    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    nbins = len(pvalues)/opts.nperbin ### could produce poor binning...
    if opts.xmin!=None:
        this_min = opts.xmin
    else:
        this_min = min( pvalues )

    bins = np.logspace( np.log10(this_min), 0, nbins)
#    ax.hist( pvalues, bins=bins, histtype="bar", color='g', log=True )
    if unsafes:
        ax.hist( [pvalues, unsafes], bins=bins, histtype="step", color=['g', 'r'], label=['safe', 'unsafe'], log=True )
    else:
        ax.hist( [pvalues], bins=bins, histtype="step", color=['g'], label=['safe'], log=True )

    ax.set_xlabel('$\Pi_i\mathrm{pvalue}_i^{1/%d}$'%(Ntrials))
    ax.set_ylabel('$\mathrm{count}$')

    ax.set_xscale("log")

    ax.grid(True, which="both")

    ax.plot( [opts.pvalue]*2, plt.gca().get_ylim(), 'k--', linewidth=2 )

    ax.legend(loc='best')

    if opts.xmin!=None:
        ax.set_xlim(xmin=opts.xmin)
    ax.set_xlim(xmax=1.0)
    ax.set_ylim(ymin=0.5)

    figname = "%s/stacked-pvalue%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "\t%s"%(figname)
    fig.savefig(figname)
    plt.close(fig)

