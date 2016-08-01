#!/usr/bin/python
usage = "histMultiPopVectors.py [--options] vectors.txt"
description = "a tool to make graphical representations of data stored in vectors.txt"
author = "reed.essick@ligo.org"

import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
plt.rcParams['text.usetex'] = True

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--oneD", default=False, action="store_true", help="plot 1D histograms of each feature")

parser.add_option("", "--twoD", default=False, action="store_true", help="plot 2D histograms of all pairs of features")

parser.add_option("", "--normD", default=False, action="store_true", help="plot histogram of the L2 norms of the vectors")
parser.add_option("", "--normThr", default=1e2, type="float", help="print the filenames of vectors associated with norms larger than --normThr. Must be used with --normD to have any effect")

parser.add_option("", "--hammingD", default=False, action="store_true", help="plot histogram of the hamming distance of vectors's participation")
parser.add_option("", "--hammingThr", default=3, type="int", help="print hte filenames of vectors associated with hamming distancances larger than --hammingThr. Must be used with --hammingD to have an effect")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_%s"%(opts.tag)

if len(args)!=1:
    raise ValueError("Please supply exactly one input argument\n%s"%(usage))
vectors = args[0]

#-------------------------------------------------

if opts.verbose:
    print "reading vectors from : %s"%(vectors)

### extract headers
file_obj = open(vectors, "r")
chans = file_obj.readline().strip().split()[1:] ### skip first column because that is the filename
Nchans = len(chans)
if opts.verbose:
    print "    found %d channels"%(Nchans)

### load vectors
vects = []
names = []
for line in file_obj:
    line = line.strip().split()
    if len(line) != Nchans+1:
        raise ValueError("inconsistent data format in : %s"%(vectors))
    vects.append( [float(_) for _ in line[1:]] )
    names.append( line[0] )
vects = np.array(vects)
if opts.verbose:
    print "    found %d samples"%(len(vects))

#-------------------------------------------------

### make 1D histograms of each feature
if opts.oneD:
    if opts.verbose:
        print "plotting 1D histograms"
    raise StandardError("--oneD not yet implemented")

### make 2D histograms of each pair of features
if opts.twoD:
    if opts.verbose:
        print "plotting 2D histograms"
    raise StandardError("--twoD not yet implemented")

### make 1D histogram of L2 norm of vectors ~ log(prod(pvalue))
if opts.normD:
    if opts.verbose:
        print "plotting histogram of L2 norms"
    fig = plt.figure()
    ax = fig.gca()

    norms = np.sum(vects**2, axis=1)

    ### prting intersting times with big norms
    truth = norms >= opts.normThr
    print "%d interesting times with norm>=%.2e"%(np.sum(truth), opts.normThr)
    for i in np.nonzero( truth )[0]:
        print "    "+names[i]

    ### manipulate for plotting
    norms = norms[(norms>0)*(norms<np.infty)]
    norms = np.log10(norms)

    bins = np.linspace(np.min(norms), np.max(norms), max(10, len(norms)))

    n, _, _ = ax.hist( norms, bins=bins, histtype="step", log=True, cumulative=-1, weights=np.ones(len(norms),dtype=float)/len(norms) )

    ax.set_ylim(ymin=min(n)*0.95, ymax=1)

    ax.plot( [np.log10(opts.normThr)]*2, ax.get_ylim(), 'k--')

    ax.set_ylabel( 'cumulative fraction of events' )
    ax.set_xlabel( '$\log_{10} \left(-\sum \log p_i \sim -\log p_\mathrm{joint}\\right)$' )

    ax.grid(True, which="both")

    figname = "%s/normHist%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "saving : %s"%(figname)
    fig.savefig( figname )
    plt.close( fig )


### make hamming distance histograms
if opts.hammingD:
    if opts.verbose:
        print "plotting histogram of hamming distances"

    fig = plt.figure()
    ax = fig.gca()

    hamms = np.sum((vects > 0).astype(int), axis=1)
    truth = hamms >= opts.hammingThr
    print "%d interesting times with hamming>=%d"%(np.sum(truth), opts.hammingThr)
    for i in np.nonzero( truth )[0]:
        print "    "+names[i]

    bins = np.linspace(np.min(hamms), np.max(hamms), max(10, len(hamms)))

    n, _, _ = ax.hist( hamms, bins=bins, histtype="step", log=True, cumulative=-1, weights=np.ones(len(hamms),dtype=float)/len(hamms) )

    ax.set_ylim(ymin=min(n)*0.95, ymax=1)

    ax.plot( [opts.hammingThr]*2, ax.get_ylim(), 'k--')

    ax.set_ylabel( 'cumulative fraction of events' )
    ax.set_xlabel( 'hamming distance' )

    ax.grid(True, which="both")

    ax.set_xscale('log')

    figname = "%s/hammingHist%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "saving : %s"%(figname)
    fig.savefig( figname )
    plt.close( fig )
