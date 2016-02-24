#!/usr/bin/python
usage = "histMultiPopVectors.py [--options] vectors.txt"
description = "a tool to make graphical representations of data stored in vectors.txt"
author = "reed.essick@ligo.org"

import numpy as np

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--oneD", default=False, action="store_true", help="plot 1D histograms")
parser.add_option("", "--twoD", default=False, action="store_true", help="plot 2D histograms")
parser.add_option("", "--normD", default=False, action="store_true", help="plot histogram of the L2 norms of the vectors")

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

    ax.hist( norms, bins=max(10, len(norms)/10), histtype="step" )

    ax.set_ylabel( 'count' )
    ax.set_xlabel( '$\sum \log p_i$ \sim \log p_\mathrm{joint}$' )

    figname = "%s/normHist%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "saving : %s"%(figname)
    fig.savefig( figname )
    plt.close( fig )

