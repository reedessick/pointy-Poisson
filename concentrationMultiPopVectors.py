#!/usr/bin/python
usage = "concentrationMultiPopVectors.py [--options] vectors.txt"
description = "a tool to make concentration diagrams of data stored in vectors.txt"
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

parser.add_option('-g', '--giniThr', default=0.9, type='float', help='a threshold on the gini index. Vectors with larger gini indexes are reported')
parser.add_option('', '--nbins', default=20, type='int', help='the number of bins for the gini index histogram')

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

### plot concentration diagrams and histogram of gini indexes

if opts.verbose:
    print "plotting concentration diagrams"

fig = plt.figure(figsize=(10,6))
axC = plt.subplot(1,2,1)
axH = plt.subplot(1,2,2)

x = (np.arange(Nchans)+1.0)/Nchans
dx = x[1]-x[0]
gini = []
keepers = []
for name, vect in zip(names, vects):
    if np.any(vect):
        y = vect[:]
        y.sort()
        y = y[::-1]
        y = np.cumsum( y**2 )
        y /= y[-1] 
    else:
        y = x
    
    axC.plot( x, y, color='b', alpha=0.5 )

    g = ( (np.sum( (y[1:]+y[:-1])*0.5 ) + 0.5*y[0]) * dx - 0.5 )*2
    if g!=g: raise ValueError

    gini.append( g ) ### approximate the integral with trapzoids and subtract off half for gini index

    if g >= opts.giniThr:
        keepers.append( name )

gini = np.array(gini)
bins = np.logspace( np.log10(1-max(gini)), 0, opts.nbins+1)
axH.hist( 1.0-gini, bins=bins )

### decorate

axC.set_xscale('log')

axC.set_xlabel('fraction of vector (len=%d)'%Nchans)
axC.set_ylabel('fraction of norm')

axH.set_xscale("log")

axH.set_xlabel('1 - (gini index)')
axH.set_ylabel('count')

axC.plot([0,1], [0,1], 'k--')

### save

figname = "%s/concentration%s.png"%(opts.output_dir, opts.tag)
if opts.verbose:
    print figname
fig.savefig(figname)
plt.close(fig)

### write keepers to file
filename = "%s/concentration%s.out"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "writing large gini names to : %s"%filename
file_obj = open(filename, "w")
for name in sorted(keepers):
    print >> file_obj, name
file_obj.close()
