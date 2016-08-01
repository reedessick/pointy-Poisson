#!/usr/bin/python
usage = "clusterMultiPopVectors.py [--options] vectors.txt"
description = "reads in vectors from an ascii file (produced by buildMultiPopVectors.py) and attempts to cluster them using k-means."
author = "reed.essick@ligo.org"

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import numpy as np

from scipy import cluster

from collections import defaultdict

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-k", "--kmeans-maxk", default=10, type="int", help="the maximum number of clusters to use in the k-means algorithm")
parser.add_option("", "--kmeans-iter", default=20, type="int", help="the number of iterations to use in the k-means algorithm")
parser.add_option("", "--kmeans-thr", default=1e-5, type="float", help="the threshold used in the k-means algorithm. This is the stopping condition based on the change of distortion in the k-means assignments")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError("Please supply either one or two input arguments\n%s"%(usage))
vectors = args[0]

if opts.tag:
    opts.tag = "_"+opts.tag

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
    v = [float(_) for _ in line[1:]]
    if np.all( np.array(v) < np.infty ):
        vects.append( v ) 
        names.append( line[0] )
    else:
        print "WARNING: found infinite pvalue for : %s"%(line[0])
vects = np.array(vects)
if opts.verbose:
    print "    found %d samples"%(len(vects))

if len(vects) < opts.kmeans_maxk:
    raise ValueError("kmeans_maxk > len(vects) => you *will* have empty clusters")

#-------------------------------------------------

### perform kmeans decomposition
codebooks = []
associations = []
for k in xrange(1, opts.kmeans_maxk+1):
    if opts.verbose:
        print "\nperforming k-means algorithm with :\n    k=%d\n    iter=%d\n    thr=%e"%(k, opts.kmeans_iter, opts.kmeans_thr)
    codebook, distortion = cluster.vq.kmeans( vects, k, iter=opts.kmeans_iter, thresh=opts.kmeans_thr )
    if opts.verbose:
        print "distortion = %.6e"%(distortion)

    codebooks.append( (codebook, distortion) )

    ### print associations between filenames and codewords
    associations.append( defaultdict( list ) )
    for filename, codeword in zip(names, cluster.vq.vq(vects, codebook)[0]):
        associations[-1][codeword].append( filename )
    for key in associations[-1].keys():
        associations[-1][key].sort()

    if opts.verbose:
        print "Associations"
        for codeword in xrange(k):
            print "codeword : %d"%(codeword)
            for filename in associations[-1][codeword]:
                print "    %s"%(filename)
   
#-------------------------------------------------
   
if opts.verbose:
    print "plotting graph!"

colors = ['b', 'r', 'g', 'c', 'm', 'y']

fig = plt.figure()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

old_order = dict( (chan, i) for i, chan in enumerate(sorted(names)) )
for k in xrange(1, opts.kmeans_maxk):
    if opts.verbose:
        print "    working on decomposition with %d clusters"%(k)

    ### figure out the new ordering
    ass = associations[k-1].items()
    ass.sort( key=lambda l: len(l[1]), reverse=True )
    i = 0
    order = {}
    for codeword, chans in ass:
        for chan in chans:
            order[chan] = i
            ax.plot( [old_order[chan], i], [k-1, k], 'k', alpha=0.5 ) ### plot the line connecting the old position to the new!
                                                           ### this probably dominates the wall-clock time...
            i += 1
    old_order = order

    ### plot the boundaries of different codewords
    i = 0
    for ind, (codeword, chans) in enumerate(ass):
        I = i+len(chans)
        ax.plot([i, I], [k, k], color=colors[ind%len(colors)], linewidth=2)
        i = I

ax.set_ylabel('number of clusters')
ax.set_xlabel('sorted associations')

ax.set_ylim(ymin=-0.5, ymax=opts.kmeans_maxk-0.5)
ax.set_xlim(xmin=-0.5, xmax=len(names)+0.5 )

ax.set_yticks(range(opts.kmeans_maxk))

figwidth = figheight = max(len(names)*0.1 + 1./len(names), opts.kmeans_maxk*0.2 + 1./opts.kmeans_maxk)
plt.setp(fig, figwidth=figwidth, figheight=figheight)

figname = "%s/kmeansGraph%s.png"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "    "+figname
fig.savefig( figname )
plt.close( fig ) 
