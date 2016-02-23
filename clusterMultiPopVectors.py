#!/usr/bin/python
usage = "clusterMultiPopVectors.py [--options] vectors.txt"
description = "reads in vectors from an ascii file (produced by buildMultiPopVectors.py) and attempts to cluster them using k-means."
author = "reed.essick@ligo.org"

import numpy as np

from scipy import cluster

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-k", "--kmeans-k", default=10, type="int", help="the number of clusters to use in the k-means algorithm")
parser.add_option("", "--kmeans-iter", default=20, type="int", help="the number of iterations to use in the k-means algorithm")
parser.add_option("", "--kmeans-thr", default=1e-5, type="float", help="the threshold used in the k-means algorithm. This is the stopping condition based on the change of distortion in the k-means assignments")

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError("Please supply exactly one input argument\n%s"%(usage))
vectors = args[0]

#-------------------------------------------------

if opts.verbose:
    print "reading vectors from : %s"%(vectors)

### extract headers
file_obj = open(vectors, "r")
chans = file_obj.readline().strip().split()
file_obj.close()
Nchans = len(chans)
if opts.verbose:
    print "    found %d channels"%(Nchans)

### load vectors
vects = np.loadtxt(vectors, skiprows=1, dtype=float)
shape = np.shape(vects)
if shape[0] != Nchans:
    raise ValueError("inconsistent data format in : %s"%(vectors))
if opts.verbose:
    print "    found %d samples"%(shape[1])

#-------------------------------------------------

### perform kmeans decomposition
if opts.verbose:
    print "performing k-means algorithm with :\n    k=%d\n    iter=%d\n    thr=%e"%(opts.kmeans_k, opts.kmeans_iter, opts.kmeans_thr)
codebook, distortion = cluster.vq.kmeans( vects, opts.kmeans_k, iter=opts.kmeans_iter, thresh=opts.kmeans_thr )

#-------------------------------------------------

### present results!
'''
this should return you a set of centroids for the cluster, which you'll then need to report somehow...
    report the non-zero (or non-trivial) dimensions of each vector?
    have this automatically define OmegaScan config files for each different cluster?
        -> need to associate rows with clusters and then define OmegaScan jobs based on gps times (require these as an input list?)
'''
