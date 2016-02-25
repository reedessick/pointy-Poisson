#!/usr/bin/python
usage = "clusterMultiPopVectors.py [--options] vectors.txt"
description = "reads in vectors from an ascii file (produced by buildMultiPopVectors.py) and attempts to cluster them using k-means."
author = "reed.essick@ligo.org"

import numpy as np

from scipy import cluster

from collections import defaultdict

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-k", "--kmeans-k", default=10, type="int", help="the number of clusters to use in the k-means algorithm")
parser.add_option("", "--kmeans-iter", default=20, type="int", help="the number of iterations to use in the k-means algorithm")
parser.add_option("", "--kmeans-thr", default=1e-5, type="float", help="the threshold used in the k-means algorithm. This is the stopping condition based on the change of distortion in the k-means assignments")

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError("Please supply either one or two input arguments\n%s"%(usage))
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

if len(vects) < opts.kmeans_k:
    raise ValueError("kmeans_k > len(vects) => you *will* have empty clusters")

#-------------------------------------------------

### perform kmeans decomposition
if opts.verbose:
    print "performing k-means algorithm with :\n    k=%d\n    iter=%d\n    thr=%e"%(opts.kmeans_k, opts.kmeans_iter, opts.kmeans_thr)
codebook, distortion = cluster.vq.kmeans( vects, opts.kmeans_k, iter=opts.kmeans_iter, thresh=opts.kmeans_thr )

#-------------------------------------------------

### present results!
print "distortion = %.6e\n"%(distortion)

for ind, codeword in enumerate(codebook):
    print "codeword : %d"%(i)
    for i in np.nonzero(codeword)[0]:
        print "    %.6e  %s"%(codeword[i], chans[i])

### print associations between filenames and codewords
associations = defaultdict( list )
for filename, codeword in zip(names, cluster.vq.vq(vects, codebook)):
    associations[codeword].append( filename )
    
print "\nAssociations\n"
for codeword in xrange(opts.kmeans_k):
    print "codeword : %d"%(codeword)
    for filename in associations[codeword]:
        print "    %s"%(filename)

#-------------------------------------------------
    
'''
have this automatically define OmegaScan config files for each different cluster?
    -> need to associate rows with clusters and then define OmegaScan jobs based on gps times (require these as an input list?)
'''
