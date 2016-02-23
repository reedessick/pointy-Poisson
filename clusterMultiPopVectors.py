#!/usr/bin/python
usage = "clusterMultiPopVectors.py [--options] vectors.txt"
description = "reads in vectors from an ascii file (produced by buildMultiPopVectors.py) and attempts to cluster them using k-means."
author = "reed.essick@ligo.org"

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

opts, args = parser.parse_args()

if len(args)!=1:
    raise ValueError("Please supply exactly one input argument\n%s"%(usage))

#-------------------------------------------------

'''
need to load in the vector, format it into an array, and then plug it into scipy.cluster.kmeans
this should return you a set of centroids for the cluster, which you'll then need to report somehow...
    report the non-zero (or non-trivial) dimensions of each vector?
    have this automatically define OmegaScan config files for each different cluster?
        -> need to associate rows with clusters and then define OmegaScan jobs based on gps times (require these as an input list?)
'''
