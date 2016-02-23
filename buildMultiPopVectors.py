#!/usr/bin/python
usage = "buildMultiPopVectors.py [--options] pointy.out pointy.out pointy.out ..."
description = "given a set of channels (which define a \"feature\" space of correlations), this script builds vectors in that space based on the pvalues extracted from pointy.out files supplied as arguments. The pvalues are stored as sqrt(-log10(p)) in the vectors so that Euclidean distance look like the product of probabilities"
author = "reed.essick@ligo.org"

import numpy as np

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-c", "--channel-list", default=None, type="string", help="a list of channels to use when defining the feature space")

parser.add_option("-o", "--output-filename", default=None, type="string")

opts, args = parser.parse_args()

if opts.channel_list==None:
    opts.channel_list = raw_inpu("--channel-list=")

#-------------------------------------------------

### read in channel lists
if opts.verbose:
    print "reading interesting channels from : %s"%(opts.channel_list)
file_obj = open(opts.channel_list, "r")
channels = sorted([_.strip() for _ in file_obj.readlines()])
file_obj.close()
N = len(channels)
if opts.verbose:
    print "    identified %d channels"%(N)
indmap = dict( (chan, i) for i, chan in enumerate(channels) )

#-------------------------------------------------
### write to file
if opts.output_filename:
    if opts.verbose:
        print "writing vectors into : %s"%(opts.output_filename)
    outfile = open(opts.output_filename, "w")
    errfile = open(opts.output_filename+".trgs", "w")
else:
    import sys
    outfile = sys.stdout
    errfile = sys.stderr

print >> outfile, "filename " + " ".join(channels)

### read in pointy.out files and store relevant pvalues
for pointy in args:
    vect = np.zeros((N,), dtype=float)
    if opts.verbose:
        print "processing : %s"%(pointy)

    print >> errfile, pointy

    file_obj = open(pointy, "r")
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
            if indmap.has_key(chan): ### only add a channel if it is in the selected set
                vect[indmap[chan]] = -np.log10( float(lines[nind].strip().split("=")[-1]) )
            ind = ind
        ind += 1
    vect = vect**0.5 ### we take the square root of the log of probabilities so that the Euclidean distance looks like a product of probabilities

    if opts.verbose:
        print "    found %d non-zero significances"%(np.sum(vect>0))

    print >> outfile, pointy+" "+" ".join("%.9e"%_ for _ in vect)

if opts.output_filename:
    outfile.close()
    errfile.close()
