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

opts, args = parser.parse_args()

if opts.tag:
    opts.tag = "_%s"%(opts.tag)

nargs = len(args)
if not nargs:
    raise ValueError("please supply at least one pointy.out file as an input argument")

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
chans = dict( [ (key, np.prod([p for p, _, _ in chans[key]])) for key in chans.keys() ] )

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
keys = sorted(chans.keys())
for key in keys:
    pvalue = chans[key]

    print >> file_obj, "%.6e\t%s"%(pvalue, key)
    if pvalue <= opts.pvalue:
        print >> file_OBJ, key

    if opts.plot:    
        pvalues.append( pvalue )

#=================================================
### make the plot
if opts.plot:
    if opts.verbose:
        print "plotting"

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = plt.subplot(1,1,1)

    nbins = len(pvalues)/opts.nperbin ### could produce poor binning...
    this_min = min( pvalues )

    bins = np.logspace( np.log10(this_min), 0, nbins)
    ax.hist( pvalues, bins=bins, histtype="bar", color='g', log=True )

    ax.set_xlabel('prod{pvalue}')
    ax.set_ylabel('count')

    ax.set_xscale("log")

    ax.grid(True, which="both")

    ax.plot( [opts.pvalue]*2, plt.gca().get_ylim(), 'k--', linewidth=2 )

    figname = "%s/stacked-pvalue%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "\t%s"%(figname)
    fig.savefig(figname)
    plt.close(fig)

