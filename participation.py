usage = "python participation.py [--options] pointy.out pointy.out pointy.out ..."
description = "picks the best pvalue for each channel from multiple pointy.out files"
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
parser.add_option("-u", "--unsafe", default=False, type="string", help="a list of unsafe channels. Only used to color the plot. We require an exact match for channel to be flagged")

parser.add_option("-o", "--output-dir", default=".", type="string")
parser.add_option("-t", "--tag", default="", type="string")

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

#=================================================
### load in list of unsafe channels
if opts.unsafe:
    if opts.verbose:
        print "reading in unsafe channels from : %s"%(opts.unsafe)
    file_obj = open(opts.unsafe, "r")
    unsafe_chans = [ line.strip() for line in file_obj if line.strip()]
    file_obj.close()
else:
    unsafe_chans = []

#=================================================
### write the summary file
filename = "%s/present%s.out"%(opts.output_dir, opts.tag)
fileNAME = "%s/participating%s.out"%(opts.output_dir, opts.tag)
if opts.verbose:
    print "writing :\n\t%s\n\t%s"%(filename, fileNAME)
file_obj = open(filename, "w")
file_OBJ = open(fileNAME, "w")

### order the lists
if opts.plot:
    pvalues = []
    unsafes = []
keys = sorted(chans.keys())
keys.sort(key=lambda key: min(p for p, l, f in chans[key]) )
for key in keys:
    chans[key].sort(key=lambda l: l[0])    
    pvalue, lines, filename = chans[key][0]
    string = " ".join(lines)

    print >> file_obj, string
    if pvalue <= opts.pvalue:
        print >> file_OBJ, string
    
    if key in unsafe_chans:
        unsafes.append( pvalue )
    else:
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

    nbins = (len(pvalues) + len(unsafes))/8 ### could produce poor binning...

    bins = np.logspace( np.log10(min(np.min(pvalues), np.min(unsafes))), 0, nbins)
    ax.hist( [pvalues, unsafes], bins=bins, histtype="barstacked", color=['g', 'r'], label=['safe', 'unsafe'], log=False)

    ax.set_xlabel('min{pvalue}')
    ax.set_ylabel('count')

    ax.set_xscale("log")

    ax.grid(True, which="both")
    ax.legend(loc='best')

    ax.plot( [opts.pvalue]*2, plt.gca().get_ylim(), 'k--', linewidth=2 )

    figname = "%s/pvalue%s.png"%(opts.output_dir, opts.tag)
    if opts.verbose:
        print "\t%s"%(figname)
    fig.savefig(figname)
    plt.close(fig)

