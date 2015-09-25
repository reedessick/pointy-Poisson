usage = "chanlist2sections.py [--options] chanlist chanlist chanlist ... "

from optparse import OptionParser

#=================================================

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--signifmin", default=0, type="float")
parser.add_option("", "--signifmax", default=100000, type="float")

parser.add_option("", "--fmin", default=0, type="float")
parser.add_option("", "--fmax", default=16384, type="float")

parser.add_option("", "--durmin", default=0, type="float")
parser.add_option("", "--durmax", default=100, type="float")

opts, args = parser.parse_args()

#=================================================

report = """[%s]
signifmin = %.6f
signifmax = %.6f
fmin = %.6f
fmax = %.6f
durmin = %.6f
durmax = %.6f
"""%("%s", opts.signifmin, opts.signifmax, opts.fmin, opts.fmax, opts.durmin, opts.durmax)

#=================================================

channels = set()
for chanlist in args:
    if opts.verbose:
        print "reading in channels from : %s"%chanlist
    file_obj = open(chanlist, "r")
    for line in file_obj:
        channels.add( line.strip() )
    file_obj.close()

for chan in sorted(channels):
    print report%chan

