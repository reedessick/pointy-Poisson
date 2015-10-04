usage = "chanlist2sections.py [--options] chanlist chanlist chanlist ... "

from ConfigParser import SafeConfigParser

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

parser.add_option("-c", "--config", default="pointed.ini", type="string")

parser.add_option("", "--ifo", default="H1", type="string" )
parser.add_option("", "--gdsdir", default="/gds-h1/dmt/triggers/", type="string" )
parser.add_option("", "--basename", default="H-KW_TRIGGERS", type="string" )
parser.add_option("", "--stride", default=32, type="int")

opts, args = parser.parse_args()

#=================================================

config = SafeConfigParser()

config.add_section("general")
config.set("general", "ifo", opts.ifo )
config.set("general", "conf", "0.68" )

config.add_section("kleinewelle")
config.set("kleinewelle", "gdsdir", opts.gdsdir )
config.set("kleinewelle", "basename", opts.basename )
config.set("kleinewelle", "stride", "%d"%opts.stride )

config.add_section("Omicron")
config.set("Omicron", "gdsdir", "" )
config.set("Omicron", "stride", "" )
config.set("Omicron", "channels", "" )

config.add_section("OfflineOmicron")
config.set("OfflineOmicron", "gdsdir", "" )
config.set("OfflineOmicron", "channels", "" )

#=================================================

channels = set()
for chanlist in args:
    if opts.verbose:
        print "reading in channels from : %s"%chanlist
    file_obj = open(chanlist, "r")
    for line in file_obj:
        channels.add( line.strip() )
    file_obj.close()
channels = sorted(channels)

#==============

config.set("kleinewelle", "channels", " ".join(channels) )

#==============

for chan in sorted(channels):
    config.add_section( chan )

    config.set( chan, "signifmin", "%.6f"%opts.signifmin )
    config.set( chan, "signifmax", "%.6f"%opts.signifmax )

    config.set( chan, "fmin", "%.6f"%opts.fmin )
    config.set( chan, "fmax", "%.6f"%opts.fmax )

    config.set( chan, "durmin", "%.6f"%opts.durmin )
    config.set( chan, "durmax", "%.6f"%opts.durmax )


config.write( opts.config )
