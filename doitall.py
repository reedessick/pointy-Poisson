#!/usr/bin/python
description = """an all-in-one pointy analysis tool"""
usage = """doitall.py [--options] gps gps gps..."""

import os
import glob
import subprocess as sp

from optparse import OptionParser

#=================================================

frqmap = {
        65536.0:[(8,128), (32,2048), (1024,4096), (2048,8192)],
        32768.0:[(8,128), (32,2048), (1024,4096), (2048,8192)],
        16384.0:[(8,128), (32,2048), (1024,4096), (2048,8192)],
         8192.0:[(8,128), (32,2048), (1024,4096), (2048,8192)],
         4096.0:[(8,128), (32,2048), (1024,4096)],
         2048.0:[(8,128), (32,2048)],
         1024.0:[(8,128), (32,1024)],
          512.0:[(8,128), (32,512)],
          256.0:[(8,128), (32,256)],
           16.0:[],
        }

#========================

KWconfig_header = """stride %.1f
basename %s-KW_%s_TRIGGERS
segname %s-KW_%s_SEGMENTS
significance 15.0
threshold 3.0
decimateFactor -1"""

#=================================================

def FrChannels2KWconfig( FrChannels, observatory, tag="C", kwstride=32 ):
    kwconfig = KWconfig_header%(kwstride, observatory, tag, observatory, tag)
    for line in FrChannels.split("\n"):
        if not line:
            continue
        chan, freq = line.split()
        freq = int(freq)
        for minF, maxF in frqmap[freq]:
            kwconfig += "\nchannel %s %d %d"%(chan, minF, maxF)
    return kwconfig

def KWchans2KWconfig( KWchans, observatory, tag="C", kwstride=32 ):
    kwconfig = KWconfig_header%(kwstride, observatory, tag, observatory, tag)
    for chan in KWchans:
        chan = chan.split("_")
        lF, hF = chan[-2:]
        if ":" in chan[1]: ### special case for some slow channels with non standard naming conventions
            chan = "%sE%s"%(chan[0], "_".join(chan[1:-2]))
        else:
            chan = "%s:%s"%(chan[0], "_".join(chan[1:-2]))
        kwconfig += "\nchannel %s %s %s"%(chan, lF, hF)
    return kwconfig

#=================================================

parser = OptionParser(description=description, usage=usage)

#====

parser.add_option("-v", "--verbose", default=False, action="store_true")

#==== options about data ranges

parser.add_option("-o", "--observatory", default=False, type="string")

parser.add_option("-t", "--frame-type", default="C", type="string", help="either \"C\" (default) or \"R\"")

parser.add_option("", "--whitelist", default=False, type="string", help="a file containing a list of channel names defining which channels should be analyzed. Only the channels in this list will be analyzed and we require an exact match. Channel names should be as they come out of the frames (less IFO prefix). Default is False, so all channels present will be allowed.")
parser.add_option("", "--error-if-absent", default=False, action="store_true", help="raise an error if one of the channels in --whitelist is not found in the frames")

parser.add_option("", "--blacklist", default=False, type="string", help="a file containing a list of channel names defining which channels should not be analyzed. Channel names should be as they come out of the frames (less IFO prefix). Default is False, so all channels present will be allowed.")

#==== options about KW processing and frame width

parser.add_option("-s", "--kwstride", default=32, type="int", help="desired stride for KW analysis. Should be a power of 2")
parser.add_option("-S", "--framestride", default=64, type="int", help="duration of frames. Should be a power of 2")

#==== options about pointy analysis

parser.add_option("", "--pointybin", default="/home/idq/fishyStatistics/pointy-Poisson/", type="string" )

parser.add_option("-w", "--coinc-window", default=1.0, type="float")
parser.add_option("-W", "--signif-window", default=5000.0, type="float")
parser.add_option("-e", "--exclude", default=0.0, type="float", help="passed on to pointy commands")

parser.add_option("", "--signifmin", default=[], action="append", type="float")
parser.add_option("", "--signifmax", default=1e6, type="float")
parser.add_option("", "--fmin", default=0, type="float")
parser.add_option("", "--fmax", default=16384, type="float")
parser.add_option("", "--durmin", default=0.0, type="float")
parser.add_option("", "--durmax", default=100.0, type="float")

#==== options about output and post-processing

parser.add_option("-p", "--pvalue", default=1.0, type="float")
parser.add_option("-u", "--unsafe", default=None, type="string")

parser.add_option("-O", "--output-dir", default=".", type="string")

#====

opts, args = parser.parse_args()

if not len(args):
    raise ValueError("please supply at least one gps time as an argument")
args = [float(arg) for arg in args]

if not opts.signifmin:
    opts.signifmin.append( 0.0 )

if not opts.observatory:
    opts.obervatory = raw_input("observatory = ")

cwd = os.getcwd()
if not os.path.exists(opts.output_dir):
    os.makedirs(opts.output_dir)

#=================================================

if opts.whitelist:
    if opts.verbose:
        print "loading whitelist from : %s"%(opts.whitelist)
    file_obj = open(opts.whitelist, "r")
    whitelist = [line.strip() for line in file_obj.readlines() if line.strip()]
    file_obj.close()

if opts.blacklist:
    if opts.verbose:
        print "loading blacklist from : %s"%(opts.blacklist)
    file_obj = open(opts.blacklist, "r")
    blacklist = [line.strip() for line in file_obj.readlines() if line.strip()]
    file_obj.close()

    ### check for conflict with whitelist
    if opts.whitelist:
        if opts.verbose:
           print "checking whitelist and blacklist for conflicts"
        for chan in whitelist:
            if chan in blacklist:
                raise ValueError("channel:%s in both whitelist and blacklist and is therefore ambiguous."%chan)

#=================================================

tag = opts.frame_type
frametype = "%s1_%s"%(opts.observatory, tag)
tagrds="%s-RDS"%(tag)

#=================================================

N = len(args)
for ind, gps in enumerate(args):
    output_dir = "%s/%s-%.6f"%(opts.output_dir, opts.observatory, gps)
    if not os.path.exists(output_dir):
        os.makedirs( output_dir )

    if opts.verbose:
        print "%d / %d : gps=%.3f"%(ind+1, N, gps)

    inlist = "%s/%s.frames-%d"%(output_dir, frametype, ind)
    inlistrds = "%s/%s-RDS.frames-%d"%(output_dir, frametype, ind)

    FrChannelscache = "%s/%s.FrChannels-%d"%(output_dir, frametype, ind)

    KWconf = "%s/%s-KW.cnf-%d"%(output_dir, frametype, ind)
    KWrdsconf = "%s/%s-KW-RDS.cnf-%d"%(output_dir, frametype, ind)

    coincChanlist = "%s/%s-coinc.chans"%(output_dir, frametype)

    ### go findeth data
    start = (int((gps-opts.coinc_window)/opts.framestride) - 1)*opts.framestride
    end = gps + opts.coinc_window
    cmd = "gw_data_find -o %s --type %s -u file -s %d -e %d"%(opts.observatory, frametype, start, end)
    if opts.verbose:
        print "\t", cmd
    lines = sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0].replace("file://localhost","") 

    if opts.verbose:
        print "\twriting : %s"%(inlist)
    file_obj = open( inlist, "w" )
    file_obj.write( lines )
    file_obj.close()

    ### find the channels present
    frame = lines.split()[0].strip()
    cmd = "FrChannels %s"%(frame)
    if opts.verbose:
        print "\t", cmd
    frchannels = sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0]
  
    if opts.whitelist:
        if opts.verbose:
            print "\tkeeping channels in whitelist"
        absent = False
        frchans = ""
        for line in frchannels.split("\n"):
            if not line: 
                continue
            chan, freq = line.split()
            chan = chan[3:] ### remove the prefix
            if chan in whitelist: ### channel is in the allowed set
                frchans += "%s\n"%(line)
            else:
                absent = True
        if opts.error_if_absent and absent:
            raise ValueError("at least one channel from whitelist was not found in the frame!")
        frchannels = frchans ### only keep the ones that were allowed

    if opts.blacklist:
        if opts.verbose:
            print "\tremoving channels in blacklist"
        frchans = ""
        for line in frchannels.split("\n"):
            if not line: 
                continue
            chan, freq = line.split()
            chan = chan[3:]
            if chan in blacklist:
                pass ### ignore this
            else:
                frchans += "%s\n"%(line)
        frchannels = frchans

    if opts.verbose:
        print "\tsetting up KW config : %s"%KWconf
    file_obj = open( KWconf, "w" )
    file_obj.write( FrChannels2KWconfig( frchannels, opts.observatory, tag=tag, kwstride=opts.kwstride ) )
    file_obj.close()

    ### launch KW process
    if opts.verbose:
        print "\tcd %s"%output_dir
    os.chdir(output_dir)

    cmd = "kleineWelleM %s -inlist %s"%(KWconf, inlist)
    out = "%s/%s.out"%(output_dir, frametype)
    err = "%s/%s.err"%(output_dir, frametype)
    if opts.verbose:
        print "\t", cmd
        print "\t\tout : %s"%(out)
        print "\t\terr : %s"%(err)
    out_obj = open( out, "w" )
    err_obj = open( err, "w" )
    proc = sp.Popen(cmd.split(), stdout=out_obj, stderr=err_obj)
    out_obj.close()
    err_obj.close()
    proc.wait()
    if opts.verbose:
        print "\tcd %s"%cwd
    os.chdir(cwd)

    ### scrape trg files to get coincident triggers
    if opts.verbose:
        print "\tfinding %s-KW_%s_TRIGGERS files"%(opts.observatory, tag)
    trgfiles = sorted(glob.glob("%s/%s-KW_%s_TRIGGERS-*/%s-KW_%s_TRIGGERS-*-%d.trg"%(output_dir, opts.observatory, tag, opts.observatory, tag, opts.kwstride) ))
    if opts.verbose:
        print "\t\tfound %d files\n\tsearching for triggers coincident to within %.3f sec"%(len(trgfiles), opts.coinc_window)
    chans = set()
    coinc_start = gps-opts.coinc_window
    coinc_stop = gps+opts.coinc_window
    for trgfile in trgfiles:
        file_obj = open(trgfile, "r" )
        for line in file_obj:
            if line[0] != "#":
                fields = line.strip().split()
                tstart = float(fields[0])
                tend = float(fields[1])
                chan = fields[-1]
                if (coinc_start <= tstart) and (tstart <= coinc_stop): ### at least some overlap with time-of-interest
                    chans.add( chan )
    chans = sorted( chans )
    if opts.verbose:
        print "\t\tfound %s coincident channels\n\tsetting up RDS KWconfig : %s"%(len(chans), KWrdsconf)
    file_obj = open(KWrdsconf, "w")
    file_obj.write( KWchans2KWconfig( chans, opts.observatory, tag=tagrds, kwstride=opts.kwstride ) )
    file_obj.close()

    if opts.verbose:
        print "\twriting : %s"%(coincChanlist)
    file_obj = open( coincChanlist, "w" )
    for chan in chans:
        print >> file_obj, chan
    file_obj.close()

    ### compute new KW config for only coincident triggers
    start = gps - opts.signif_window
    end = gps + opts.signif_window
    cmd = "gw_data_find -o %s --type %s -u file -s %d -e %d"%(opts.observatory, frametype, start, end)
    if opts.verbose:
        print "\t", cmd
    lines = sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0].replace("file://localhost","")

    if opts.verbose:
        print "\twriting : %s"%(inlistrds)
    file_obj = open( inlistrds, "w" )
    file_obj.write( lines )
    file_obj.close()

    ### launch KW process
    if opts.verbose:
        print "\tcd %s"%output_dir
    os.chdir(output_dir)

    cmd = "kleineWelleM %s -inlist %s"%(KWrdsconf, inlistrds)
    out = "%s/%s-RDS.out"%(output_dir, frametype)
    err = "%s/%s-RDS.err"%(output_dir, frametype)
    if opts.verbose:
        print "\t", cmd
        print "\t\tout : %s"%(out)
        print "\t\terr : %s"%(err)
    out_obj = open( out, "w" )
    err_obj = open( err, "w" )
    proc = sp.Popen(cmd.split(), stdout=out_obj, stderr=err_obj)
    out_obj.close()
    err_obj.close()
    proc.wait()
    if opts.verbose:
        print "\tcd %s"%cwd
    os.chdir(cwd)

    ### format data so that pointy call works
    gdsdir = "%s/%s-KW_%s_TRIGGERS"%(output_dir, opts.observatory, tagrds)
    if not os.path.exists( gdsdir ):
        os.makedirs( gdsdir )

    cmd = "mv %s-*/ %s/"%(gdsdir, gdsdir)   
    os.system( cmd )

    ### set up pointy config file
    procs = []
    outs = []
    for signifmin in opts.signifmin:

        pointyconf = "%s/%s-RDS-%d.ini"%(output_dir, frametype, signifmin)

        cmd = "python %s/chanlist2config.py --ifo %s1 --gdsdir %s --basename %s --stride %d --signifmin %f --signifmax %f --fmin %f --fmax %f --durmin %f --durmax %f -c %s %s"%(opts.pointybin, opts.observatory, output_dir, "%s-KW_%s_TRIGGERS"%(opts.observatory, tagrds), opts.kwstride, signifmin, opts.signifmax, opts.fmin, opts.fmax, opts.durmin, opts.durmax, pointyconf, coincChanlist)
        if opts.verbose:
            print "\t", cmd
        sp.Popen(cmd.split()).wait()

        ### set up pointy command
        cmd = "python %s/pointed.py -c %s -w %.6f -e %.6f %.6f"%(opts.pointybin, pointyconf, opts.signif_window, opts.exclude, gps)
        out = "%s/%s-%s-pointy-%d.out"%(output_dir, opts.observatory, tagrds, signifmin)
        err = "%s/%s-%s-pointy-%d.err"%(output_dir, opts.observatory, tagrds, signifmin)
        outs.append( out )
        if opts.verbose:
            print "\t", cmd
            print "\t\tout : %s"%(out)
            print "\t\terr : %s"%(err)
        out_obj = open( out, "w" )
        err_obj = open( err, "w" )
        proc = sp.Popen(cmd.split(), stdout=out_obj, stderr=err_obj )
        out_obj.close()
        err_obj.close()

        procs.append( proc )

    for proc in procs:
        proc.wait()

    ### define participation job
    cmd = "python %s/participation.py -p %.6e -P -o %s %s"%(opts.pointybin, opts.pvalue, output_dir, " ".join(outs))
    if opts.unsafe:
        cmd += " -u %s"%(opts.unsafe)
    if opts.verbose:
        print "\t", cmd
    sp.Popen(cmd.split()).wait()


