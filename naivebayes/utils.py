__doc__ = "a module for utilisties associated with the naive bayes classifier based on the point statistic"
__author__ = "reed.essick@ligo.org"

#-------------------------------------------------

import os
import subprocess as sp

import numpy as np

#-------------------------------------------------

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

KWconfig_header = """stride %.1f
basename %s-KW_%s_TRIGGERS
segname %s-KW_%s_SEGMENTS
significance 15.0
threshold 3.0
decimateFactor -1"""

#------------------------

def chanlist2KWconfig( chanlist, observatory, tag="C", kwstride=32 ):
    """
    generate a KW config
    """
    kwconfig = KWconfig_header%(kwstride, observatory, tag, observatory, tag)
    for chan, freq in chanlist:
        freq = int(freq)
        for minF, maxF in frqmap[freq]:
            kwconfig += "\nchannel %s %d %d"%(chan, minF, maxF)
    return kwconfig

def seg2inlist(start, end, frametype, observatory, verbose=False)
    """
    go find data and generate an inlist for KW trigger production
    """
    cmd = "gw_data_find -o %s --type %s -u file -s %d -e %d"%(observatory, frametype, start, stop) ### add some buffer to be safe

    if verbose:
        print( cmd )
    return sp.Popen(cmd.split(), stdout=sp.PIPE).communicate()[0].replace("file://localhost","").split('\n')

def genkwtrgs(configpath, inlist, frametype, directory='.', verbose=False):
    """
    schedule and manage KW trig production
    """
    cwd = os.getcwd()
    cmd = "kleineWelleM %s -inlist %s"%(kwconfig_path, inlist_path)
    out = "%s/%s.out"%(directory, frametype)
    err = "%s/%s.err"%(directory, frametype)
    if verbose:
        print( "    cd %s"%directory )
        print "    "+cmd
        print "    out : "+out
        print "    err : "+err
    out_obj = open( out, "w" )
    err_obj = open( err, "w" )
    os.chdir(directory)
    proc = sp.Popen(cmd.split(), stdout=out_obj, stderr=err_obj)
    out_obj.close()
    err_obj.close()
    proc.wait()

    if verbose:
        print "    cd "+cwd
    os.chdir(cwd)

#-------------------------------------------------

def chanlist2chans(path):
    with open(chanlist_path, 'r') as obj:
        chanlist = [line.strip().split() for line in obj.readlines()]
    return chanlist

#------------------------

def path2pdict(path, exclude=[]):
    """
    parse the pointy.out files into (channel, pvalue) pairs
    """
    pvalues = []
    with open(path, 'r') as file_obj:

        line = out_obj.readline()
        while line:
            if "channel" in line:
                chan = line.strip().split("=")[-1]
                if chan not in exclude:
                    for _ in xrange(8):
                        line = out_obj.readline() ### skip down to the interesting line
                    pvalues.append( (chan, float(line.strip().split('=')[1])) )
            line = out_obj.readline()
    return dict(pvalues)

#------------------------

def paths2trgs(paths, thresholds, exclude=[], start=-np.infty, stop=+np.infty, verbose=False):
    """
    read in KW trg files into the appropriate structures
    """
    trgs = dict()
    for trg_path in paths:
        if verbose:
            print('reading: '+trg_path)
        with open(trg_path, 'r') as obj:
            cols = obj.readline().strip()[1:].split() ### FIXME may be fragile...
            toFloat = [_ for _ in cols if _!='channel']
            for line in obj:
                fields = dict(zip(cols, line.strip().split()))

                if fields['channel'] not in exclude:
                    if not trgs.has_key(fields['channel']):
                        trgs[fields['channel']] = dict((thr, []) for thr in opts.kwsignif_thr) ### set up this KW channel

                    for key in toFloat:
                        fields[key] = float(fields[key])
                    ### bin data by kwsignif_thr
                    for thr in opts.kwsignif_thr:
                        if fields['significance'] >= thr:
                            trgs[fields.pop('channel')][thr].append( (fields['time'], fields['stop_time']-fields['start_time']) ) ### add it only to one set 
                            break ### do not double count
                    else:
                        pass ### only remember things that we'll count later on...
                             ### because of the default list for opts.kwsignif_thr, we should retain everything, though

    # map into structured arrays for easy access later
    # filter triggers by gpsstart, gpsstop
    dtype=[('time', 'float'), ('duration', 'float')]
    for chan, val in trgs.items():
        for thr, times in val.items():
            times = np.array(times, dtype=dtype)
            t = times['time']
            val[thr] = times[(gpsstart<=t)*(t<gpsstop)] ### keep only the triggers within the requested window
    return trgs

def trgs2timeseries(gps, trgs, window, thresholds, reference_gps=None, verbose=False, plot_start_dur=[], save=True, directory='.', tag=''):
    """
    generate p-value timeseries from triggers
    """
    dur = len(gps)*(gps[1]-gps[0]) ### FIXME, possibly fragile?

    # generate gps time samples for reference
    if reference_gps!=None:
        tref = int(reference_gps)
    else:
        tref = int(gps[0])
    N = len(gps)

    if plot_start_dur:
        plotting_index = max(1, int(N/10000)) ### only plot 1 out of every "plotting_index" points to avoid breaking matplotlib

    ### get boundaries for rate estimation window
    minima = gps-window
    minima[minima<start] = start
    maxima = gps+window
    maxima[maxima>stop] = stop

    norms = maxima-minima ### normalization for rate estimate
                          ### we divide by this, but should be pretty confident no element==0 based on reasonable input parameters...
    if np.any(norms==0):
        raise ValueError, 'some values of norms==0, so you will want to re-think your arguments'

    ### temporary helper arrays to speed things up by avoiding reallocating memory
    tmp = np.empty((2, N), dtype=float)
    lpvl = np.zeros(N, dtype=float) ### a holder for log(pval) as a function of time
    counts = np.zeros(N, dtype=float)
    dts = np.empty(N, dtype='float') ### a holder used later
    truth = np.empty(N, dtype='bool') ### same here. We're trying hard to avoid re-allocating memory within the loop
    inds = np.arange(N) ### used when setting minimum duration

    # compute the biggest windows we'd ever have
    tmp[0,:] = maxima-gps
    tmp[1,:] = gps-minima
    starting_dts = np.max(tmp, axis=0) ### the biggest window we have at each gps time...

    nbays = np.zeros_like(gps, float) ### holder to the total naive bayes statistic

    ###----------- rate pointy-Poisson p-value ESTIMATION
    if verbose:
        print( "estimating p-values of aux triggers" )

    for chan, data in trgs.items():
        if verbose:
            print( "    "+chan )

        lpvl[:] = 0. ### reset this to zero
        counts[:] = 0.
        dts[:] = starting_dts[:]

        ### iterate over KWsignif thr and extremize lnpvl
        for thr in thresholds:
            times = data[thr]
            N = len(times)
            if N: ### find the minimum dt, otherwise stick with starting_dts
                ### find the mid-points between triggers and then bracket the whole thing with +/-infty
                bounds = [-np.infty] + list(0.5*(times['time'][1:]+times['time'][:-1])) + [np.infty]

                for m, M, trg in zip(bounds[:-1], bounds[1:], times):
                    t = trg['time']     ### extract crap
                    counts += (minima<=t)*(t<=maxima) ### add to rate estimate if within the windows

                    ### NOTE: 
                    ###     while we iterate directly over triggers, we only analyze the few gps samples relevant for that trigger instead of the entire array
                    truth[:] = (m<=gps)*(gps<M) ### select off only the times relevant for this triggers
                    if np.any(truth): ### at times, neighboring triggers can be so close there are not relevant samples between them...
                                      ### only perform complicated minimum logic when we have something meaningful to do
                        min_dt = opts.frac_of_dur*trg['duration']
                        tmp[0,truth] = np.abs(gps[truth]-t)

                        tmp[1,truth] = min_dt
                        tmp[0,truth] = np.max(tmp[:,truth], axis=0)
#                        tmp[0,inds[tmp[0,truth]<min_dt]] = min_dt ### FIXME this doesn't work for some reason...

                        tmp[1,truth] = dts[truth]
                        dts[truth] = np.min(tmp[:,truth], axis=0)

            truth[:] = counts > 0 ### FIXME: only update things for which we have a non-zero rate estimate...
                                  ### this is not great, but it's the best I can do for now. 
                                  ### Should really marginalize over rates, which will completely remove this issue
            tmp[0,truth] = lpvl[truth]
            tmp[1,truth] = np.log10(1 - np.exp(-2*dts[truth]*counts[truth]/norms[truth])) ### still occassionally getting -infty here...
            lpvl[truth] = np.min(tmp[:,truth], axis=0) ### minimize at each gps time

        ### write time-series to disk
        if save:
            path = "%s/%s%s-%d-%d.npy"%(directory, chan, tag, start, dur)
            if verbose:
                print( "      "+path )
            np.save(path, lpvl)

        ### plot time-series
        if plot_start_dur:
            fig = plt.figure()
            ax = fig.gca()

            ax.plot((gps-tref)[::plotting_index], lpvl[::plotting_index])

            ax.set_xlabel( '$t-%d$'%tref )
            ax.set_ylabel( '$\log_{10} p_i$' )
            ax.grid(True, which='both')
            ax.set_title(chan.replace('_','\_'))

            if reference_gps!=None:
                ax.plot(ax.get_xlim(), [np.interp(reference_gps, gps, lpvl)]*2, 'k--')
                ax.plot([reference_gps-tref]*2, ax.get_ylim(), 'k--')

            for start, dur in plot_start_dur:
                ax.set_xlim(xmin=start-tref, xmax=stop-tref)
                path = "%s/%s%s-%d-%d.png"%(directory, chan, tag, start, dur)
                if verbose:
                    print( "      "+path )
                fig.savefig(path)
            plt.close(fig)

        ### add this to the running total
        nbays += lpvl

    return nbays
