#!/usr/bin/env python2


import xyzstats, xyzfile
import argparse, abc, functools
from contextlib import closing
import numpy as np, math
import re
import os, os.path, glob
from simw import running_avg
from smooth import smooth

import logging, sys

import numpy
numpy.seterr('raise')

verbosity = False
cut = 0

def status(s):
    global verbosity
    if verbosity:
        print(s, file=sys.stderr)

Plotters = {}
def Plotter(name, help, *args, **kw):
    def Decorator(func):
        Plotters[name] = (help, func, args, kw)
        return func
    return Decorator

def summarize(name, ts, vals, npts = None):
    vals = np.array(vals)
    mean, std = vals.mean(), vals.std()
    stdperc = 100.0 *std / mean if mean != 0 else 0
    print("%s: %.3g+-%.3g (%.3g%%) - N: %d" % (name, mean, std, stdperc, len(ts)))
     
def summarizeErr(name, vals, npts):
    import math
    vals = np.array(vals)
    mean, std = vals.mean(), vals.std()
    stdperc = 100.0 *std / mean if mean != 0 else 0
    err = std / math.sqrt(npts + 1)
    errperc = 100.0 *err / mean if mean != 0 else 0
    print("%s: %.3g+-%.3g (%.3g%%); err %.3g (%.3g%%) - N: %d" 
            % (name, mean, std, stdperc, err, errperc, npts))

def runfiles(fnames):
    import nsort
    fnames = nsort.nsort(fnames)
    for n,f in enumerate(fnames):
        base, sep, ext = f.rpartition('.')
        statfname = base + sep + 'stats'
        if len(fnames) > 1 or verbosity:
            print('Running', f, '(%d / %d)' % (n+1, len(opts.files)))
        with xyzstats.statkeeper(f, statfname, cut=cut) as sk:
            yield sk

@Plotter('tT', 'Temperature vs. Time')
@Plotter('tTg',  'Temperature vs. Time, grouped', group=True)
def TempTime(fnames, plot, group=False):
    for sk in runfiles(fnames):
        times = [t/1000.0 for t in sk.times]
        status('Getting temperatures...')
        Ts = sk.temp()
        summarize('T',times, Ts)
        if plot:
            status('plotting...')
            #~ plt.clf()
            plt.plot(times, Ts)
            if not group: plt.show()
        else:
            status('Not plotting.')
    if plot and group: plt.show()

@Plotter('TRg', 'Temp. vs. Rg')
@Plotter('TRgerr', 'Temp. vs. Rg; use Rg autocorr for Rg error bars', err=True)
@Plotter('TRgerr', 'Temp. vs. Rg; use ISF for Rg error bars', err='ISF')
def TempTime(fnames, plot, err=False):
    Allvals = []
    for sk in runfiles(fnames):
        times = [t/1000.0 for t in sk.times]
        Ts = sk.temp()
        Rgs = sk.Rg()
        Tmean, Terr = np.mean(Ts), np.std(Ts)
        summarize('T',times, Ts)
        
        if err:
            status('getting relaxation')
            relaxt = sk.relax_acorr() if err is not 'ISF' else sk.relax_ISF()
            status('getting sim_Rg')
            sim_Rg = Rgmean,Rgstd,Rgerr,npoints = sk.sim_Rg(relaxt)            
            if not npoints > 0 or not err > 0:
                summarize('Rg:',times, Rgs)
                continue
            else:
                print('npoints:', npoints, npoints > 0)
                summarizeErr('Rg:', Rgs, npoints)
        else:
            summarize('Rg',times, Rgs)
            Rgmean, Rgerr = np.mean(Rgs), np.std(Rgs)
        Allvals.append((Tmean, Terr, Rgmean, Rgerr))
    if plot:
        status('plotting...')
        #~ plt.clf()
        Trow, Terr, Rgrow, Rgerr = list(zip(*Allvals))
        plt.errorbar(Trow, Rgrow, xerr=Terr, yerr=Rgerr, fmt='b.')
        plt.xlim([0, None])
        plt.ylim([0, None])
        plt.show()
    else:
        status('Not plotting.')

@Plotter('end', 'end-to-end distance (using alpha carbons) vs. time')
def EndtoEnd(fnames, plot):
    for sk in runfiles(fnames):
        times = [t/1000.0 for t in sk.times]
        dists = sk.endtoend()
        summarize('end to end',times, dists)
        if plot:
            status('plotting...')
            plt.plot(times, dists)
            plt.show()
        else:
            status('Not plotting.')

@Plotter('endhist', 'end-to-end distance (using alpha carbons) vs. time')
def EndHist(fnames, plot):
    for sk in runfiles(fnames):
        dists = sk.endtoend()
        summarize('end to end',sk.times, dists)
        if plot:
            status('plotting...')
            plt.hist(dists, 50)
            plt.show()
        else:
            status('Not plotting.')

@Plotter('Rg', 'Radius of Gyration vs. time', ISF=False)
@Plotter('Rgisf', 'Radius of Gyration vs. time, autocorrelation from ISF', ISF=True)
def Rg(fnames, plot, ISF = True):
    for sk in runfiles(fnames):
        gyradii = sk.Rg()
        #~ Rg_avg = xyzstats.average(Rg)
        status('getting times')
        ts = [(t-sk.times[0])/1000.0 for t in sk.times]
        status('getting relaxation')
        relaxt = sk.relax_acorr() if not ISF else sk.relax_ISF()
        status('getting sim_Rg')
        sim_Rg = Rg,std,err,npoints = sk.sim_Rg(relaxt)
        if npoints is None: npoints = 0
        if np.isnan([Rg,std,err]).any() or std <= 1e-8:
            print('NANS!!!')
            continue
        acts, acs = sk.ISF_Rg() if ISF else sk.autocorr()
        acts = [t/1000.0 for t in acts]
        
        if not err > 0:
            summarize('Rg:',ts, gyradii)
        else:
            summarizeErr('Rg:', gyradii, npoints)
        if plot:
            status('plotting...')
            pairplot('Rgs from ' + sk.fname, ts, gyradii, acts, acs,
                'Time (Thousands)', 'Rg', 'Autocorrelation')
            plt.show()
        else:
            status('Not plotting.')
 
def pairplot(title, xs1, ys1, xs2, ys2, xlabel, ylabel1,ylabel2, saveto=None):
    plt.clf()
    fig=plt.gcf()
    ax1 = fig.add_subplot(111)
    ax1.plot(xs1, ys1, color='b')
    #~ ax1.plot(ts, running_avg(gyradii), color=(.6,.6,.6))
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel1, color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax2 = ax1.twinx()
    ax2.plot(xs2, ys2, color='g')
    ax2.set_ylabel(ylabel2, color='g')
    for tl in ax2.get_yticklabels():
        tl.set_color('g')
    if saveto is not None:
        plt.savefig(saveto)

if __name__ == '__main__':    
    parser = argparse.ArgumentParser()
    plotchoices = list(Plotters.keys())
    plothelp = '\n'.join(['Options:'] + [k + ': ' + v[0] for k,v in list(Plotters.items())])
    parser.add_argument('choice', type=str, choices = plotchoices,
        help=plothelp)
    parser.add_argument('files', nargs='+', type=str)
    parser.add_argument('-x', dest='plot', help='Do not plot, just gather data', action='store_false')
    parser.add_argument('-v', dest='verbose', action='store_true')
    parser.add_argument('-c', dest='cut', type=float, default=0.0,
        help='Fraction of frames to ignore (0 < cut < 1) or amount of time to ignore (cut > 1)')

    opts = parser.parse_args()
    verbosity = opts.verbose
    cut = opts.cut
        
    if len(opts.files) == 1 and os.path.isdir(opts.files[0]):
        d = opts.files[0]
        opts.files = glob.glob(d + '/*.xyz')

    if opts.plot:
        from matplotlib import pylab as plt

    funchelp, func, args, kw = Plotters[opts.choice]
    #~ print("To Plot:", opts.choice, '--', funchelp)
    #~ print("Files:", *opts.files)
    
    func(opts.files, opts.plot, *args, **kw)
