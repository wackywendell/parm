# encoding: UTF-8
from __future__ import print_function

from sim import *
import math
import numpy as np

def geometric(tot, n):
    """Returns n integers less than (with last equal to) tot, such that
    it is as close to a geometric series as possible."""
    lst = [1]
    curval = 1
    while curval < tot and len(lst) < n:
        newval = int(.5 + curval * (float(tot) / curval) ** (1./(n-len(lst))))
        if newval == curval:
            newval += 1
        if newval > tot:
            newval = tot
        curval = newval
        lst.append(newval)
    return lst

average = np.mean

def average_squared(lst):
    l = np.array(lst)
    return np.sqrt(np.mean(l*l))

def avg(lst):
    return(sum(lst[1:], lst[0]) / float(len(lst)))

for name in list(locals().keys()):
    if '_swigregister' in name:
        del locals()[name]

def calc_Rg(residues):
    locs = [r['CA'].x for r in residues]
    center = avg(locs)
    locsqs = [(loc).sq() for loc in locs]
    return math.sqrt(average(locsqs))

def samp_std(lst):
    return np.std(lst) * math.sqrt(float(len(lst)) / (len(lst)-1))

def autocorr(lst):
    lst = np.array(lst, dtype=float)
    mu = np.mean(lst)
    var = np.var(lst)
    N = len(lst)
    def nthpair(n):
        (l1, l2) = (lst[n:], lst[:-n]) if n>0 else (lst, lst)
        l1, l2 = l1 -mu, l2-mu
        return float((N-n))/N * (l1*l2) / var
        #~ l1, l2 = l1 - l1.mean(), l2 - l2.mean()
        #~ return l1 * l2 / float(samp_std(l1) * samp_std(l2))
    lsts = [nthpair(n) for n in range(len(lst))]
    return [np.mean(l) for l in lsts]

def ISF(arrlst, scale, maxavg=None, ntimes=None):
    newarrs = [np.exp(a*(2j*math.pi/scale)) for a in arrlst]
    outlst = []
    ns = [0] + geometric(len(newarrs), ntimes)[:-1] if ntimes else range(len(newarrs))
    for n in ns:
        skip = max(1,(len(newarrs)-n) // maxavg) if maxavg is not None else 1
        pairs = zip(newarrs, newarrs[n:])[::skip]
        print(n,'pairs:',len(pairs))
        exps = [p1 / p2 for p1,p2 in pairs]
        outlst.append(abs(np.mean(exps)))
    return ns, outlst

def running_avg(lst):
    lst = np.array(lst)
    return np.cumsum(lst) / (np.arange(len(lst)) + 1)
    #~ return [np.mean(lst[:n]) for n in range(len(lst))]

def rand_walk(n, *dims):
    return np.cumsum(np.random.rand(n, *dims)-.5, 0)

class Stat:
    def __init__(self, name, func, shortname=None):
        self.name = name
        self.func = func
        self.shortname = shortname
        if not shortname:
            self.shortname = name
        self.values = []
        self.ts = []
    
    def update(self, t):
        self.ts.append(t)
        self.values.append(self.func())
    
    @staticmethod
    def plots(*stats):
        import numpy as np
        import matplotlib.pyplot as plt
        ymin = 0
        ymax = 0
        tmins = []
        tmaxs = []
        stdcol = '.75'
        labels = []
        valstart = .1
        for s in stats:
            startloc = int(len(s.values)*valstart)
            l = np.array(s.values[startloc:])
            ts = np.array(s.ts[startloc:])
            tmins.append(ts[0])
            tmaxs.append(ts[-1])
                
            mean = l.mean(0)
            std = l.std()
            sovm = 0
            if mean != 0:
                sovm = abs(std/mean)
            label = u'$\overline {{{}}} = {:.3f}$'.format(s.shortname,mean)
            stdlabel = '$\sigma={0:.3f}\;({1:.4f})$'.format(std, sovm)
            plt.axhline(mean, color=stdcol, linewidth=2)
            pstd=plt.axhline(mean-std, color=stdcol)
            plt.axhline(mean+std, color=stdcol)
            (p,)=plt.plot(ts, l, '.-')
            ymin = min(ymin, l.min())
            ymax = max(ymax, l.max())
            labels.append((p,label))
            labels.append((pstd,stdlabel))
        
        plt.axis([min(tmins), max(tmaxs), ymin, ymax ])
        plt.title(", ".join([s.name for s in stats]))
        plt.legend(*zip(*labels))
        plt.show()
    
    def plot(self):
        self.plots(self)
    
    def mean(self):
        import numpy as np
        return np.mean(self.values,0)
        
    def std(self):
        import numpy as np
        return np.std(self.values)
    
    def __str__(self):
        mean = self.mean()
        std = self.std()
        sovm = 0
        if mean != 0:
            sovm = abs(std*100.0/mean)
        return "{}={:.4g}+-{:.4g}%".format(self.name, mean, sovm)
