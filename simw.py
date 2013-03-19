# encoding: UTF-8


from sim import *
import math
import numpy as np

Matrix = Matr

def uniq(lst):
    """Returns the elements of lst in the same order, with no duplicates."""
    l = []
    for n in lst:
        if n not in l: l.append(n)
    return l

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
    if len(lst) == 1: return lst[0]
    return(sum(lst[1:], lst[0]) / float(len(lst)))

for name in list(locals().keys()):
    if '_swigregister' in name:
        del locals()[name]

def calc_Rg(residues):
    locs = [r['CA'].x for r in residues]
    center = avg(locs)
    locsqs = [(loc-center).sq() for loc in locs]
    return math.sqrt(average(locsqs))

def samp_std(lst):
    return np.std(lst) * math.sqrt(float(len(lst)) / (len(lst)-1))

def autocorr(lst):
    """Returns $<x(t)x(t+\tau)>_t$. 
    tau is the index of the output array, from 0 to len(x)-1, and
    x(t) is the normalized input: x(t) = (y(t) - \mu) / \sigma, where
        y(t) is the input to this function, and \mu and \sigma are the 
        mean and std. dev."""
    x = np.array(lst, dtype=float)
    y = (x - x.mean()) / x.std()
    return np.correlate(y, y, mode='full')[-len(x):] / len(x)
    
    # this is fully equivalent to the following:
    
    lst = np.array(lst, dtype=float)
    if len(lst) == 0: return []
    elif len(lst) == 1: return [1]
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
    return np.array([np.mean(l) for l in lsts], dtype=float)

def autocorrs(xs, ns=1000, times=None):
    if not np.iterable(ns):
        ns = np.logspace(0,np.log10(len(xs)-1), int(ns))
    ns = np.unique(np.array(ns, dtype=int))
    acs = np.array([np.corrcoef(xs[:-i], xs[i:])[0,1] for i in ns])
    ts = ns if times is None else np.array(times)[ns]
    return ts, acs

def relax(ns, times=None, acorrs=None, cutfirst = 1.0/math.e, cutlast=1.0):
    t=0
    if times is None:
        times = np.arange(len(ns))
    elif not np.iterable(times):
        times = np.arange(len(ns)) * times
    if acorrs is None:
        acorrs = np.array(autocorr(ns))
    acorr = np.array(list(zip(np.arange(len(acorrs)), acorrs)))
    try:
        tabove = (max((t for t,v in acorr if abs(v) > cutlast))
                if cutlast < max(abs(acorrs)) else 0)
        t0 = min((t for t,v in acorr if abs(v) < cutfirst and t >= tabove))
    except ValueError:
        # It never gets in the acceptable region
        return None
    t = max([times[t0] - times[0], t])
    return t

def ISF(arrlst, scale, maxavg=None, ntimes=None):
    newarrs = [np.exp(a*(2j*math.pi/scale)) for a in arrlst]
    outlst = []
    ns = [0] + geometric(len(newarrs), ntimes)[:-1] if ntimes else list(range(len(newarrs)))
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

def Rij(reslist, i, j):
    return Vec(reslist[i]['CA'].x-reslist[j]['CA'].x).mag()

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
            label = '$\overline {{{}}} = {:.3f}$'.format(s.shortname,mean)
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
        plt.legend(*list(zip(*labels)))
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

class collecStat(Stat):
    def __init__(self, collec, residues):
        Stat.__init__(self, self.name, self.func, self.shortname)
        self.residues = residues
        self.collec = collec

class Rg(collecStat):
    name = 'Radius of Gyration'
    shortname = 'Rg'
    def func(self):
        return calc_Rg(self.residues)

class Temperature(collecStat):
    name = 'Temperature'
    shortname = 'T'
    def func(self):
        return self.collec.temp()

class Energy(collecStat):
    name = 'Energy'
    shortname = 'E'
    def func(self):
        return self.collec.energy()


class Res:
    def __init__(self, anames, masses):
        self.atomvec = atomvec(masses)
        self.names = anames
    
    def __getitem__(self, name):
        if isinstance(name, int):
            idx = name
            name = self.names[idx]
        elif name not in self.names:
            raise KeyError(name)
        else:
            idx = self.names.index(name)
        a = self.atomvec[idx]
        a.name = name
        a.element = name[0]
        return a
    
    def __iter__(self):
        for n, a in zip(self.names, self.atomvec):
            a.name = n
            a.element = n[0]
            yield a
    
    def __len__(self): return len(self.names)
    
    def __contains__(self, name):
        if isinstance(name, str):
            return name in self.names
        return name in self.atomvec
    
    def get_id(self, *args):
        return self.atomvec.get_id(*args)
