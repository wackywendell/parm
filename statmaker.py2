from __future__ import print_function
from xyzfile import XYZreader
import gzip, os, os.path
from numpy import array, mean, loadtxt, genfromtxt, savetxt
import math, numpy
from simw import Vec, autocorr, calc_Rg, geometric, average, average_squared
import simw as sim
import simpdb
import logging
from functools import wraps
import csv

mydir = os.path.expanduser('~/idp/')
pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'

class NotCalculableError(Exception): pass

class Statmaker:
    _statdict={}
    def __init__(self, cut=None, H=False, memoize=True,
                pdbfname = pdbfile, loadfname = loadfile):
        self.cut = cut
        self.H = H
        self.memoize = memoize
        self.pdbfile = pdbfname
        self.loadfile = loadfname
        for k,stat in self._statdict.items():
            if hasattr(self, k):
                raise ValueError("Already filled slot " + str(k))
            setattr(self, k, stat(self))
    
    def _cut(self, lst, ts=None):
        if 0 < self.cut < 1:
            numcut = int(len(lst) * self.cut)
            return lst[numcut:]
        elif self.cut >= 1:
            raise ValueError("cut %r larger than 1" % self.cut)
        return lst
    
    @classmethod
    def register(cls, name, stat):
        if name in cls._statdict and cls._statdict[name] != stat:
            raise ValueError("cls._statdict cannot register " + repr(stat) + 
                    ", key " + repr(name) + " already taken by " +
                    repr(cls._statdict[name]))
        cls._statdict[name] = stat
    
    @property
    def residues(self):
        if not hasattr(self, '_residues'):
            logging.info('building residues')
            self._residues = simpdb.Resvec.from_pdb('aS', 
                    self.pdbfile, loadfile, H=self.H, numchains=1)
            #~ atomn = len(f0.locarray)
            #~ #reslens = list(np.cumsum([len(r) for r in self._residues]))
            #~ reslens = list(np.cumsum([len(r) for r in self._residues]))
            #~ if atomn not in reslens:
                #~ raise ValueError('Could not find %d in reslens: %s' % (atomn, reslens))
            #~ indx = reslens.index(atomn)+1
            #~ if indx <= len(self._residues):
                #~ self._residues = self._residues[:indx]
                #~ if hasattr(self, '_atoms'): del self._atoms
            logging.info('residues built')
        return self._residues
    
    @property
    def atoms(self):
        if not hasattr(self, '_atoms'):
            self._atoms = [a for r in self.residues for a in r]
        return self._atoms
        

class Statistic:
    class __metaclass__(type):
        def __new__(mcs, name, bases, cdict):
            cls = type.__new__(mcs, name, bases, cdict)
            # don't register this class or others starting with 'Stat'
            if name[:4] == 'Stat': return cls
            if name[0] == '_': return cls
            if not hasattr(cls, name):
                cls.name = name
            Statmaker.register(cls.name, cls)
            return cls
    
    def __init__(self, statmaker, memoize='auto'):
        self.maker = statmaker
        self.memoize = statmaker.memoize if memoize == 'auto' else memoize
        self.memodict = {}
    
    @property
    def residues(self):
        return self.maker.residues
    
    @property
    def atoms(self):
        return self.maker.atoms
    
    @wraps(Statmaker._cut)
    def _cut(self, *args, **kw):
        return self.maker._cut(*args, **kw)
    
    @property
    def H(self):
        return self.maker.H
    
    def fname(self, xyzfilename, *args, **kw):
        direc, xyzbase = os.path.split(xyzfilename)
        bname, dot, ext = xyzbase.rpartition('.')
        addbit = self.addtofname(*args, **kw)
        return os.path.join(direc, bname + '-' + addbit + '.txt.gz')
    
    def addtofname(self, *args, **kw):
        return (self.name + '-'.join(str(a) for a in args)
                + '-'.join([k + '-' + str(v) for k,v in sorted(kw.items())]))
    
    def _memoize(self, k, v):
        if self.memoize:
            self.memodict[k] = v
        return v
    
    def _hash_args(self, xyzfilename, *args, **kw):
        return (xyzfilename, args, tuple(sorted(kw.items())))
    
    def __call__(self, xyzfilename, *args, **kw):
        if not isinstance(xyzfilename, str):
            xyzfilename = str(xyzfilename)
        memokey = self._hash_args(xyzfilename, *args, **kw)
        if self.memoize and memokey in self.memodict:
            return self.memodict[memokey]
        return self._memoize(memokey, 
                    self.run(xyzfilename, *args, **kw))
        
        
    def run(self, xyzfilename, *args, **kw):
        fname = self.fname(xyzfilename, *args, **kw)
        
        failed=[]
        
        if fname is None:
            failed.append('Statfile given as "None"')
        elif not os.path.exists(fname):
            failed.append('Statfile ' + repr(fname) + ' does not exist')
        elif (os.path.exists(xyzfilename) and 
                    os.stat(xyzfilename).st_mtime > os.stat(fname).st_mtime):
            failed.append('Statfile too young')
            os.remove(fname)
        elif os.path.getsize(fname) > 0:
            with gzip.open(fname, 'r') as f:
                f.seek(1)
                if f.tell() == 1:
                    f.seek(0)
                    try:
                        return self._from_file(f, *args, **kw)
                    except NotCalculableError:
                        failed.append('From Statfile ' + repr(fname))
                else:
                    failed.append('Statfile ' + repr(fname) + ' is empty')
                    
        else:
            failed.append('Statfile ' + repr(fname) + ' does not exist')
        
        #print(*failed)
        if fname is None:
            try:
                return self._derived(xyzfilename, None, *args, **kw)
            except NotCalculableError:
                failed.append('Derived' + repr(xyzfilename))
        
        else:
            if self.readoutfile and os.path.exists(fname):
                with gzip.open(fname, 'r') as outf:
                    outf.seek(1)
                    if outf.tell() == 1:
                        outf.seek(0)
                        readoutdat = self._read_dat(outf, *args, **kw)
                    else:
                        readoutdat = None
                args = (readoutdat,) + args
            elif self.readoutfile:
                readoutdat = None
                args = (readoutdat,) + args
            else: readoutdat = None
            
            with gzip.open(fname, 'w') as outf:
                try:
                    return self._derived(xyzfilename, outf, *args, **kw)
                except NotCalculableError:
                    failed.append('Derived' + repr(xyzfilename))
                
                if not os.path.exists(xyzfilename):
                    failed.append('xyz file %s does not exist' % xyzfilename)
                else:
                    with open(xyzfilename, 'r') as f:
                        xyzreader = XYZreader(f)
                        try:
                            vs = (v for v in xyzreader.vdicts())
                            return self._from_vdicts(vs, outf, *args, **kw)
                        except NotCalculableError:
                            failed.append('From vdicts ' + repr(xyzfilename))
                            f.seek(0)
                        
                        try:
                            return self._from_frames(iter(xyzreader), outf, *args, **kw)
                        except NotCalculableError:
                            failed.append('From frames ' + repr(xyzfilename))
                            f.seek(0)
                        
                        try:
                            return self._from_atoms(
                                    xyzreader.into(self.atoms, False),
                                    outf, *args, **kw)
                        except NotCalculableError:
                            failed.append('From atoms ' + repr(xyzfilename))
            
            raise ValueError('Tried ' + str(len(failed)) + 'methods: ' + 
                        '; '.join(failed))
    
    
    def _derived(self, xyzf, outf, *args, **kw):
        raise NotCalculableError
    
    def _from_frames(self, frames, outf, *args, **kw):
        raise NotCalculableError
    
    def _from_atoms(self, frameiter, outf, *args, **kw):
        raise NotCalculableError
    
    def _from_vdicts(self, vs, outf, *args, **kw):
        raise NotCalculableError
    
    readoutfile = False # whether to read the outf before calculating
                        # will be placed in the first arg, comes after outf
                        #e.g. _derived(self, outf, outfdat, *args, **kw)
    
    def _read_dat(self, outf, *args, **kw):
        raise NotImplementedError
    
    def _from_file(self, f, *args, **kw):
        raise NotCalculableError

########################################################################
########################################################################
########################################################################

def _readtsv_old(f, dtype=None):
    r = csv.reader(f,delimiter='\t')
    try:
        ks, vals = next(r), iter(r)
    except StopIteration:
        # empty / less than one line file
        return dict()
    vals = array([array(line, dtype=int) for line in vals if line], dtype=int)
    
    if dtype is not None:
        return (dict(zip(ks, numpy.array(list(vals),dtype=dtype).T))
                    if dtype is not None
                    else dict(zip(ks, zip(*vals))))
    
def _readtsv(f, dtype=None):
    d = (
        genfromtxt(f, dtype=dtype, names=True) if dtype is not None
        else genfromtxt(f, names=True)
        )
    return dict([(n,d[n]) for n in d.dtype.names])

def _writetsv(f, d, fmt):
    w = csv.writer(f,delimiter='\t')
    ks, cols = zip(*sorted(d.items()))
    w.writerow(map(str, ks))
    if all(isinstance(c, numpy.ndarray) for c in cols):
        rows = array(cols).T
        delimiter='\t'
        savetxt(f, rows, fmt=fmt, delimiter='\t')
        return
    rs = [[fmt % n for n in r] for r in zip(*cols)]
    w.writerows(rs)

class StatRunner(Statistic):
    """A Running statistics class."""
    def _from_file(self, f, *args, **kw):
        arr = self._cut(array(
                [map(float, l.strip().split('\t')) for l in f]
            ))
        shape = arr.shape
        if len(shape) == 2 and 1 in shape: return arr.flatten()
        return arr
        
    floatformat = '%.6g'
    
    def _from_vdicts(self, vs, outf, *args, **kw):
        vdictname = self.vdictname if hasattr(self, 'vdictname') else self.name
        if vdictname is None:
            raise NotCalculableError
        try:
            vals = array([float(v[vdictname]) for v in vs])
        except KeyError:
            raise NotCalculableError
        if outf is not None:
            savetxt(outf, vals.T, fmt=self.floatformat, delimiter='\t')
        return self._cut(vals)

class StatTable(StatRunner):
    vdictname=None
    readoutfile=True
    dtype=None
    
    def addtofname(self, *args, **kw): return self.name
    
    def keys(self, *args, **kw):
        raise NotImplementedError
    
    def _read_dat(self, f, *args, **kw):
        fkeys = f.readline().strip().split('\t')
        if len(fkeys) == 1 and fkeys[0] is '': return None
        arr = loadtxt(f, dtype=self.dtype)
        if len(fkeys) == 1 and arr.ndim == 1:
            return {fkeys[0]:arr}
        return dict(zip(fkeys,arr.T))
        
    def _from_file(self, f, *args, **kw):
        fkeys = f.readline().strip().split('\t')
        vkeys, strkeys = self.keys(*args, **kw)
        if not all(k in fkeys for k in strkeys): raise NotCalculableError
        
        arr = loadtxt(f, dtype=self.dtype)
        if len(fkeys) == 1 and arr.ndim == 1:
            d = {fkeys[0]:arr}
        else: d = dict(zip(fkeys,arr.T))
        
        cutd = dict([
                    (vkey, self._cut(d[skey]))
                    for vkey, skey in zip(vkeys, strkeys)
                ])
        return cutd
    
    def _from_frame(self, f, outfdat, *args, **kw):
        raise NotImplementedError
    
    def _from_atoms(self, frames, outf, outfdat, *args, **kw):
        vkeys, strkeys = self.keys(*args, **kw)
        valsbyframe = [self._from_frame(f, outfdat, *args, **kw) for f in frames]
        vals = zip(*valsbyframe)
        outd = dict(zip(strkeys,vals))
        firstoutd = dict() if outfdat is None else outfdat
        firstoutd.update(outd)
        _writetsv(outf,firstoutd, self.floatformat)
        return dict(zip(vkeys, self._cut(vals)))

class T(StatRunner): pass
class E(StatRunner): pass
class Times(StatRunner): vdictname='time'
class Rg(StatRunner): pass
class relax(Statistic):
    """Relaxation time, calculated using the autocorrelation function 
    of the Rg timeseries.
    
    The relaxation time is calculated as the first time at which the
    autocorrelation function is below cutfirst and stays below cutlast.
    """
    def addtofname(self, cutfirst = 1.0/math.e, cutlast=1.0):
        return "relax-%.2f-%.2f" % (cutfirst, cutlast)
    
    def _from_file(self, f, cutfirst = 1.0/math.e, cutlast=1.0):
        return float(f.readline().strip())
    
    def _derived(self, xyzf, outf, cutfirst = 1.0/math.e, cutlast=1.0):
        Rgs = self.maker.Rg(xyzf)
        Times = self.maker.Times(xyzf)
        acorrs = numpy.array(sim.autocorr(Rgs))
        acorr = numpy.array(zip(numpy.arange(len(acorrs)), acorrs))
        try:
            tabove = (max((t for t,v in acorr if abs(v) > cutlast))
                    if cutlast < max(abs(acorrs)) else 0)
            t0 = min((t for t,v in acorr if abs(v) < cutfirst and t >= tabove))
        except ValueError:
            # It never gets in the acceptable region
            ts1 = [t for t,v in acorr if abs(v) > cutlast]
            #~ print(ts1[:3], ts1[-3:])
            #if ts1:
                #ts2 = [t for t,v in acorr if v < cutfirst and t >= max(ts1)]
                #print("ERROR in statmaker.relax.run", max(ts1), ts2[:3], ts2[-3:])
            return None
        t = Times[t0] - Times[0]
        print(t, file=outf)
        return t

class MonomerDist(StatRunner):
    """Average distance between monomers, calculated from COMs"""
    vdictname=None
    floatformat='%.4g'
    
    def _from_frame(self):
        COMs = [r.getCOM() for r in self.residues]
        return array([(c2 - c1).mag() for c1,c2 in zip(COMs[:-1], COMs[1:])], dtype=float)
    
    def _from_atoms(self, frames, outf):
        vals = array([self._from_frame() for f in frames], dtype=float)
        if outf is not None:
            savetxt(outf, vals, fmt=self.floatformat, delimiter='\t')
        return self._cut(vals)

class MonomerDist2(MonomerDist):
    """Average distance between monomers, calculated from center of RGs"""
    vdictname=None
    
    def _from_frame(self):
        locs = [r.getCORg() for r in self.residues]
        return array([(c2 - c1).mag() for c1,c2 in zip(locs[:-1], locs[1:])], dtype=float)

class MonomerDists(MonomerDist):
    vdictname=None
    floatformat='%.4g'
    
    def _from_frame(self):
        COMs = [r.getCOM() for r in self.residues]
        
        return array([(c2 - c1).mag() 
                for n,c1 in enumerate(COMs)
                for c2 in COMs[n+2:]
                ], dtype=float)

class relaxN(Statistic):
    def fname(self, xyzf, *args, **kw): return None
    def _derived(self, xyzf, outf, cutfirst = 1.0/math.e, cutlast=1.0):
        relax = self.maker.relax(xyzf, cutfirst, cutlast)
        Times = self.maker.Times(xyzf)
        return 0.0 if relax is None else float(Times[-1]) / float(relax) 

class Rijs(StatTable):
    vdictname=None
    dtype=float
    floatformat='%.3f'
    
    def addtofname(self, ijs): return 'Rijs'
    
    def keys(self, ijs):
        return list(ijs), ['%d-%d' % ij for ij in ijs]
    
    def _hash_args(self, xyzfilename, ijs):
        return (xyzfilename, tuple(map(tuple, ijs)))
    
    def _from_frame(self, frame, outfdat, ijs):
        apairs = [(self.residues[i-1]['CA'],self.residues[j-1]['CA']) for (i,j) in ijs]
        return array([(a1.x - a2.x).mag() for a1,a2 in apairs], dtype=float)
        

class Rijext(StatTable):
    vdictname=None
    dtype=float
    floatformat='%.3f'
    
    def addtofname(self, ijs, dist): return 'Rijext-%.2f' % dist
    
    def keys(self, ijs, dist):
        return list(ijs), ['%d-%d' % ij for ij in ijs]
    
    def _hash_args(self, xyzfilename, ijs, dist):
        return (xyzfilename, tuple(map(tuple, ijs)), dist)
    
    @staticmethod
    def getLoc(res, dist):
        vCOM = res.getCOM()
        vdist = (vCOM - res['CA'].x).norm() * dist
        return vCOM + (vdist * dist)
            
    
    def _from_frame(self, frame, outfdat, ijs, dist):
        rpairs = [(self.residues[i-1],self.residues[j-1]) for (i,j) in ijs]
        return array([(self.getLoc(r1,dist) - self.getLoc(r2, dist)).mag() 
                    for r1,r2 in rpairs], dtype=float)    

class Dihedral(StatRunner):
    vdictname=None
    floatformat='%.4f'
    def addtofname(self, ang):
        return 'dihedral-%s' % ang
    
    def _single_from_residues(self, ang):
        res = self.residues
        restriplets = zip(res, res[1:], res[2:])
        return array([res.dihedral(ang, p, n) for p,res,n in restriplets])
    
    def _from_atoms(self, frames, outf, ang):
        vals = array([self._single_from_residues(ang) for f in frames], dtype=float)
        if outf is not None:
            savetxt(outf, vals, fmt=self.floatformat, delimiter='\t')
        return self._cut(vals)

default_ijs = [(9,130), (33,130), (54,130), (72,130), (33,72), (9,54),
        (72,92), (54,72), (9,72), (9,33), (54,92), (92,130)]
        
class ETeffsAS(StatRunner):
    def addtofname(self, R0):
        return 'ETeffs-%.2f' % R0
        
    def _from_file(self, f, R0):
        lines = [l.strip().split('\t') for l in f]
        return dict([((int(i),int(j)),float(v)) for i,j,v in lines])
    
    @staticmethod
    def from_Rijs(Rij, R0):
        Rijs = array(Rij, dtype=float)
        return mean(1.0/(1.0 + (Rijs / R0)**6))
    
    def _derived(self, xyzf, outf, R0):
        rijd = self.maker.Rijs(xyzf, default_ijs)
        d = dict([(ij,self.from_Rijs(rijs, R0)) for (ij,rijs) in rijd.items()])
        for ij,v in d.items():
            i,j = ij
            print(i, j, v, sep='\t', file=outf)
        return d

########################################################################
########################################################################
########################################################################

def filefinder(dir,  *names, **args):
    """
    extra kw:
    regexp
    matchall=False
    """
    import re, fpath
    from decimal import Decimal
    from namespace import Namespace
    matchall = args.get('matchall',True)
    if 'matchall' in args: del args['matchall']
    
    dir = fpath.Dir(dir)
    xchildren = [f for f in dir.children() if f.extension == 'xyz' or f[-1][-7:] == '.txt.gz']
    
    regexpr = (args['regexp'] if 'regexp' in args else
    '-'.join([n + '((?:[0-9]*\.?[0-9]*)|(?:inf))' for n in names]))
    if 'regexp' in args: del args['regexp']
    regex = re.compile(regexpr)
    matches = [(list(regex.finditer(f[-1])), f) for f in xchildren]
    badmatches = [f for ms,f in matches if len(ms) == 0]
    if badmatches and matchall:
        raise ValueError('%s could not match %s' % (regexpr, badmatches[0][-1]))
    matches = [(ms[0],f) for ms,f in matches if len(ms) > 0]
    filetomatches = {}
    for m, f in list(matches):
        if None in m.groups():
            print('Groups:', m.groups())
            raise ValueError('%s could not entirely match %s' % (regexpr, f[-1]))
        elif m.start() != 0:
            print('Groups:', m.groups())
            raise ValueError('%s matched %s at the wrong place' % (regexpr, f[-1]))
        if f.extension == 'xyz':
            filetomatches[f] = m
        else:
            xyzf = fpath.File(f[:-1] + (f[-1][:m.end()] + '.xyz'))
            #~ print("filefinder:", f, xyzf)
            if xyzf not in filetomatches:
                filetomatches[xyzf] = m
    
    groups = sorted([zip(names, map(Decimal, [m or 'nan' for m in mtch.groups()]))    
        + [('xyz',f)] for f,mtch in filetomatches.items()])
        
    pdicts = [dict(lst) for lst in groups]
    for p in pdicts:
        p.update(args)
    return [Namespace(d) for d in pdicts]

def groupdicts(dlst, key):
    bigdict = dict()
    for d in dlst:
        curval = d[key]
        keylist = bigdict.get(curval, [])
        keylist.append(d)
        bigdict[curval] = keylist
    return bigdict
