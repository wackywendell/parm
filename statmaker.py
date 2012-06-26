from __future__ import print_function
from xyzfile import XYZreader
import gzip, os, os.path
from numpy import array, mean
import math, numpy
from simw import Vec, autocorr, calc_Rg, geometric, average, average_squared
import simw as sim
import simpdb
import logging
from functools import wraps

mydir = os.path.expanduser('~/idp/')
pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'

class NotCalculableError(Exception): pass

class Statmaker:
    _statdict={}
    def __init__(self, cut=None, H=False, memoize=False,
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
            atomn = len(f0.locarray)
            #reslens = list(np.cumsum([len(r) for r in self._residues]))
            reslens = list(np.cumsum([len(r) for r in self._residues]))
            if atomn not in reslens:
                raise ValueError('Could not find %d in reslens: %s' % (atomn, reslens))
            indx = reslens.index(atomn)+1
            if indx <= len(self._residues):
                self._residues = self._residues[:indx]
                #~ print('Only %d residues' % len(self._residues))
                if hasattr(self, '_atoms'): del self._atoms
                #~ print(len(self.atoms), sum(len(r) for r in self.residues))
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
    
    def __call__(self, xyzfilename, *args, **kw):
        xyzfilename = str(xyzfilename)
        memokey = (xyzfilename, args, kw)
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
        
        
        with gzip.open(fname, 'w') as outf:
            try:
                return self._from_derived(self, xyzfilename, outf, *args, **kw)
            except NotCalculableError:
                failed.append('Derived' + repr(xyzfilename))
        
            with open(xyzfilename, 'r') as f:
                xyzreader = XYZreader(f)
                try:
                    vs = (v for v in xyzreader.vdicts())
                    return self._from_vdicts(vs, outf, *args, **kw)
                except NotCalculableError:
                    failed.append('From vdicts ' + repr(xyzfilename))
                    f.seek(0)
                except TypeError:
                    #temporary
                    print('args:',args)
                    print('kw:',kw)
                    raise
                
                try:
                    return self._from_frames(iter(xyzreader), outf=outf, *args, **kw)
                except NotCalculableError:
                    failed.append('From frames ' + repr(xyzfilename))
                    f.seek(0)
                
                try:
                    return self._from_atoms(
                            xyzreader.into(self.atoms, False),
                            outf=outf, *args, **kw)
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
    
    def _from_file(self, f, *args, **kw):
        raise NotCalculableError

########################################################################
########################################################################
########################################################################

class StatRunner(Statistic):
    """A Running statistics class."""
    def _from_file(self, f, *args, **kw):
        arr = self._cut(array(
                [map(float, l.strip().split('\t')) for l in f]
            ))
        shape = arr.shape
        if len(shape) == 2 and 1 in shape:
            return arr.flatten()
    
    def _from_vdicts(self, vs, outf, *args, **kw):
        vdictname = self.vdictname if hasattr(self, 'vdictname') else self.name
        if vdictname is None:
            raise NotCalculableError
        vals = array([float(v[vdictname]) for v in vs])
        if outf is not None:
            for v in vals:
                print(v, file=outf)
        return self._cut(vals)
        
class T(StatRunner): pass
class E(StatRunner): name='E'
class Times(StatRunner): vdictname='time'
class Rg(StatRunner): name='Rg'
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
        acorrs = sim.autocorr(Rgs)
        dt = mean(Times-Times[0])
        acorr = zip(numpy.arange(len(acorrs)), acorrs)
        try:
            tabove = (max((t for t,v in acorr if abs(v) > cutlast))
                    if cutlast < max(abs(acorrs)) else dt)
            t0 = min((t for t,v in acorr if abs(v) < cutfirst and t >= tabove))
        except ValueError:
            # It never gets in the acceptable region
            ts1 = [t for t,v in acorr if v > cutlast]
            #~ print(ts1[:3], ts1[-3:])
            if ts1:
                ts2 = [t for t,v in acorr if v < cutfirst and t >= max(ts1)]
                print("ERROR in statmaker.relax.run", max(ts1), ts2[:3], ts2[-3:])
            return None
        print(t0*dt, file=outf)
        return t0*dt

class Rij(StatRunner):
    def _from_vdicts(self, vs, outf, i,j):
        vdictname = '%s%d-%d' % (self.name,i,j)
        vals = array([float(v[vdictname]) for v in vs])
        for v in vals:
            print(v, file=outf)
        return self._cut(vals)

class Dihedral(StatRunner):
    vdictname=None
    def addtofname(self):
        return 'dihedral-%s' % self.ang
    
    def _single_from_residues(self, ang):
        res = self.residues
        restriplets = zip(res, res[1:], res[2:])
        return array([res.dihedral(ang, p, n) for p,res,n in restriplets])
    
    def _from_atoms(self, frames, outf, ang):
        vals = array([self._single_from_residues(ang) for f in frames], dtype=float)
        for v in vals:
            print('\t'.join(map(str,v)), file=outf)
        return self._cut(vals, cut)

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
    '-'.join([n + '((?:[.0-9]*)|(?:inf))' for n in names]) + r'.*')
    if 'regexp' in args: del args['regexp']
    regex = re.compile(regexpr)
    matches = [(regex.match(f[-1]), f) for f in xchildren]
    filetomatches = {}
    for m, f in list(matches):
        if m is None:
            if not matchall:
                matches.remove((m, f))
                continue
            raise ValueError('%s could not match %s' % (regexpr, f[-1]))
        elif None in m.groups():
            print('Groups:', m.groups())
            raise ValueError('%s could not entirely match %s' % (regexpr, f[-1]))
        elif m.start() != 0:
            print('Groups:', m.groups())
            raise ValueError('%s matched %s at the wrong place' % (regexpr, f[-1]))
        if f.extension == 'xyz':
            filetomatches[f] = m
        else:
            xyzf = fpath.File(f[:-1] + (f[-1][:m.end()] + '.xyz')
            if xyzf not in filetomatches:
                filetomatches[xyzf] = m
    
    groups = sorted([zip(names, map(Decimal, [m or 'nan' for m in mtch.groups()]))    
        + [('xyz',f)] for f,mtch in filetomatches.items()])
        
    pdicts = [dict(lst) for lst in groups]
    for p in pdicts:
        p.update(args)
    return [Namespace(d) for d in pdicts]
