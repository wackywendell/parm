# encoding: UTF8


import math, gzip, io, sys, os, os.path, logging, itertools
#from contextlib import nested
import numpy
#~ from scipy.spatial.distance import cdist
from decimal import Decimal, Context
def D(arg): return Context(prec=15).create_decimal(arg)
    
from namespace import Namespace
import simw as sim
from simw import autocorr
import xyzfile, simpdb

def gzopentxt(fname, *args, encoding='utf8'):
    """Open a gzip file for writing, and write str objects using encoding.
    
    Used for python 3.0, 3.1; not needed for python 3.2, but it works.
    Very hackish."""
    gzf = gzip.open(fname, *args)
    gzf.readable = lambda: False
    gzf.seekable = lambda: False
    gzf.writable = lambda: True
    if not hasattr(gzf, 'closed'):
        gzf.closed = False
    return io.TextIOWrapper(gzf, encoding=encoding)

def get_basename(f):
    basename, sep, ext = str(f).rpartition('.')
    if ext == 'xyz':
        return basename
    elif ext == 'gz':
        basename, sep, ext = str(basename).rpartition('.')
        return basename
    return basename

class Stat(object):
    AllStats = {}
    def __init__(self, calcfunc, headersfunc=None, name=None, add=True):
        self.name = calcfunc.__name__ if name is None else name
        self.calculate = calcfunc
        if headersfunc is not None: self.headers = headersfunc
        
        if not add: return
        if name in self.AllStats and self.AllStats[self.name] != self:
            raise ValueError("Stat cannot register " + repr(calcfunc) + 
                    ", key " + repr(name) + " already taken by " +
                    repr(self.AllStats[name]))
        self.AllStats[self.name] = self

    def calculate(self, collec, reslists):
        raise NotImplementedError
    
    def headers(self, reslists):
        return [self.name]
    
    def headerfunc(self, newfunc):
        self.headers = newfunc

def RStat(func):
    name = func.__name__
    def headers(reslists):
        return ([name] if len(reslists) <= 1 else
            ['%s_%d' % (name, n) for n in range(len(reslists))])
    
    s = Stat(func, headers)
    #@s.headerfunc
    #s.headers = headers
    return s

@Stat
def T(collec, reslists):
    return [collec.temp()]

@Stat
def E(collec, reslists):
    return [collec.energy()]

@Stat
def L(collec, reslists):
    return [collec.angmomentum().mag()]

@Stat
def Kinetic(collec, reslists):
    return [collec.kinetic()]

@Stat
def COMV(collec, reslists):
    return [collec.comv().mag()]

@RStat
def Rg(collec, reslists):
    return [sim.calc_Rg(reslist) for reslist in reslists]

@Stat
def Rgall(collec, reslists):
    allres = [r for rl in reslists for r in rl]
    return [sim.calc_Rg(allres)]

@Stat
def comd(collec, reslists):
    mgroups = [sim.metagroup(rl) for rl in reslists]
    colleccom = collec.com()
    coms = [g.com() for g in mgroups]
    ds = [(c1 - c2).mag() for n,c1 in enumerate(coms) 
                    for m,c2 in enumerate(coms[:n])]
    return [numpy.mean(ds)]

#~ def Rijs(ijs, name='Rijs', add=True):
    #~ def Rijs(collec, reslists):
        #~ return [(reslist[i-1]['CA'].x - reslist[j-1]['CA'].x).mag()
                    #~ for reslist in reslists for i,j in ijs]
    #~ def headers(reslists):
        #~ ijstrs = ['R%d_%d' % ij for ij in ijs]
        #~ if len(reslists) == 1: return ijstrs
        #~ else: return [ijstr + '_' + str(n) 
            #~ for n,reslist in enumerate(reslists) for ijstr in ijstrs]
    #~ return Stat(Rijs, headers, name=name, add=add)

#Rijs([(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
#            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)])


class TableWriter(object):
    def __init__(self, f, headers=None, gzip=True, sep='\t'):
        self.f = f
        self.gzip = gzip
        self.sep = sep
        self.headers = headers
        self.headers_written = False
        if headers is not None:
            assert all([isinstance(s, str) for s in headers]), headers
    
    def write_row(self, cols):
        txt = '\t'.join(cols)
        if not self.headers_written:
            self.headers_written = True
            mode = 'w'
            if self.headers is not None:
                txt = '\t'.join(self.headers) + '\n' + txt
        else:
            txt = '\n' + txt
            mode = 'a'
        
        if not isinstance(self.f, str):
            self.f.write(txt)
        else:
            with gzopentxt(self.f, mode) as f:
                f.write(txt)
            
class SimpleStat(object):
    def __init__(self, fname, collec, reslists, headers=None):
        self.fname = fname
        self.headers = headers
        self.collec = collec
        self.reslists = reslists
        self.table = []
        
    def update(self, t, write=False): raise NotImplementedError
    
    def write_table(self, headers=None, fmt='.6g'):
        with gzip.open(self.fname, 'w') as f:
            tw = TableWriter(f, headers=headers)
            for row in self.table:
                tw.write_row([format(n, fmt) for n in row])
    
    def write(self):
        self.write_table()

class Rijs(SimpleStat, TableWriter):
    ext = '.rij.gz'
    def __init__(self, fname, collec, reslists, rijlist, gzip=True,
                    contin=False):
        self.rijlist = rijlist
        headers = self.get_headers(reslists)
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
        if contin: self.headers_written=True
    
    def get_headers(self, reslists):
        if len(reslists) == 1:
            return ['time'] + ['R%d_%d' % (i,j)
                    for n,rl in enumerate(reslists)
                    for i,j in self.rijlist
                ]
        
        return ['time'] + ['CA%d_%d_%d' % (i,j,n+1) 
                    for n,rl in enumerate(reslists)
                    for i,j in self.rijlist
                ]
    
    def update(self, t, write=True):
        vals = [(reslist[i-1]['CA'].x - reslist[j-1]['CA'].x).mag()
                    for reslist in self.reslists for i,j in self.rijlist]
        
        statstrs = ['%.2f' % t] + [format(v, '.4f') for v in vals]
        
        assert len(statstrs) == len(self.headers), ("%d vs. %d" % 
                    (len(statstrs), len(self.headers)))
        self.write_row(statstrs)
    
    @staticmethod
    def read(fname):
        with gzip.open(fname, 'rb') as f:
            cols = f.readline().decode('ascii').strip().split('\t')
            table = numpy.genfromtxt(f, dtype=float)
        ts = table[:,0]
        ns = [k[1:].split('_') for k in cols[1:]]
        if all(len(n) == 2 for n in ns):
            istrs, jstrs = list(zip(*ns))
            ivals, jvals = list(map(int, istrs)), list(map(int, jstrs))
            rijd = dict([((i,j),table[:, cols.index(h)]) 
                    for i,j,h in zip(ivals, jvals, cols[1:])])
            return ts, rijd
        else:
            hls = ['%r (%d)' % (h,len(n)) for n,h in zip(ns, cols[1:])
                    if len(n) != 2]
            raise ValueError('Could not handle headers: ' + ('; '.join(hls)))
    
    def write(self):
        pass

class PairDists(SimpleStat):
    ext = '.psd.gz'
    def __init__(self, fname, collec, reslists, contin = False):
        nres = len(reslists[0])
        assert all(len(r) == nres for r in reslists)
        SimpleStat.__init__(self, fname, collec, reslists)
        if contin: self.load()
        else:
            self.n = 0
            self.table = numpy.zeros((nres*2, nres))
    
    def load(self, fname = None):
        with gzip.open(fname or self.fname, 'r') as f:
            self.n = int(f.readline().strip())
            self.table = numpy.genfromtxt(f)
        
    def update(self, t, write=False):
        nres = len(self.reslists[0])
        
        v1 = numpy.zeros((nres, 3), dtype=float)
        v2 = numpy.zeros((nres, 3), dtype=float)
        for rl1, rl2 in itertools.combinations(self.reslists, 2):
            for n, r in enumerate(rl1):
                v1[n] = r['CA'].x
            for n, r in enumerate(rl2):
                v2[n] = r['CA'].x
            # print("Shapes:", v1.shape, v2.shape)
            #dotprods = numpy.einsum('ij,kj',v1,v2)
            #~ dotprods = cdist(v1,v2)
            d0 = numpy.subtract.outer(v1[:,0], v2[:,0])
            d1 = numpy.subtract.outer(v1[:,1], v2[:,1])
            d2 = numpy.subtract.outer(v1[:,2], v2[:,2])
            dotprods = numpy.sqrt(numpy.sum((d0**2,d1**2,d2**2),0))

            #~ print("dotprod shape:", dotprods.shape)
            self.table[:nres] += dotprods
            self.table[nres:] += dotprods**2
            
        #~ for rl1, rl2 in itertools.combinations(self.reslists, 2):
            #~ for n1, r1 in enumerate(rl1):
                #~ for n2, r2 in enumerate(rl2):
                    #~ d = (r1['CA'].x - r2['CA'].x).mag()
                    #~ self.table[n1,n2] += d
                    #~ self.table[n1+nres,n2] += d**2
        self.n += 1
        if write: self.write()
    
    def write(self):
        self.write_table([str(self.n)], '.3f')
    
    @staticmethod
    def meansstds(fname):
        with gzip.open(fname, 'r') as f:
            n = int(f.readline().strip())
            table = numpy.genfromtxt(f)
        nrow, ncol = table.shape
        assert nrow == ncol*2
        sums, sumsq = table[:ncol], table[ncol:]
        means = sums / n
        stds = numpy.sqrt(sumsq / n - means**2)
        return means,stds

class DihedralHist(SimpleStat):
    ext = '.dih.gz'
    def __init__(self, fname, collec, reslists, contin = 0, sz=3000):
        SimpleStat.__init__(self, fname, collec, reslists)
        if contin: self.load(contin)
        else: 
            self.n = 0
            self.table = numpy.zeros((sz,sz))
            self.sz = sz
    
    def load(self, fname = None):
        with gzip.open(fname or self.fname, 'r') as f:
            self.n = int(f.readline().strip())
            self.table = numpy.genfromtxt(f)
        self.sz = len(self.table)
        
    def update(self, t, write=False):
        nres = len(self.reslists[0])
        for res in self.reslists:
            for p,res,n in zip(res, res[1:], res[2:]):
                phi, psi = res.dihedral('phi', p, n), res.dihedral('psi', p, n)
                phix = (phi + numpy.pi) * self.sz / (2 * numpy.pi)
                psix = (psi + numpy.pi) * self.sz / (2 * numpy.pi)
                if phix >= self.sz:
                    print("PHIX TOO BIG: %.6f (%.4g) -> %f" % (phi, phi-numpy.pi, phix))
                    phix = self.sz-1
                if psix >= self.sz:
                    print("PSIX TOO BIG: %.6f (%.4g) -> %f" % (psi, psi-numpy.pi, psix))
                    psix = self.sz-1
                if phix < 0:
                    print("PHIX TOO SMALL: %.6f -> %f" % (phi, phix))
                    phix = 0
                if psix < 0:
                    print("PSIX TOO SMALL: %.6f -> %f" % (psi, psix))
                    psix = 0
                self.table[int(phix),int(psix)] += 1
        self.n += 1
        if write: self.write()
    
    def write(self):
        self.write_table([str(self.n)], '.6f')
    
    @staticmethod
    def probs(fname):
        with gzip.open(fname, 'r') as f:
            n = f.readline().strip()
            table = numpy.genfromtxt(f)
        return table
    
class StatWriter(TableWriter):
    ext = '.tsv.gz'
    def __init__(self, fname, collec, reslists, stats=None,
                    gzip=True, contin=False):
        if stats == None:
            names, stats = list(zip(*sorted(Stat.AllStats.items())))
        headers = ['time'] + [h for s in stats for h in s.headers(reslists)]
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
        self.stats = stats
        if contin: self.headers_written=True
    
    def add_stat(self, s):
        self.headers += s.headers(self.reslists)
        self.stats += (s,)
        
    def update(self, t, write=True):
        def tostr(x):
            if isinstance(x, float): return '%.6g' % x
            return str(x)
        statstrs = [str(t)] + [tostr(x) 
                        for s in self.stats
                        for x in s.calculate(self.collec, self.reslists)]
        assert len(statstrs) == len(self.headers)
        self.write_row(statstrs)

class Dihedrals(SimpleStat, TableWriter):
    ext = '.dhl.gz'
    def __init__(self, fname, collec, reslists, gzip=True,
                    contin=False):
        headers = self.get_headers(reslists)
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
        if contin: self.headers_written=True
    
    @staticmethod
    def get_headers(reslists):
        if len(reslists) == 1:
            return ['time'] + ['%s%d' % (ang,m+1)
                    for ang in ('phi','psi')
                    for n,rl in enumerate(reslists)
                    for m,r in enumerate(rl[:-1])
                ]
                    
        return ['time'] + ['%s%d_%d' % (ang,m+1,n+1) 
                    for ang in ('phi','psi')
                    for n,rl in enumerate(reslists)
                    for m,r in enumerate(rl[:-1])
                ]
    
    @staticmethod
    def _get_dihedrals(rl):
        """Note that psis are for 0:N-1, and phis for 1:N. You need
        to do phis[:-1], psis[1:] to get matching pairs..."""
        phis,psis = simpdb.Resvec.phipsis(rl)
        return [v for angvals in (phis,psis) for v in angvals]
    
    def update(self, t, write=True):
        vals = list(zip(*[simpdb.Resvec.phipsis(rl) for rl in self.reslists]))
        statstrs = ['%.2f' % t] + [format(v, '.5f')
                    for angvals in vals
                    for rlvals in angvals
                    for v in rlvals]
        
        assert len(statstrs) == len(self.headers), (
                    "%d vs. %d" % (len(statstrs), len(self.headers)))
        self.write_row(statstrs)
    
    def write(self):
        pass
    
    @classmethod
    def readFile(cls, fname):
        """Note that psis are for 0:N-1, and phis for 1:N. You need
        to do phis[:-1], psis[1:] to get matching pairs..."""
        with gzip.open(fname, 'r') as f:
            cols = f.readline().strip().split('\t')
            table = numpy.genfromtxt(f, dtype=float)
        ts = table[:,0]
        if '_' not in cols[1]:
            idx = cols.index('psi1')
            phis,psis = table[:,1:idx], table[:,idx:]
            return ts, phis, psis
        
        startindxs = [n for n,hdr in enumerate(cols) 
            if hdr.split('_')[0] in ('phi1','psi1')]
        vecs = ([table[:,n:m]
            for n,m in zip(startindxs, startindxs[1:] + [None])])
        heads = ([cols[n:m]
            for n,m in zip(startindxs, startindxs[1:] + [None])])
        # now heads (vecs) is a list of headers (values) corresponding
        # to all phis or all psis in one chain
        
        phih, phiv = list(zip(*[(hs, vs) for hs,vs in zip(heads, vecs)
                    if all('phi' in h for h in hs)]))
        assert phih[0][0] == 'phi1_1', phih
        assert phih[0][1] == 'phi2_1', phih
        assert phih[0][10] == 'phi11_1', phih
        assert phih[1][0] == 'phi1_2', phih
        assert phih[1][1] == 'phi2_2', phih
        assert phih[1][10] == 'phi11_2', phih
        
        psih, psiv = list(zip(*[(hs, vs) for hs,vs in zip(heads, vecs)
                    if all('psi' in h for h in hs)]))
        assert psih[1][10] == 'psi11_2', psih
        
        assert len(phih) + len(psih) == len(heads)
        
        return ts, phiv, psiv

class CALocs(SimpleStat, TableWriter):
    ext = '.cax.gz'
    def __init__(self, fname, collec, reslists, gzip=True,
                    contin=False):
        headers = self.get_headers(reslists)
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
        if contin: self.headers_written=True
    
    @staticmethod
    def get_headers(reslists):
        if len(reslists) == 1:
            return ['time'] + ['CA%s%d' % (coord,m+1)
                    for n,rl in enumerate(reslists)
                    for m,r in enumerate(rl)
                    for coord in 'xyz'
                ]
        
        return ['time'] + ['CA%s%d_%d' % (coord,m+1,n+1) 
                    for n,rl in enumerate(reslists)
                    for m,r in enumerate(rl)
                    for coord in 'xyz'
                ]
    
    def update(self, t, write=True):
        com = self.collec.com()
        vals = [(r['CA'].x - com).get(coord)
                    for rl in self.reslists
                    for r in rl
                    for coord in (0,1,2)]
            
        statstrs = ['%.2f' % t] + [format(v, '.3f') for v in vals]
        
        assert len(statstrs) == len(self.headers), ("%d vs. %d" % 
                    (len(statstrs), len(self.headers)))
        self.write_row(statstrs)
    
    def write(self):
        pass
    
    @classmethod
    def readFile(cls, fname):
        with gzip.open(fname, 'r') as f:
            cols = f.readline().decode('ascii').strip().split('\t')
            table = numpy.genfromtxt(f, dtype=float)
        ts = table.T[0]
        if '_' not in cols[1]:
            xs,ys,zs = table[:,1::3], table[:,2::3], table[:,3::3]
            #del table
            return ts, numpy.array([xs,ys,zs]).T
        startindxs = [n for n,hdr in enumerate(cols) 
            if hdr.split('_')[0] == 'CAx1']
        vecs = ([numpy.array([table[:,n:m:3], table[:,n+1:m:3],
                                table[:,n+2:m:3]]).transpose((1,2,0))
            for n,m in zip(startindxs, startindxs[1:] + [None])])
        #del table
        return ts, vecs
    
    @classmethod
    def cdists(cls, fname):
        from scipy.spatial.distance import cdist
        ts, vecs = cls.readFile(fname)
        assert len(vecs) == 2, "Must have only two residue pairs to make cdicts"
        for t, v1, v2 in zip(ts, *vecs):
            yield t, cdist(v1, v2)
    
    @classmethod
    def Rgdist(cls, fname):
        ts, vecs = cls.readFile(fname)
        assert len(vecs) == 2, "Must have only two residue pairs to make cdicts"
        v1,v2 = vecs
        ds = numpy.mean(v1,1) - numpy.mean(v2,1)
        return ts, numpy.sqrt(numpy.mean(ds**2, 1))
        
    
class StatManager(object):
    def __init__(self, fname, collec, reslists, stats=None, 
            groups=None, gzip=True, contin=False):
        self.collec = collec
        self.reslists = reslists
        self.basename, sep, ext = fname.rpartition('.')
        assert ext != ''
        if ext == 'gz': self.basename, sep, ext = self.basename.rpartition('.')
        else: logging.warning('Unknown ext: ' + repr(ext))
        self.stats = stats
        if stats is None:
            self.stats = StatWriter(fname, collec,
                reslists, gzip=gzip, contin=contin)
        else:
            self.stats = StatWriter(fname, collec,
                reslists, stats=stats, gzip=gzip, contin=contin)
        self.groups = []
        for g in groups: self.add_group(g, contin=contin)
        
    def add_group(self, g, contin=False):
        if isinstance(g, SimpleStat): groups.append(g)
        elif isinstance(g, tuple):
            args = g[1:]
            g = g[0]
            self.groups.append(
                g(self.basename + g.ext, self.collec, self.reslists,
                *args, contin=contin))
        elif g is None:
            return
        else: self.groups.append(
                g(self.basename + g.ext, self.collec, self.reslists, contin=contin))
    
    def update(self, t, write=True):
        for s in [self.stats] + self.groups:
            s.update(t, write=write)
    
    def write(self):
        for s in self.groups: s.write()
        
    

########################################################################
def loadTable(fname, cut=0):
    with gzip.open(str(fname), mode='rb') as f:
        cols = [n.strip() for n in f.readline().decode('ascii').strip().split('\t')]
        arr = numpy.genfromtxt(f)
        #arr = numpy.genfromtxt(f, names=True)
    arrlen = len(arr)
    if 0.0 < cut < 1.0:
        cutlen = int(arrlen * cut)
        arr = arr[cutlen:]
    elif cut >= 1.0:
        raise NotImplementedError
    colnames = [(n if '-' not in n else ('R' + n.replace('-','_')))
                    for n in cols]
    d = Namespace(dict(zip(colnames, arr.T)))
    d.npts = arrlen
    return d
    #names = arr.dtype.names
    #return [(nm, arr[nm]) for nm in names]

#~ def load(d, cut=0):
    #~ d2 = loadTable(fname, cut=cut)
    #~ d.update(d2)
    #~ psdf = (f[:-1] + (f[-1][:-7] + '.psd.gz')).transform()
    #~ if 
        #~ d.psdmeans, d.psdstds = pairDists.meansstds(str(d.f)[:-7] + '.psd.gz', d.npts)
    
#~ def getRijs(d):
    #~ if isinstance(d, str):
        #~ d = loadTable(d)
    #~ ks = [k for k in dir(d) if k[0] == 'R' and '_' in k]
    #~ Rs = [k[1:].split('_') for k in ks]
    #~ Rs = [(int(i), int(j)) for i,j in Rs]
    #~ return dict([(ij, d[k]) for ij,k in zip(Rs,ks)])
            #~ 
#~ def getRij(d, ij):
    #~ return d['R%d_%d' % ij]
    #~ 
#~ def FRET(d,R0=54):
    #~ Rijs = getRijs(d)
    #~ if 'FRET' not in d: d.FRET = dict([(ij,numpy.mean(1.0/(1.0+(v/float(R0))**6)))
                        #~ for ij,v in getRijs(d).items()])
    #~ return d.FRET

def relax(Rgsets, Times, cutfirst = 1.0/math.e, cutlast=1.0):
    t=0
    for Rgs in Rgsets:
        acorrs = numpy.array(autocorr(Rgs))
        acorr = numpy.array(list(zip(numpy.arange(len(acorrs)), acorrs)))
        try:
            tabove = (max((t for t,v in acorr if abs(v) > cutlast))
                    if cutlast < max(abs(acorrs)) else 0)
            t0 = min((t for t,v in acorr if abs(v) < cutfirst and t >= tabove))
        except ValueError:
            # It never gets in the acceptable region
            return None
        t = max([Times[t0] - Times[0], t])
    return t

_Rgdist_ext = '.rgd.gz'
def Rgdist(d, read=True, write=True):
    bname = (str(d.dir + d.basename)
            if hasattr(d, 'basename') or 'basename' in d
            else get_basename(str(d)))
    fname = bname + _Rgdist_ext
    CAfname = bname + CALocs.ext
    if read:
        try:
            with gzip.open(fname, 'r') as f:
                cols = [n.strip() for n in f.readline().strip().split('\t')]
                arr = numpy.genfromtxt(f, delimiter='\t').T
            return arr[0], arr[1:]
        except Exception as e:
            #print("problem reading file", repr(fname))
            #raise e
            pass
    
    # couldn't get it from file
    ts, vecs = CALocs.readFile(str(CAfname))
    cols = ['time']
    nvecs = len(vecs)
    arr = numpy.zeros((1 + (nvecs * (nvecs-1))/2, len(ts)), float)
    arr[0] = ts
    for n, pair in enumerate(itertools.combinations(enumerate(vecs), 2)):
        #vs1, vs2 = pair
        #n1, v1 = vs1
        #n2, v2 = vs2
        (n1,v1), (n2, v2) = pair
        cols.append('%d_%d' % (n1+1, n2+1))
        ds = numpy.mean(v1,1) - numpy.mean(v2,1)
        arr[n+1] = numpy.sqrt(numpy.mean(ds**2, 1))
        
    assert len(cols) == len(arr)
    
    with gzip.open(fname, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        numpy.savetxt(f, arr.T, fmt='%.4f', delimiter='\t')
    return arr[0], arr[1:]


########################################################################
def filefinder(dir,  *names, **args):
    """
    extra kw:
    regexp          : Define the regular expression for groups
    matchall=False  : Match the entire expression
    types           : What types to apply (default: all Decimal)
    """
    import re, fpath
    matchall = args.get('matchall',True)
    if 'matchall' in args: del args['matchall']
    types = args.get('types', [Decimal] * len(names))
    if 'types' in args: del args['types']
    
    dir = fpath.Dir(dir)
    children = [f for f in dir.children() if f[-1][-7:] == '.tsv.gz']
    
    if 'regexp' in args:
        args['regexp'].replace(r'\F', r'((?:[0-9]*\.?[0-9]*)|(?:inf))')
    regexpr = (args['regexp'] if 'regexp' in args else
    '-'.join([n + r'((?:[0-9]*\.?[0-9]*)|(?:inf))' for n in names]))
    if 'regexp' in args: del args['regexp']
    regex = re.compile(regexpr)
    matches = [(list(regex.finditer(f[-1])), f) for f in children]
    badmatches = [f for ms,f in matches if len(ms) == 0]
    if badmatches and matchall:
        raise ValueError('%s could not match %s' % (regexpr, badmatches[0][-1]))
    matches = [(ms[0],f) for ms,f in matches if len(ms) > 0]
    
    
    
    groups = sorted([list(zip(names, [t(m or 'nan') for t,m in zip(types,mtch.groups())]))    
        + [('f',f)] for mtch,f in matches])
    
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

class StatGroup(Namespace):
    def __init__(self, ns, cut=0):
        dict.__init__(self, ns)
        import fpath
        d = loadTable(self.f, cut=cut)
        self.f = fpath.File(self.f)
        self.cut = cut
        for k in d:
            if k in self:
                self[k + '0'] = self[k]
        self.update(d)
        self.dir = fpath.Dir(self.f[:-1])
        self.basename, ext, null = self.f[-1].rpartition('.tsv.gz')
        assert null == ''
        assert ext == '.tsv.gz'
    
    def relax(self, cutfirst=1.0/math.e, cutlast=1.0/math.e):
        if not hasattr(self, '_relax'):
            self._relax = dict()
        cutkey = '%.8g-%.8g' % (cutfirst, cutlast)
        if cutkey in self._relax:
            return self._relax[cutkey]
        rgs, ts = ([self.Rg], self.time) if 'Rg' in self else (
        [self[k] for k in sorted(self) if k[:3] == 'Rg_'], self.time)
        rx = relax(rgs, ts, cutfirst, cutlast)
        self._relax[cutkey] = rx
        return rx
    
    def psd(self):
        if '_psd' not in self:
            fname = self.dir + (self.basename + '.psd.gz')
            self['_psd'] = PairDists.meansstds(str(fname))
        return self._psd

    def dih(self, npts=1):
        if '_dih' not in self:
            fname = self.dir + (self.basename + DihedralHist.ext)
            self['_dih'] = DihedralHist.probs(str(fname))
        if npts == 1:
            return self._dih
        m,n = self._dih.shape
        #print(m,n,m/npts,n/npts,m*n,(m/npts)*(n/npts)*npts*npts)
        return self._dih.reshape((m/npts,npts,n/npts,npts)).sum(3).sum(1)
    
    def dhl(self):
        """Note that psis are for 0:N-1, and phis for 1:N. You need
        to do phis[:-1], psis[1:] to get matching pairs..."""
        fname = self.dir + (self.basename + Dihedrals.ext)
        return Dihedrals.readFile(str(fname))
    
    def cax(self):
        fname = self.dir + (self.basename + CALocs.ext)
        ts, cax = CALocs.readFile(str(fname))
        indx = ts.searchsorted(self.time[0])
        return ts[indx:], cax[:, indx:, :]
        
    def CAstats(self):
        from scipy.io import savemat, loadmat
        statfname = self.dir + (self.basename + '.castats.mat')
        fname = self.dir + (self.basename + CALocs.ext)
        if (os.path.exists(str(statfname))): 
            if (not os.path.exists(str(fname)) or 
            os.path.getmtime(str(statfname)) >= os.path.getmtime(str(fname))):
                d = loadmat(str(statfname))
                return d['caxs'], d['bonds'], d['angs'], d['dihs']
            else:
                logging.warning("Ignoring " + str(statfname))
        
        ts, caxs = self.cax()
        caxswapped = numpy.swapaxes(caxs,0,1)
        array = numpy.array
        Vec = sim.Vec
        bonds = array([array([(Vec(*r1) - Vec(*r2)).mag()
                for r1,r2 in zip(caxbyt, caxbyt[1:])])
                for caxbyt in caxswapped])
        angs = array([array([Vec.angle(Vec(*r1), Vec(*r2), Vec(*r3))
                for r1,r2,r3 in zip(caxbyt, caxbyt[1:], caxbyt[2:])])
                for caxbyt in caxswapped])
        dihs = array([array([Vec.dihedral(Vec(*r1), Vec(*r2), Vec(*r3), Vec(*r4))
                for r1,r2,r3,r4 in zip(caxbyt, caxbyt[1:], caxbyt[2:], caxbyt[3:])])
                for caxbyt in caxswapped])
        savemat(str(statfname), {'caxs':caxs, 'bonds':bonds,'angs':angs, 'dihs':dihs})
        return caxs, bonds, angs, dihs
    
    def cdists(self):
        fname = self.dir + (self.basename + CALocs.ext)
        return CALocs.cdists(str(fname))
    
    def Rgdist(self):
        fname = self.dir + (self.basename + CALocs.ext)
        return CALocs.Rgdist(str(fname))
    
    def Rijs(self):
        #print('CALCULATING RIJS', self.basename)
        #raise ValueError
        fname = self.dir + (self.basename + Rijs.ext)
        ts, d =  Rijs.read(str(fname))
        indx = ts.searchsorted(self.time[0])
        
        def cutout(arr):
            return arr[indx:]
        for ij in sorted(d):
            d[ij] = cutout(d[ij])
        return cutout(ts), d
        
    def FRET(self, R0=54):
        return self.FRETerr()[0]
    
    def FRETerr(self, R0=54):
        if '_FRETerr' in self: 
            return self['_FRETerr']
        
        cutstr = ('cut: %.6f' % self.cut).strip()
        fname = self.dir + (self.basename + '.etf.gz')
        rijfname = self.dir + (self.basename + Rijs.ext)
        if (os.path.exists(str(fname)) and 
            os.path.exists(str(rijfname)) and
            os.path.getmtime(str(rijfname)) <= os.path.getmtime(str(fname))):
            with gzip.open(str(fname), 'r') as f:
                cutline = f.readline().strip()
                if cutline == cutstr:
                    myijs = [tuple(map(int, c.split('_'))) for c in 
                                f.readline().strip().split('\t')]
                    myETeffs, myETerrs = numpy.genfromtxt(f, dtype=float)
            
            #print('(returned) ', end='')
            if cutline == cutstr:
                self['_FRETerr'] = (dict(zip(myijs, myETeffs)), dict(zip(myijs, myETerrs)))
                return self['_FRETerr']
        
        cuts = 1.0/math.e, 1.0/math.e
        rx = self.relax(*cuts)
        N = (self.time[-1] - self.time[0]) / rx
        
        ts, rijs = self.Rijs()
        myijs = sorted(rijs.keys())
        ETeff_momentary = [1.0/(1.0+(rijs[ij]/float(R0))**6) for ij in myijs]
        myETeffs = [numpy.mean(ET) for ET in ETeff_momentary]
        myETerrs = [numpy.std(ET) / numpy.sqrt(N) for ET in ETeff_momentary]
        
        self['_FRETerr'] = (dict(zip(myijs, myETeffs)), dict(zip(myijs, myETerrs)))
        
        cols = ['%d_%d' % ij for ij in myijs]
        effstrs = ['%.8f' % v for v in myETeffs]
        errstrs = ['%.8f' % v for v in myETerrs]
        with gzopentxt(str(fname), 'wb') as f:
            print(cutstr, file=f)
            print(*cols, sep='\t', file=f)
            print(*effstrs, sep='\t', file=f)
            print(*errstrs, sep='\t', file=f)
        
        return self['_FRETerr']
    
    def _old_FRETerr(self, R0=54):
        cuts = 1.0/math.e, 1.0/math.e
        rx = self.relax(*cuts)
        N = (self.time[-1] - self.time[0]) / rx
        if '_FRETerr' not in self: 
            Rijs = self.Rijs()
            self['_FRETerr'] = dict(
                [(ij,numpy.std(1.0/(1.0+(v/float(R0))**6)) / numpy.sqrt(N))
                                for ij,v in list(Rijs.items())])
        return self['_FRETerr']

def _default_str(d):
    if 'Rg' in d:
        t0, tf = d.time[0]/1000, d.time[-1]/1000
        t = tf - t0
        rlx = d.relax(1.0/math.e, 1.0/math.e)/1000.0
        return '%8.2f - %8.2f (%6.1f t₀, %6.1f x), Rg = %6.2f +- %5.2f' % (
            t0, tf, rlx, tf/rlx, numpy.mean(d.Rg), numpy.std(d.Rg))
    #return '%8.2f - %8.2f (%6.1f)' % (d.time[0]/1000, d.time[-1]/1000, relax(d, .1, .1)/1000.0)
    rlx = d.relax(1.0/math.e, 1.0/math.e)/1000.0
    return '%8.2f - %8.2f (%6.1f)' % (d.time[0]/1000, d.time[-1]/1000, rlx)
    
def loadAll(fs, cut=0, strfunc=_default_str, func=None):
    fstrs = [f.f[-1][:-7] for f in fs]
    maxlen = max(len(s) for s in fstrs)
    def load(f):
        print(format(f.f[-1][:-7], str(maxlen) + 's'), end=': ')
        sys.stdout.flush()
        d = StatGroup(f, cut=cut)
        if func is not None: func(d)
        print(strfunc(d))
        return d
    return [load(f) for f in fs]

########################################################################
mydir = os.path.expanduser('~/idp/')
loadfile= mydir + 'blengths/stats.pkl'
aSfile = mydir + 'pdb/aS.pdb'
taufile = mydir + 'pdb/tau_phyre2.pdb'
aSijs = [(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)]
aSijs += [(1, 140), (12, 127), (14, 125), (21, 118), 
            (21, 111), (36, 75), (38, 77), (43, 82), (57, 133),
            (59, 135), (63, 139), (9, 130), (33, 130), (54, 130),
            (72, 130), (92, 130), (33, 72), (9, 54), (72, 92),
            (54, 72), (9, 72), (9, 33), (54, 92)]
tauijs = [ (291, 322), (291, 354), (354, 433), (103, 184),
        ( 17, 103), (184, 291), (244, 354), (322, 433), 
        (291, 433), (103, 291), ( 17, 291), ( 17, 433)]

class Updater(object):
    _n_to_res = {1013: ([(aSfile, aSijs, False)]), 
                 2026: ([(aSfile, aSijs, False), (aSfile, aSijs, False)]), 
                 2016: ([(aSfile, aSijs, True)]),
                 3217: ([(taufile, tauijs, False)]),
                 } 
    
    def __init__(self):
        self._residues = {}
        self._collecs = {}
    
    def residues(self, n):
        if n in self._residues: return self._residues[n]
        if n not in self._n_to_res:
            raise KeyError("Could not recognize %d atoms" % n)
        reslists = self._residues[n] = [simpdb.Resvec.from_pdb('sample', 
                    pdbfile, loadfile, H=H, numchains=1)
                    for pdbfile, ijset, H in self._n_to_res[n]]
        return reslists
    
    def collec(self, n):
        if n in self._collecs: return self._collecs[n]
        avecs = sim.avector([r for rlist in self.residues(n) for r in rlist])
        return sim.StaticCollec(avecs)
    
    def run(self, xyzfname, stattypes, stats=None):
        print('%16s : %s' % (xyzfname, "Reading xyz..."))
        with open(xyzfname, 'r') as f:
            reader = xyzfile.XYZreader(f)
            for frame in reader:
                if 'stage' in frame.values: continue
                firstframe = frame
            frames = reader
                     #xyzfile.Frames([frame for frame in reader 
                     #       if 'stage' not in frame.values])
        
        print('%16s : %s' % (xyzfname, "Setting up residues..."))
        n = len(firstframe.locarray)
        pdbfile, ijset, H = self._n_to_res[n][0]
        
        collec = self.collec(n)
        reslists = self.residues(n)
        basename, ext, null = xyzfname.rpartition('.xyz')
        assert ext == '.xyz'
        assert null == ''
        
        atoms = [a for rlist in reslists for r in rlist for a in r]
        
        statfs = [basename + stat.ext for stat in stattypes]
        
        if stats is None:
            stats = StatWriter(basename + StatWriter.ext, collec, reslists)
        
        statgs = [(stat(f, collec, reslists) if stat is not Rijs
                                else Rijs(f, collec, reslists, ijset))
                for stat,f in zip(stattypes, statfs)]
        with open(xyzfname, 'r') as f:
            for t, vdict in xyzfile.XYZreader(f).into(atoms):
                if 'stage' in vdict: continue
                if stats is not False: stats.update(t)
                for stat in statgs:
                    stat.update(t, write=False)
                if (int(t) % 10000) == 0: print('%16s : %8d' % (xyzfname, t/1000))
        
        print('%16s : %s' % (xyzfname, "Writing..."))
        for stat in statgs:
            stat.write()
        print('%16s : %s' % (xyzfname, "Done."))

updater = Updater()

def update(xyzfname, stats=False, *args):
    basename, ext, null = xyzfname.rpartition('.xyz')
    assert ext == '.xyz'
    assert null == ''
    statgs = []
    for a in args:
        if not type(a) == type:
            raise TypeError("%s is not a class" % a)
        if not issubclass(a, SimpleStat):
            raise TypeError("%s is not a subclass of SimpleStat" % a)
        statgs.append(a)
    updater.run(xyzfname, statgs, stats=stats)

if __name__ == '__main__':
    import sys, argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r',dest='get_rijs', action='store_true')
    parser.add_argument('-d',dest='get_dihedrals', action='store_true')
    parser.add_argument('-p',dest='get_pairs', action='store_true')
    parser.add_argument('-c',dest='get_CA_locs', action='store_true')
    parser.add_argument('-s',dest='get_stats', action='store_true')
    parser.add_argument('-R', dest='Rgdist', action='store_true')
    parser.add_argument('files', nargs='+', type=str)
    opts = parser.parse_args()
    
    statgs = []
    if opts.get_dihedrals: statgs.append(Dihedrals)
    if opts.get_pairs: statgs.append(PairDists)
    if opts.get_CA_locs: statgs.append(CALocs)
    if opts.get_rijs: statgs.append(Rijs)
    for f in opts.files:
        if statgs or opts.get_stats:
            update(f, (None if opts.get_stats else False), *statgs)
        if opts.Rgdist:
            Rgdist(f, read=False)
