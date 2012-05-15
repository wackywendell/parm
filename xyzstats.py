from __future__ import print_function

import sys, os.path, shelve, gdbm, itertools, collections, math, cmath
import numpy as np
from Bio.PDB import PDBParser

import simpdb
from simw import Vec, autocorr, calc_Rg, geometric, average, average_squared
import simw as sim
from xyzfile import XYZreader, Frame, Frames
from functools import wraps
#~ import numpy as np

import logging

bondspring = 1000
anglespring = 1000
LJepsilon = 1
mydir = os.path.expanduser('~/idp/')
pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'

def _use_key(key):
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return self._shelf_out(key) or self._shelf_in(key, func(self, *args, **kwargs))
        return wrapper
    return decorator

def _key_and_cut(key):
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            vals = self._shelf_out(key)
            if vals is None: vals = self._shelf_in(key, func(self, *args, **kwargs))
            return self._cut(vals)
        return wrapper
    return decorator

class statkeeper:
    def __init__(self, xyzfname, shelfname='auto', pdbfname = pdbfile, 
                    loadfname = loadfile, H=False, cut=0, ignoreblank=True, open=True):
        if xyzfname is not None: xyzfname = unicode(xyzfname)
        self.xyzfname = xyzfname
        self.fname = xyzfname if xyzfname else shelfname
        if shelfname == 'auto':
            base, sep, _ = xyzfname.rpartition('.')
            self.shelfname = base + sep + 'stats'
        elif shelfname is not None:
            self.shelfname = unicode(shelfname)
        self.cut = cut
        # build on demand
        self._resfunc = lambda: simpdb.Resvec.from_pdb('aS', pdbfname, loadfname, H=H, numchains=1)
        if open: self.open()
        else:
            self.reader = None
            self.shelf = None
        
    def open(self):
        self.reader = XYZreader(open(self.xyzfname,'r')) if self.xyzfname else None
        
        if self.shelfname:
            try:
                self.shelf = shelve.Shelf(gdbm.open(self.shelfname, 'c'))
            except gdbm.error:
                print('Could not open', self.shelfname, file=sys.stderr)
                raise
        else:
            self.shelf = None
            
            #~ mkey = '_md5'
        if self.shelfname and self.reader: # check that they match
            mkey = '_size'
            rsize = self.reader.size()
            #~ print("Size:", self.reader.size(), 'Ignoring:', 
                        #~ self.reader.size() <= 1 and ignoreblank)
            if rsize <= 1 and ignoreblank:
                self.reader = None
                return
            
            elif rsize <= 10000:
                raise NotImplementedError('Fix me')
            
            elif mkey in self.shelf and self.shelf[mkey] == rsize: #.md5():
                logging.info('shelf exists, size matches')
                # So... ISF and autocorr depend on cut. Ugh.
                #~ if ((cut and not 'cut' in self.shelf) or 
                    #~ ('cut' in self.shelf and cut != self.shelf['cut'])):
                    #~ badks = [k for k in self.shelf if 
                            #~ ('ISF' in k or 'autocorr' in k or k == 'cut')]
                    #~ # they use a different cut
                    #~ for k in badks:
                        #~ del self.shelf[k]
                    
            else:
                if '_md5' not in self.shelf and mkey not in self.shelf:
                    logging.info('New shelf')
                else: logging.info('Size does not match')
                self.shelf.clear()
                self.shelf[mkey] = self.reader.size() #.md5()
                self.shelf.sync()
    
    def _cut(self, lst):
        """Cuts the first portion off a list, based on self.cut."""
        if 0 < self.cut < 1:
            numcut = int(len(lst)* self.cut)
            #~ print(numcut,'/',len(lst), 'frames cut', lst[0], lst[numcut], lst[-1])
            return lst[numcut:]
        elif self.cut >= 1:
            numcut = len(lst)- len(self.times)
            #~ print(numcut, 'frames cut.')
            return lst[numcut:]
        return lst
    
    def sync(self):
        if hasattr(self, 'shelf') and self.shelf: 
            self.shelf.sync()

    def close(self):
        if self.shelf is not None: 
            self.shelf.close()
            self.shelf = None
        if self.reader is not None:
            self.reader.close()
            self.reader = None
        for attr in ['_residues','_atoms','_times','_frames', '_vdicts','_collec']:
            if hasattr(self, attr):
                delattr(self, attr)
    
    def __del__(self):
        if hasattr(self, 'shelf') and self.shelf: 
            self.shelf.close()
        if hasattr(self, 'reader') and self.reader:
            self.reader.close()
    
    def __enter__(self):
        if self.reader is None and self.shelf is None:
            self.open()
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()
    
    @property
    def residues(self):
        if not hasattr(self, '_residues'):
            logging.info('building residues')
            self._residues = self._resfunc()
            f0 = self.frames()[0]
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
    
    def frames(self, usecut=True):
        if hasattr(self, '_frames'):
            return self._frames
        logging.info('importing frames')
        out = self._frames = self.reader.all()
        if usecut and 0 < self.cut < 1:
            numcut = int(len(self._frames)* self.cut)
            out = Frames(self._frames[numcut:])
            #~ print(numcut, 'frames cut')
        elif usecut and self.cut >= 1:
            oldlen = len(self._frames)
            out = Frames([f for f in self._frames if f.time > self.cut])
            #~ print(oldlen - len(self._frames), 'frames cut.')
        #~ else:
            #~ print('cut:', self.cut)
        logging.info('frames imported')
        return out
    
    @property
    def vdicts(self):
        key = '_vdicts'
        if not hasattr(self, '_vdicts'):
            logging.info('importing vdicts')
            self._vdicts = (self._shelf_out(key) or
                        self._shelf_in(key, self.reader.vdicts()))
                
            logging.info('vdicts imported')
        return self._vdicts
    
    @property
    def times(self):
        if not hasattr(self, '_times') or self._times is None:
            out = self._shelf_out('times')
            self._times = np.array(out if out is not None else [int(d['time']) for d in self.vdicts])
            self._shelf_in('times', self._times)
        if self.cut < 1:
            return np.array(self._cut(self._times))
        else:
            return np.array([t for t in self._times if t > self.cut])
    
    @property
    def collec(self):
        if not hasattr(self, '_collec'):
            self._collec = sim.StaticCollec(self.residues)
        return self._collec
    
    def add_bonds(self, springk):
        self._bonds = simpdb.make_bonds(self.residues, springk)
        self.collec.addInteraction(self._bonds)
        
    def add_angles(self, springk):
        self._angles = simpdb.make_angles(self.residues, springk)
        self.collec.addInteraction(self._angles)
        
    def add_LJ(self, eps, cutoff=2.0):
        self._LJ, self._neighbors = simpdb.make_LJ(self.residues, LJepsilon, 2.0)
        self.collec.addInteraction(self._LJ)
        self.collec.addTracker(self._neighbors)
    
    def add_dihedral(self, springk):
        self._dihedral = simpdb.make_dihedrals(self.residues, springk)
        self.collec.addInteraction(self._dihedral)
    
    def _shelf_out(self, key):
        o=object()
        if self.shelf is not None:
            val = self.shelf.get(key, o)
            if val is not o:
                logging.debug('found %s', key)
                return self.shelf[key]
        logging.info('key %s not found in shelf', key)
        return None
    
    def _shelf_in(self, key, val):
        if self.shelf is not None:
            self.shelf[key] = val
            logging.debug('putting ' + key)
        else:
            logging.debug('shelf does not exist %r' % self.shelf)
        return val
    
    @_key_and_cut('endtoend')
    def endtoend(self):
        a1, a2 = self.residues[0]['CA'], self.residues[-1]['CA']
        return [(a1.x - a2.x).mag() for t in self.frames(False).into(self.atoms)]
            
    @_key_and_cut('gyradius')
    def gyradius(self):
        return [self.collec.gyradius() for t in self.frames(False).into(self.atoms)]
            
    @staticmethod
    def getVal(key, frame, atoms, func, *args, **kw):
        if hasattr(frame, 'vals') and key in frame.vals:
            return frame.vals[key]
        frame.into(atoms)
        return func(*args, **kw)
    
    def getVals(self, vkey, func, conv=lambda x:x, *args, **kw):
        o = object()
        values = [d.get(vkey, o) for d in self.vdicts]
        values = [(v if v is o else conv(v)) for v in values]
        if o not in values:
            #~ print('Found', vkey, '::', values[:20], '...')
            return np.array(values)
        
        #~ print("getVals Couldn't always find", vkey)
        def getsingle(v, f):
            if v is not o:
                return v
            f.into(self.atoms)
            return func(*args, **kw)
        
        return [getsingle(v,f) for v,f in zip(values, self.frames(False))]
    
    @_key_and_cut('Rg')
    def Rg(self):
        def getRg():
            r = calc_Rg(self.residues)
            #~ print('Rg:', r)
        vals = self.getVals('Rg', getRg, float)
        return np.array(vals, dtype=float)
        
        #~ return self.getVals('Rg', lambda: calc_Rg(self.residues), float)
        #~ val = self._shelf_out('Rg')
        #~ if val: return val
        #~ return [calc_Rg(self.residues) for t in self.frames.into(self.atoms)]
            
            #~ locs = [r['CA'].x for r in self.residues]
            #~ center = sum(locs, Vec(0,0,0)) / len(locs)
            #~ locsqs = [(loc - center).sq() for loc in locs]
            #~ Rgs.append(math.sqrt(average(locsqs)))
        #~ return self._shelf_in('Rg', Rgs)
    
    @_key_and_cut('T')
    def temp(self):
        return self.getVals('T', lambda: self.collec.temp(), float)
        #~ return [float(self.getVal('T', f, self.atoms, self.collec.temp))
                        #~ for f in self.frames]
    
    @_key_and_cut('E')
    def energy(self):
        return self.getVals('E', lambda: self.collec.energy(), float)
        #~ return [float(self.getVal('E', f, self.atoms, self.collec.energy))
                        #~ for f in self.frames]
    
    def Rij(self, i, j):
        key = 'Rij%d-%d' % (i,j)
        func = lambda: self.getVals(key, lambda: sim.Rij(self.residues, i, j), float)
        vals = self._shelf_out(key)
        if vals is None: vals = self._shelf_in(key, func())
        return self._cut(vals)
    
    def _calcdihedrals(self, ang):
        restriplets = zip(self.residues, self.residues[1:], self.residues[2:])
        return np.array([res.dihedral(ang, p, n) for p,res,n in restriplets])
    
    @_key_and_cut('allPhis')
    def getPhis(self):
        return np.array([self._calcdihedrals('phi') 
                        for t in self.frames(False).into(self.atoms)])
        
    @_key_and_cut('allPsis')
    def getPsis(self):
        return np.array([self._calcdihedrals('psi')
                        for t in self.frames(False).into(self.atoms)])
        
    @_key_and_cut('allOmegas')
    def getOmegas(self):
        return np.array([self._calcdihedrals('omega') 
                        for t in self.frames(False).into(self.atoms)])
        
                
        
    def getDihedral(self, ang, resnum):
        key = 'dihedral-%s-%d' % (ang,resnum)
        def getfunc():
            def calcfunc():
                res1, res2, res3 = self.residues[resnum-1:resnum+2]
                return res2.dihedral(ang, res1, res3)
            return self.getVals(key, calcfunc, float)
        vals = self._shelf_out(key) or self._shelf_in(key, getfunc())
        return self._cut(vals)
    
    def getOmega(self, resnum):
        return self.getDihedral('omega', resnum)
    
    def getPhi(self, resnum):
        return self.getDihedral('phi', resnum)
    
    def getPsi(self, resnum):
        return self.getDihedral('psi', resnum)
    
    def ETeff(self, i, j, R0):
        Rijs = np.array(self.Rij(i, j), dtype=float)
        return np.mean(1/(1 + (Rijs / R0)**6))
        
    #~ @_key_and_cut('Rg')
    #~ def Rg(self):
        #~ return [float(self.getVal('Rg', f, self.atoms, lambda: calc_Rg(self.residues)))
                        #~ for f in self.frames]
    
    def ISF(self, scale, maxavg=200, ntimes=200):
        """Calculate the intermediate scattering function.
        Maxavg=n means that for a given dt, windows will be spaced evenly
        such that between n and 2n windows are used.
        ntimes refers to the number of points taken; they will be approximately
        a geometric series.
        """
        key = 'ISF-%.3g' % scale
        val = self._shelf_out(key)
        if not val: val = dict()
        if maxavg is None: maxavg = len(self.frames)
        if ntimes is None: ntimes = len(self.frames)
        valkey = (maxavg, ntimes)
        for mavg, nt in val:
            if mavg >= maxavg and nt >= ntimes:
                return val[(mavg,nt)]
        if val:
            logging.info('autocorr ignoring:', '(%d,%d)'%valkey, val.keys())

        if not self.times:
            return ([],[])
        
        locs = [f.locarray for f in self.frames]
        times = [t - self.times[0] for t in self.times]
        
        ns, isfpoints = sim.ISF(locs, scale, maxavg, ntimes)
        times = [times[n] for n in ns]
        
        val[valkey] = (times, isfpoints)
        self._shelf_in(key, val)
        return val[valkey]
    
    def ISF_Rg(self, maxavg=200, ntimes=200):
        Rgs = self.Rg()
        return self.ISF(np.mean(Rgs), maxavg, ntimes)
        
    def autocorr(self):
        #~ key = 'autocorr-Rg'
        #~ val = self._shelf_out(key)
        #~ if val: return val
        
        Rgs = np.array(self.Rg(), dtype=float)
        if self.times is None or len(self.times) == 0:
            return ([],[])
        t0 = self.times[0]
        dts = [t - t0 for t in self.times]
        return (np.array(dts, dtype=float), np.array(autocorr(Rgs), dtype=float))
        #~ return self._shelf_in(key, (dts, autocorr(Rgs)))
    
    def relax_acorr(self, cutfirst = 1.0/math.e, cutlast=1.0):
        """Relaxation time, calculated using the autocorrelation function 
        of the Rg timeseries.
        
        The relaxation time is calculated as the first time at which the
        autocorrelation function is below cutfirst and stays below cutlast.
        """
        dts, vals = acorr = self.autocorr()
        #~ print('acorr:', acorr)
        try:
            tabove = (max((t for t,v in zip(*acorr) if abs(v) > cutlast))
                    if cutlast < max(abs(vals)) else dts[0])
            t0 = min((t for t,v in zip(*acorr) if abs(v) < cutfirst and t >= tabove))
        except ValueError:
            # It never gets in the acceptable region
            ts1 = [t for t,v in zip(*acorr) if v > cutlast]
            #~ print(ts1[:3], ts1[-3:])
            if ts1:
                ts2 = [t for t,v in zip(*acorr) if v < cutfirst and t >= max(ts1)]
                print(max(ts1), ts2[:3], ts2[-3:])
            return None
        return t0
    
    def N_acorr(self, cutfirst = 1.0/math.e, cutlast=1.0):
        key = 'N-%.4f-%.4f' % (cutfirst, cutlast)
        val = self._shelf_out(key)
        if val: return val
        relax = self.relax_acorr(cutfirst, cutlast)
        if relax is None:
            return self._shelf_in(key, None)
        return self._shelf_in(key, float(self.timespan()) / float(relax))
    
    def relax_ISF(self, scale=None, maxavg=200, ntimes=200):
        """Relaxation time, calculated using the intermediate scattering
         function. If scale is None, average Rg is used."""
        isf = (self.ISF(scale, maxavg, ntimes) if scale 
                                else self.ISF_Rg(maxavg, ntimes))
        ts = [(t,v) for t,v in zip(*isf) if v < 1.0/math.e]
        if not ts:
            return None
        t0,v0 = min(ts)
        return t0
    
    
    def timespan(self):
        return self.times[-1] - self.times[0]
    
    def sim_Rg(self, relaxt):
        """Returns (average Rg, std, error, Npoints). Error is determined by 
        std(Rgs) / sqrt(Npoints), where Npoints is the length of the 
        simulation over the relaxation time.
        
        If 1 relaxation time was not reached, error is returned as -std(Rgs).
        """
        Rgs = self.Rg()
        Rg = average(Rgs)
        #~ relaxt = self.relax_Rg()
        std = float(np.std(Rgs))
        
        if relaxt is None or relaxt <= 0 or relaxt >= self.timespan():
            #~ logging.debug('Rg_avg: %.3f (%.3f); std: %.3f',Rg,average(self.gyradius()),std)
            return Rg, std, -std, None
        
        Npoints = self.timespan() / relaxt
        #~ logging.debug('Rg_avg: %.3f (%.3f); std: %.3f; N: %.3f',
                            #~ Rg,average(self.gyradius()),std,Npoints)
        #~ logging.debug('Rg_avg: %.3f (%.3f); std: %.3f; N: %.3f; err: %.3f',
                #~ Rg,average(self.gyradius()),std,Npoints,std/math.sqrt(Npoints-1))
        return Rg, std, std/math.sqrt(Npoints-1), Npoints
    
    @_key_and_cut('bondstds')
    def bond_stds(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        bonds = simpdb.make_bonds(self.residues, 0)
        return [bonds.std_dists() for t in self.frames.into(self.atoms)]
    
    @_key_and_cut('bondmeans')
    def bond_means(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        bonds = simpdb.make_bonds(self.residues, 0)
        return [bonds.mean_dists() for t in self.frames.into(self.atoms)]
    
    def bond_std(self):
        """Return bond distance std. dev. averaged over all atoms and frames"""
        return average_squared(self.bond_stds())
    
    @_key_and_cut('anglestds')
    def angle_stds(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        angles = simpdb.make_angles(self.residues, 0)
        return [angles.std_dists() for t in self.frames.into(self.atoms)]
    
    @_key_and_cut('anglemeans')
    def angle_means(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        angles = simpdb.make_angles(self.residues, 0)
        return [angles.mean_dists() for t in self.frames.into(self.atoms)]
    
    def angle_std(self):
        """Return bond distance std. dev. averaged over all atoms and frames"""
        return average_squared(self.angle_stds())
    
    @classmethod
    def runfiles(cls, fnames, func, *args, **kwargs):
        for fname in fnames:
            sk = cls(fname)
            func(sk, fname, *args, **kwargs)
    
    def old_ISF(self, lengthscale, ntimes=None, maxcount=100):
        key = 'autocorr-%.3g' % lengthscale
        val = self._shelf_out(key)
        #~ print('autocorr key',key, bool(val))
        if ntimes is None:
            ntimes = self.frames[-1].indx - self.frames[0].indx
        if val and (len(val) >= ntimes or len(val) >= len(self.times)-1):
            return val
        elif val:
            logging.info('autocorr ignoring:', len(val), '<', ntimes, '(%d)' % len(self.times))
            print('val existed.')
            print(bool(val), type(val), len(val), ntimes, len(val) >= ntimes, len(val) >= len(self.times)-1)
        #~ print('Running autocorr.')
        nframes = self.frames[-1].indx - self.frames[0].indx
        nindices = geometric(nframes, int(ntimes))
        ntimes = [self.frames[i].time - self.frames[0].time for i in nindices]
        corrs = [(t,self._autocorr_n(i, lengthscale, maxcount)) for i,t in zip(nindices, ntimes)]
        
        return self._shelf_in(key, dict(corrs))
        
        exponents = collections.defaultdict(list)
        #~ print('a:',v2)
        infostep = len(self.frames) / 10
        infoloc = infostep
        for (frame1, frame2) in itertools.combinations(self.frames,2):
            dindx = frame2.indx - frame1.indx
            if (ntimes is not None and 
                    abs(dindx) not in ntimes):
                continue
            if frame1.indx > infoloc:
                infoloc += infostep
                logging.info('autocorr %d %%', int(frame1.indx * 100.0 /len(self.frames)))
            #~ print(frame1.time, frame2.time, dt)
            locs1, locs2 = frame1.locarray, frame2.locarray
            #~ print(locs1, locs2)
            diff = (locs2 - locs1) * (2j * math.pi / lengthscale)
            exps = np.exp(diff)
            exponents[dindx].append(abs(np.mean(exps)))
            #~ for (loc1, loc2) in zip(frame1, frame2):
                #~ v1.set(*loc1)
                #~ v2.set(*loc2)
                #~ v2 -= v1
                #~ for x in (v2.X, v2.Y, v2.Z):
                    #~ exponents[dt] += cmath.exp(1j*x / lengthscale)
                #~ counts[dt] += 3
        
        logging.info('autocorr done')
        
        retdict = dict()
        for dindx, exp in exponents.items():
            dt = self.frames[dindx].time - self.frames[0].time
            retdict[dt] = np.mean(exp)
        
        return self._shelf_in(key, retdict)
    
    def _autocorr_n(self, n, lengthscale, maxcount=100):
        k = (2j * math.pi / lengthscale)
        fpairs = zip(self.frames, self.frames[n:])
        f1,f2 = fpairs[0]
        assert f2.indx - f1.indx == n
        if len(fpairs) > maxcount:
            oldlen = len(fpairs)
            skipn = len(fpairs) // maxcount
            fpairs = fpairs[::skipn]
            logging.info('skipping... %d; len: %d -> %d; maxcount %d, skip %d',
                n, oldlen, len(fpairs) ,maxcount, skipn)
        return np.mean([np.exp((f2.locarray - f1.locarray)*k) for (f1,f2) in fpairs])
    
    def Rg_autocorr_old(self, ntimes=None):
        key = 'autocorr-Rg'
        val = self._shelf_out(key)
        if ntimes is None:
            ntimes = self.frames[-1].indx - self.frames[0].indx
        if val and (len(val) >= ntimes or len(val) >= len(self.times)-1):
            return val
        elif val:
            logging.info('Rg_autocorr ignoring: %d < %d (%d)', len(val), ntimes, len(self.times))
        nframes = self.frames[-1].indx - self.frames[0].indx
        logging.debug('Rg_autocorr ntimes %d nframes %.3f',ntimes,nframes)
        ntimes = geometric(nframes, int(ntimes))
        
        ts, indxs = [f.time for f in self.frames], [f.indx for f in self.frames] 
        
        Rgs = np.array(self.Rg())
        Rgdtdict = collections.defaultdict(list)
        Rgts = zip(Rgs, ts, indxs)
        mean = float(np.mean(Rgs))
        var = float(np.mean(Rgs*Rgs) - mean*mean)
        logging.info('<Rg> = %.3f, var = %.3f, std = %.3f', mean, var, np.std(Rgs))
        logging.debug('Rg_autocorr ntimes indices: ' + str(ntimes))
        
        dtdict = dict()
        seen = []
        for (Rgt1, Rgt2) in itertools.combinations(Rgts,2):
            R1, t1, i1 = Rgt1
            R2, t2, i2 = Rgt2
            dindx = i2 - i1
            if (ntimes is not None and 
                    abs(dindx) not in ntimes):
                continue
            if t1 not in seen:
                #~ logging.info('Rg_autocorr at %.3f out of %.3f', t1, self.times[-1])
                seen.append(t1)
            Rgdtdict[dindx].append(R1*R2)
            dtdict[dindx] = t2 - t1
        
        corrdict = dict([(dtdict[dindx],float(np.mean(v))) for (dindx,v) in Rgdtdict.items()])
        #~ print("corrdict:", sorted(corrdict.items()))
        
        corrdict = dict([(dt,(v-mean*mean) / var) for (dt,v) in corrdict.items()])
        
        return self._shelf_in(key, corrdict)



def groupdicts(dlst, key):
    bigdict = dict()
    for d in dlst:
        curval = d[key]
        keylist = bigdict.get(curval, [])
        keylist.append(d)
        bigdict[curval] = keylist
    return bigdict

def filefinder(dir, *names, **args):
    import re, fpath
    from decimal import Decimal
    from namespace import Namespace
    cut = args.get('cut', 0)
    
    dir = fpath.Dir(dir)
    xchildren = [(f, f[:-1] + (f[-1].rpartition('.')[0] + '.stats'))
                    for f in dir.children() if f.extension == 'xyz']
    schildren = [(f[:-1] + (f[-1].rpartition('.')[0] + '.xyz'), f)
                    for f in dir.children() if f.extension == 'stats']
    cpairs = list(set(xchildren + schildren))
    
    regexpr = '-'.join([n + '((?:[.0-9]*)|(?:inf))' for n in names]) + r'\.xyz'
    regex = re.compile(regexpr)
    matches = [(regex.match(xf[-1]), xf, sf) for xf, sf in cpairs]
    for m, xf,sf in matches:
        if m is None:
            raise ValueError('%s could not match %s' % (regexpr, xf[-1]))
        elif None in m.groups():
            print('Groups:', m.groups())
            raise ValueError('%s could not entirely match %s' % (regexpr, xf[-1]))
    groups = sorted(      [zip(names, map(Decimal, mtch.groups()))    
                        + [('xyz',fpath.File(xf)), ('stats',fpath.File(sf))] 
                        for mtch, xf, sf in matches])
    pdicts = [dict(lst) for lst in groups]
    for p in pdicts:
        if not p['xyz'].exists(): p['xyz'] = None
        if not p['stats'].exists(): p['stats'] = None
        p['sk'] = statkeeper(p['xyz'], p['stats'] or 'auto', cut=cut, open=False)
        p.update(args)
    return [Namespace(d) for d in pdicts]

#~ import matplotlib.pyplot as plt
#~ def dplot(dct, func=plt.plot, *args, **kwargs):
    #~ ts, vs = zip(*sorted(dct.items()))
    #~ func(ts, vs, *args, **kwargs)
    
def showcorrs(fnames):
    import matplotlib.pyplot as plt
    for xyzf, statf, T in fnames:
        print('opening', xyzf, 'shelf:', statf)
        sk = statkeeper(xyzf, statf)
        gr = sk.gyradius()
        gr = sum([g for t,g in gr]) / len(gr)
        print('gyradius:', gr)
        acorr = sk.autocorr(gr, 100)
        relax=sk.relaxation(gr, 100)
        print('relaxation:', relax)
        if relax: relax = '%.0f' % relax
        else: relax = '?'
        label='T=%.3g; Rg=%.2f; t=%s' % (T, gr, relax)
        dplot(acorr, label=label, func=plt.semilogx)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    import sys
    for fname in sys.argv[1:]:
        try:
            with statkeeper(fname) as sk:
                print(fname, end='')
                sys.stdout.flush()
                t0 = sk.times[-1]
                print(', t =', t0/1e3)
        except Exception as e:
            print(e)
    exit()
    
    
    #OLD
    
    #~ logging.getLogger().setLevel(logging.DEBUG)
    
    #~ nums = ('0.2','0.5','1','3','5','10','20')
    #~ basenames = [mydir + 'data/T%s-200K' % s for s in nums]
    nums = ('0.2','0.5','1','3','5','10','20','40','60','80','100')
    basenames = [mydir + 'data2/T%s-1000K' % s for s in nums]
    
    fnames = [(bn + '.xyz', bn + '.stats', float(n)) for bn,n in zip(basenames, nums)]
    
    
    Ts = []
    relaxts = []
    for xyzf, statf, T in fnames:
        print('Running', xyzf)
        sk = statkeeper(xyzf, statf)
        #~ for k in ('bondstds','anglestds','bondmeans','anglemeans'):
        for k in ('anglemeans',):
            if k in sk.shelf: del sk.shelf[k]
        print('bondstd:', sk.bond_std(), average(sk.bond_stds()))
        print('angle std:', sk.angle_std(), average(sk.angle_stds()))
        #~ bmeans, ameans = sk.bond_means(), sk.angle_means()
        #~ print('bond mean:', bmeans[0], average(bmeans))
        #~ print('angle mean:', ameans[0], average(ameans))
        gr = sk.gyradius()
        gr = average(gr)
        print('gyradius:', gr)
        acorr = sk.autocorr(gr, 100)
        relax=sk.relaxation(gr, 100)
        print('relax:', relax)
        if relax is None: relax = 0
        Ts.append(T)
        relaxts.append(relax)
    
    
    plt.plot(Ts, relaxts)
    plt.xlabel('Temperature (unitless)')
    plt.ylabel('Relaxation time (unitless)')
    plt.show()
