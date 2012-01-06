from __future__ import print_function

import sys, os.path, shelve, itertools, collections, math, cmath
import numpy as np
from Bio.PDB import PDBParser

import simpdb
from simw import Vec, autocorr, calc_Rg, geometric, average, average_squared
import simw as sim
from xyzfile import XYZreader
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

class statkeeper:
    def __init__(self, xyzfname, shelfname=None, pdbfname = pdbfile, loadfname = loadfile):
        self.reader = XYZreader(open(xyzfname,'r'))
        self._resfunc = lambda: simpdb.Resvec.from_pdb('aS', pdbfname, loadfname, numchains=1)
        self.collec = sim.StaticCollec(self.residues)
        if shelfname: 
            self.shelf = shelve.open(shelfname)
            
            szkey = '_size'
            if szkey in self.shelf and self.shelf[szkey] == self.reader.size():
                logging.info('shelf exists, size matches')
            else:
                if szkey not in self.shelf:
                    logging.info('New shelf')
                else: logging.info('Size does not match')
                self.shelf.clear()
                self.shelf[szkey] = self.reader.size()
                self.shelf.sync()
            #~ mkey = '_md5'
            #~ if mkey in self.shelf and self.shelf[mkey] == self.reader.md5():
                #~ logging.info('shelf exists, md5 matches')
            #~ else:
                #~ if '_md5' not in self.shelf:
                    #~ logging.info('New shelf')
                #~ else: logging.info('MD5 does not match')
                #~ self.shelf.clear()
                #~ self.shelf[mkey] = self.reader.md5()
                #~ self.shelf.sync()
        else: self.shelf = None
    
    def close(self):
        if self.shelf: self.shelf.close()
        if self.reader: self.reader.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self.close()
    
    @property
    def residues(self):
        if not hasattr(self, '_residues'):
            logging.info('building residues')
            self._residues = self._resfunc()
            logging.info('residues built')
        return self._residues
    
    @property
    def atoms(self):
        if not hasattr(self, '_atoms'):
            self._atoms = [a for r in self.residues for a in r]
        return self._atoms
    
    @property
    def frames(self):
        if not hasattr(self, '_frames'):
            logging.info('importing frames')
            self._frames = self.reader.all()
            logging.info('frames imported')
        return self._frames
    
    @property
    def times(self):
        if not hasattr(self, '_times'):
            self._times = self._shelf_out('times') or [f.time for f in self.frames]
        return self._shelf_in('times', self._times)
    
    def _shelf_out(self, key):
        if self.shelf is not None and key in self.shelf:
            logging.debug('found %s', key)
            return self.shelf[key]
        else:
            logging.info('key %s not found in shelf', key)
            return None
    
    def _shelf_in(self, key, val):
        if self.shelf is not None:
            self.shelf[key] = val
            logging.debug('putting ' + key)
        else:
            logging.debug('shelf does not exist %r' % self.shelf)
        return val
    
    def Rg(self):
        val = self._shelf_out('Rg')
        if val: return val
        Rgs = [calc_Rg(self.residues) for t in self.frames.into(self.atoms)]
            
            #~ locs = [r['CA'].x for r in self.residues]
            #~ center = sum(locs, Vec(0,0,0)) / len(locs)
            #~ locsqs = [(loc - center).sq() for loc in locs]
            #~ Rgs.append(math.sqrt(average(locsqs)))
        return self._shelf_in('Rg', Rgs)
            
    def gyradius(self):
        return self._shelf_out('gyradius') or self._shelf_in('gyradius',
            [self.collec.gyradius() for t in self.frames.into(self.atoms)])
    
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
        key = 'autocorr-Rg'
        val = self._shelf_out(key)
        if val: return val
        print('autocorr')
        
        Rgs = self.Rg()
        if not self.times:
            return ([],[])
        t0 = self.times[0]
        dts = [t - t0 for t in self.times]
        return self._shelf_in(key, (dts, autocorr(Rgs)))
    
    def relax_acorr(self):
        """Relaxation time, calculated using the autocorrelation function 
        of the Rg timeseries."""
        acorr = self.autocorr()
        #~ print('acorr:', acorr)
        ts = [(t,v) for t,v in zip(*acorr) if v < 1.0/math.e]
        if not ts:
            return None
        t0,v0 = min(ts)
        return t0
    
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
        return self._shelf_out('timespan') or self._shelf_in('timespan',
                                    self.times[-1] - self.times[0])
    
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
    
    @_use_key('bondstds')
    def bond_stds(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        bonds = simpdb.make_bonds(self.residues, 0)
        return [bonds.std_dists() for t in self.frames.into(self.atoms)]
    
    @_use_key('bondmeans')
    def bond_means(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        bonds = simpdb.make_bonds(self.residues, 0)
        return [bonds.mean_dists() for t in self.frames.into(self.atoms)]
    
    def bond_std(self):
        """Return bond distance std. dev. averaged over all atoms and frames"""
        return average_squared(self.bond_stds())
    
    @_use_key('anglestds')
    def angle_stds(self):
        """Return bond distance std. deviations as a list, averaged per frame"""
        angles = simpdb.make_angles(self.residues, 0)
        return [angles.std_dists() for t in self.frames.into(self.atoms)]
    
    @_use_key('anglemeans')
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
    

import matplotlib.pyplot as plt
def dplot(dct, func=plt.plot, *args, **kwargs):
    ts, vs = zip(*sorted(dct.items()))
    func(ts, vs, *args, **kwargs)
    
def showcorrs(fnames):
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
