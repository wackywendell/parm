#!/usr/bin/env python
from __future__ import print_function
import argparse, glob, os, os.path, sys

#~ import xyzstats
from statmaker import Statmaker, default_ijs
from FRETvals import distribijs

parser = argparse.ArgumentParser()
parser.add_argument('-a',dest='get_angles', action='store_true')
parser.add_argument('-r',dest='get_Rijs', action='store_true')
parser.add_argument('-J',dest='get_Rijdistribs', action='store_true')
parser.add_argument('-R',dest='get_Rijexts', action='store_true')
parser.add_argument('files', nargs='+', type=str)

opts = parser.parse_args()
        
if len(opts.files) == 1 and os.path.isdir(opts.files[0]):
    d = opts.files[0]
    opts.files = glob.glob(d + '/*.xyz')
    print('Found', len(opts.files), 'files.')

ijs = [(9,130), (33,130), (54,130), (72,130), (33,72), (9,54),
        (72,92), (54,72), (9,72), (9,33), (54,92), (92,130)]
        
def sprint(*args, **kw):
    if 'end' not in kw: kw['end'] = ''
    print(*args, **kw)
    sys.stdout.flush()

#~ def statupdate(sk):
    #~ sk.Rg()
    #~ sprint('Rg - ')
    #~ sk.N_acorr()
    #~ sk.N_acorr(.2)
    #~ sk.N_acorr(.2, .1)
    #~ sk.N_acorr(.1, .2)
    #~ sk.N_acorr(.25, .15)
    #~ sprint('N - ')
    #~ sk.temp()
    #~ sk.times
    #~ sk.energy()
    #~ sprint('TEt - ')
    #~ if opts.get_Rijs:
        #~ for i,j in ijs:
            #~ sk.Rij(i,j)
            #~ sk.ETeff(i,j,54.0)
        #~ sprint('Rij - ')
    #~ if opts.get_angles:
        #~ sk.getPhis()
        #~ sprint('phi - ')
        #~ sk.getPsis()
        #~ sprint('psi - ')

MAXPTS=100

def statupdate(stats, f):
    numpts=len(stats.Rg(f))
    sprint('Rg - ')
    stats.relax(f)
    stats.relax(f, .2)
    stats.relax(f, .2, .1)
    stats.relax(f, .2, .2)
    #stats.relax(f, .1, .2)
    #stats.relax(f, .25, .15)
    sprint('relax - ')
    assert len(stats.T(f)) - numpts <= MAXPTS
    assert len(stats.E(f)) - numpts <= MAXPTS
    assert len(stats.Times(f)) - numpts <= MAXPTS
    sprint('TEt - ')
    if opts.get_Rijexts:
        for n in [0,1,2,3,4,5,6,8]:
            rijs = stats.Rijext(f, default_ijs, n)
            #print(rijs)
            assert all([len(v) - numpts <= MAXPTS for v in rijs.values()]), (
                    "Rgs of len %d, Rijexts of len %s" % (numpts, [len(v) for v in rijs.values()]))
    if opts.get_Rijdistribs:
        rijs = stats.Rijs(f, distribijs)
        assert all([len(v) - numpts <= MAXPTS for v in rijs.values()]), (
                    "Rgs of len %d, Rijs of len %s" % (numpts, [len(v) for v in rijs.values()]))
        stats.ETeffsAS(f, 54)
        sprint('Rij - ')
        opts.get_Rijs = False
    if opts.get_Rijs or opts.get_Rijexts:
        rijs = stats.Rijs(f, default_ijs)
        assert all([len(v) - numpts <= MAXPTS for v in rijs.values()]), (
                    "Rgs of len %d, Rijs of len %s" % (numpts, [len(v) for v in rijs.values()]))
        stats.ETeffsAS(f, 54)
        sprint('Rij - ')
    if opts.get_angles:
        stats.Dihedral(f, 'phi')
        sprint('phi - ')
        stats.Dihedral(f, 'psi')
        sprint('psi - ')


statsUA = Statmaker(memoize=False, H=False)
statsAA = Statmaker(memoize=False, H=True)
for f in opts.files:
    sprint(f, ': ', end='')
    with open(f, 'r') as readf:
        l = readf.readline().strip()
    H = (l == '2016')
    if not H and l != '1013':
        raise ValueError('Length %s not understood' % l)
    st = statsAA if H else statsUA
    
    try:
        statupdate(st, f)
    except Exception as e:
        print()
        raise
        sprint('Error with',f,':', e, 'Retrying... ')
        statfs = glob.glob(f[:-4] + '*.txt.gz')
        for sf in statfs:
            os.remove(sf)
        statupdate(st, f)
    
    sprint('Done.',end='\n')
