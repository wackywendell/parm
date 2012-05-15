#!/usr/bin/env python
from __future__ import print_function

import argparse, glob, os, os.path, sys

import xyzstats

parser = argparse.ArgumentParser()
parser.add_argument('-a',dest='get_angles', action='store_false')
parser.add_argument('-r',dest='get_Rijs', action='store_false')
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
    

for f in opts.files:
    sprint(f, ': ', end='')
    with open(f, 'r') as readf:
        l = readf.readline().strip()
    H = (l =='2016')
    if not H and l != '1013':
        raise ValueError('Length %s not understood' % l)
    with xyzstats.statkeeper(f, H=H) as sk:
        sk.Rg()
        sprint('Rg - ')
        sk.N_acorr()
        sk.N_acorr(.2)
        sk.N_acorr(.2, .1)
        sk.N_acorr(.1, .2)
        sk.N_acorr(.25, .15)
        sprint('N - ')
        sk.temp()
        sk.times
        sk.energy()
        sprint('TEt - ')
        if opts.get_Rijs:
            for i,j in ijs:
                sk.Rij(i,j)
                sk.ETeff(i,j,54.0)
            sprint('Rij - ')
        if opts.get_angles:
            sk.getPhis()
            sprint('phi - ')
            sk.getPsis()
            sprint('psi - ')
            
    sprint('Done.',end='\n')
