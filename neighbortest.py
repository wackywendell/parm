#!/usr/bin/env python
# encoding: UTF-8

print("Importing...")

from simw import *
from math import sqrt
import sys, simpdb
from xyzfile import XYZwriter, XYZreader
from Bio.PDB import PDBParser
import numpy as np
from os.path import expanduser
from optparse import OptionParser

def average(lst):
    return sum(lst) / len(lst)

def average_squared(lst):
    return sqrt(sum([s*s for s in lst]) / len(lst))

mydir = expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1)
parser.add_option('-D', '--damping', type=float, dest='damping', default=.5)
parser.add_option('--dsteps', type=int, dest='dsteps', default=1000)
parser.add_option('-t', '--time', type=float, dest='time', default=1)
parser.add_option('-d', '--dt', type=float, dest='dt', default=.01)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=10)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test/T{T}-{t}K-{L}.xyz'))
parser.add_option('-L', '--label', dest = 'label', default=None)
parser.add_option('-X', dest='simplified', action='store_true')
parser.add_option('-n', '--neighborlist', type=float, default=1.4)
opts,args = parser.parse_args()

steps = int(opts.time * 1000 / opts.dt + .5)
bondspring = 5000
anglespring = 30
LJepsilon = 1
chargescreen = 9 # 9 angstrom Debye length
chargek = 6.954

pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'
label = opts.label if opts.label else (
    'simplified' if opts.simplified else 'neighbors')
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'),
    t=format(opts.time, '.4g'), N=opts.numres, L=label, d=opts.dt)
startxyz = XYZreader(open(opts.startfile,'r')) if opts.startfile else None
showsteps = int(float(steps) / opts.showsteps+.5)

print('steps:', steps, opts.showsteps, showsteps, showsteps*opts.dt)

#~ szs = [1.4,1.8,2.0]
szs = [1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.2,2.5,3.0,3.5,4.0]
szdict = dict((s,[]) for s in szs)

def avg(lst):
    return float(sum(lst)) / len(lst)

####################
print("Importing PDB...")
avecs = simpdb.Resvec.from_pdb('aS', pdbfile, loadfile, numchains=1)

if opts.numres > 0:
    oldavecs = avecs
    avecs = avecs[:opts.numres]
    print("%d Residues found, %d used, with %d atoms. Making structure..."
                % (len(oldavecs), len(avecs), sum(len(r) for r in avecs)))

print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))

print(len(list(simpdb.Resvec.all_bonds(avecs))), 'bonds, ',
    len(list(simpdb.Resvec.all_angles(avecs))), 'angles')

bonds = simpdb.make_bonds(avecs, bondspring)
angles = simpdb.make_angles(avecs, anglespring)
#~ charges = simpdb.make_charges(avecs, chargescreen, k=chargek)
LJ,neighbors = simpdb.make_LJ(avecs, LJepsilon, opts.neighborlist)
trackers = tvector([neighbors])
    
#~ interactions = ivector([bonds, angles, LJ, charges])
interactions = ivector([bonds, angles, LJ])

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

collec = collectionSol(opts.dt, opts.damping, opts.temp, atomgroups, interactions, trackers)
collec.seed(0)
collec.setForces()

if opts.damping == 0:
    print("Damping for", opts.dsteps, "steps")
    collec.changeT(.5, opts.temp)
    for i in range(opts.dsteps):
        collec.timestep()
    collec.changeT(opts.damping, opts.temp)

from datetime import datetime

lst = []
t = 0        
try:
    sys.stdout.flush()
    for i in range(opts.showsteps):
        starttime = datetime.now()
        curlim = (i+1) * opts.time * 1000 / opts.showsteps - opts.dt/2
        while t * opts.dt < curlim:
            collec.timestep()
            t+=1
            
        tottime = datetime.now() - starttime
        totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)
        lst.append(float(t) / totsecs)
        lstavg = avg(lst)
        std = np.std(lst)
        npoints = len(lst)
        if npoints > 1:
            print(lstavg, std, std / sqrt(npoints-1))
        else:
            print(lstavg, std)
except KeyboardInterrupt:
    exit()

dat = (avg(lst), np.std(lst), np.std(lst) / sqrt(len(lst)-1))
with open('data/ntest.txt', 'a') as f:
    print("%.1f\t%.1f" % (opts.temp, opts.neighborlist), "%.3f\t%.3f\t%.4f" % dat,
            sep='\t', file = f)
