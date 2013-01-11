#!/usr/bin/env python
# encoding: UTF-8

from simw import *
import simpdb
from optparse import OptionParser

parser = OptionParser("Runs a simple simulation.")

opts, args = parser.parse_args()
assert len(args) == 2
args = [a.upper() for a in args]
for arg in args:
    assert (arg in simpdb.res31 or arg in simpdb.res13), "Residue %s not understood" % arg

class Res:
    def __init__(self, anames):
        self.atomvec = atomvec([simpdb.masses[a[0]] for a in anames])
        self.names = anames
    
    def __getitem__(self, name):
        if name not in self.names:
            raise KeyError(name)
        idx = self.names.index(name)
        return self.atomvec[idx]
    
    def __iter__(self):
        return zip(self.names, self.atomvec)

def loadfile(arg, vstart = Vec(0,0,0)):
    fname = 'pdb/locs-' + simpdb.res13.get(arg, arg) + '.tsv'
    atomlocs = {}
    for l in open(fname, 'r'):
        if not l.strip(): continue
        atomname, x,y,z = l.strip().split('\t')
        atomlocs[atomname] = Vec(float(x), float(y), float(z)) + vstart
    anames = sorted(atomlocs)
    res = Res(anames)
    for aname, atom in res:
        atom.x = atomlocs[aname]

## Make bonds


## Make Lennard-Jones (sigma is the minimum, and also the cutoff)
LJsizes = {}
for l in open('pdb/alicesizes.tsv'):
    if not l.strip(): continue
    r,aname,s = l.strip().split('\t')
    LJsizes[(r,aname)] = float(s)

maxsigma = max(LJsizes.values()) * 2
outerradius = maxsigma * 1.8
neighbors = neighborlist(innerradius, outerradius)
LJ = LJgroup(neighbors)



residues = [loadfile(arg) for arg in args]

maxsigma = max(LJdict.values()) * factor
innerradius = maxsigma
outerradius = innerradius * neighborcutoff
neighbors = neighborlist(innerradius, outerradius)
LJ = LJgroup(neighbors)

