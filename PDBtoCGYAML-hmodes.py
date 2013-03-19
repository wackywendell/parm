#!/usr/bin/env python
# encoding: UTF-8
"""
Reads a PDB and converts it to a YAML file suitable for input.

Covers only the geometry:
atoms:
- type: $residue_name
  atoms:
    - name: $name
      mass: $mass
      size: $size
      radius: $radius
      xyz: [$x, $y, $z]
    ...
- ...
bonds:
- [$residue_number, $atom_name, $residue_number, 
                $atom_name, $length, $std.dev.]
- ...
angles:
- [$residue_number, $atom_name, $residue_number, $atom_name, 
                $residue_number, $atom_name, $angle, $std.dev.]
- ...
dihedrals:
- [$residue_number, $atom_name, $residue_number, $atom_name, 
                $residue_number, $atom_name, $residue_number, $atom_name, 
                [coscoeff1, coscoeff2, ..], [sincoeff1, sincoeff2, ...]]
- ...
"""

import sys, os, os.path, math
from itertools import count, chain
from optparse import OptionParser
from collections import defaultdict, OrderedDict
from simpdb import res31
import simpdb
from math import sqrt

from Bio.PDB import PDBParser
pdbparser = PDBParser()
import yaml2 as yaml

mydir = os.path.expanduser('~/idp/')
parser = OptionParser()
parser.add_option('-s', '--hset', type='choice', default='hydroavg',
        choices = ['monera', 'sharpcorrected','KD','aWW','hydroavg'])
parser.add_option('-r', '--hrescale', type='choice', default='None',
        choices=['None', 'max', 'minmax', 'zeroone'])
parser.add_option('-j', '--hjoin', type='choice', default='arithmetic',
        choices=['arithmetic', 'arithmeticzero','geometriczero','max','maxzero'])

opts, args = parser.parse_args()

def round5(n):
    try:
        return 0. if n == 0 else round(n, int(4 - math.log10(abs(n))))
    except ValueError:
        print(n)
        raise
    

aS = pdbparser.get_structure('Protein', args[0])
chains = list(aS.get_chains())#[:opts.numchains]
residues = list(chain.from_iterable(
                rchain.child_list for rchain in chains))#[:opts.numres]

atoms, bonds, angles, dihedrals = [], [], [], []

def linetotypes(l, *types, stripout='\n', delim='\t'):
    splitline = l.strip(stripout).split(delim)
    if len(splitline) != len(types):
        raise ValueError("Line %r has length %d, whereas types %r has length %d"
            % (l, len(splitline), types, len(types)))
    return [t(w) for t,w in zip(types, splitline)]

dihedralcoeffs = [2.1701+0.0000j, 0.7294+0.0152j, -0.4723-0.7341j,
    -0.0022+0.0097j, 0.0272+0.1575j, -0.0295+0.0063j, -0.0298-0.0243j]
# Where that corresponds to a potential of the form 
# V(theta) = (sum over n from 0 to 6) Re(A_n) * cos(n theta) + Im(A_n) * sin(n theta).
dihedralcos = [n.real for n in dihedralcoeffs]
dihedralsin = [n.imag for n in dihedralcoeffs]

masses = {'ALA': 71.0779, 'ARG': 157.19362, 'ASN': 114.10264,
 'ASP': 114.07946, 'CYS': 103.1429, 'GLN': 128.12922,
 'GLU': 128.10604, 'GLY': 57.05132, 'HIS': 137.13928,
 'ILE': 113.15764, 'LEU': 113.15764, 'LYS': 129.18022,
 'MET': 131.19606, 'PHE': 147.17386, 'PRO': 97.11518,
 'SER': 87.0773, 'THR': 101.10388, 'TRP': 186.2099,
 'TYR': 163.17326, 'VAL': 99.13106}

#~ radii =  {'ALA': 4.2054, 'ARG': 4.2455, 'ASN': 3.8945, 'ASP': 4.056,
         #~ 'CYS': 3.7313, 'GLN': 4.8582, 'GLU': 4.3658, 'GLY': 3.7121,
         #~ 'HIS': 5.3274, 'ILE': 3.5582, 'LEU': 5.0538, 'LYS': 4.2056,
         #~ 'MET': 5.4487, 'PHE': 4.448, 'PRO': 4.0976, 'SER': 4.346,
         #~ 'THR': 4.1937, 'TRP': 5.5958, 'TYR': 4.949, 'VAL': 4.0486}

#radius = 3.185

#from pairsCOM:
#radius, n, epsilon = 4.21771939946, 1.27173585374, 1.88531484851
radius, n, epsilon = 4., 2., 1.
sigcut = 16. / radius

BBmass = 2.0
BBradius = 1.5

bondlen, bondstd = 3.988, 0.024
BBbondlen, BBbondstd = bondlen /2, bondstd / math.sqrt(2.0)
angmean, angstd = 2.058, 0.189
BBangmean, BBangstd = math.pi, 0.01

# comes from using 293K, angstrom, and e_charge as standard units
# (e_charge^2)/(4pi 80 electric_constant angstrom boltzmann (293K))
#                    ^ 80 comes from relative permittivity of water

chargek = 7.1288713
charges = dict(ARG=1.0, LYS=1.0, HIS=.1, ASP=-1, GLU=-1.0)

########################################################################
# Get hydrophobicity values
assert hasattr(simpdb, opts.hset)
hset = getattr(simpdb, opts.hset)
for r in res31: assert r in hset

Hnames, Hcoeffs = zip(*sorted(hset.items()))
# uses σ=radius for the definition of α to work; assumes α=1, 
# lets CGbb.py figure that out
Hcoeffs = [h * chargek/radius for h in Hcoeffs]

#-----------------------------------------------------------------------
# Rescale
if opts.hrescale.lower() == 'none':
    Hcoefflist = Hcoeffs
elif opts.hrescale.lower() == 'max':
    mx = max([abs(h) for h in Hcoeffs])
    Hcoefflist = [h / mx for h in Hcoeffs]
elif opts.hrescale.lower() == 'minmax':
    mn, mx = min(Hcoeffs), max(Hcoeffs)
    Hcoefflist = [(h - mn) / (mx - mn) for h in Hcoeffs]
elif opts.hrescale.lower() == 'zeroone':
    mx = max([abs(h) for h in Hcoeffs])
    Hcoefflist = [(h + mx) / (2*mx) for h in Hcoeffs]
else:
    raise NotImplementedError("Unknown scaling format %r" % opts.hrescale)

#-----------------------------------------------------------------------
# Join into a table

if opts.hjoin.lower() == 'arithmetic':
    Hcoefflists = [[(H1 + H2)/2 for H1 in Hcoefflist] for H2 in Hcoefflist]
elif opts.hjoin.lower() == 'arithmeticzero':
    Hcoefflists = [[max((H1 + H2)/2,0) for H1 in Hcoefflist] for H2 in Hcoefflist]
elif opts.hjoin.lower() == 'geometriczero':
    Hcoefflists = [[sqrt(max(H1,0) * max(H2,0)) for H1 in Hcoefflist] for H2 in Hcoefflist]
elif opts.hjoin.lower() == 'max':
    Hcoefflists = [[max(H1,H2) for H1 in Hcoefflist] for H2 in Hcoefflist]
elif opts.hjoin.lower() == 'maxzero':
    Hcoefflists = [[max(H1,H2,0) for H1 in Hcoefflist] for H2 in Hcoefflist]
else:
    raise NotImplementedError("Unknown joining format %r" % opts.hrescale)

Hcoefflists = [[round5(H) for H in lst] for lst in Hcoefflists]

########################################################################
# Atom sizes
xs, ys, zs = zip(*[res['CA'].get_coord() for res in residues])
def cmean(coords): return float(sum(coords)) / len(residues)
comx, comy, comz = cmean(xs), cmean(ys), cmean(zs)

for res in residues:
    resdict = OrderedDict()
    resdict['type'] = name = res.resname
    resdict['atoms'] = OrderedDict()
    atomdict = OrderedDict()
    atomdict['mass'] = m = round5(masses[name] - BBmass)
    #atomdict['radius'] = radius
    Hindx = Hnames.index(name)
    atomdict['LJish'] = OrderedDict(
        [('epsilons', Hcoefflists[Hindx]), ('repeps', epsilon), ('sigma', radius),
        ('n', n), ('indx', Hindx), ('cut', sigcut)])
    atomdict['LJAttractRepulse'] = OrderedDict(
        [('epsilons', [-1.,0.]), ('sigma', BBradius),
                                ('indx', 1), ('cut', 1.0)])
                                 # we cut at r = sigma, its repulsive
                                 # eps is negative to be repulsive
    #myatoms = [a for a in residue.child_list if a.name[0] != 'H']
    if name in charges: atomdict['Charge'] = charges[name]
    x,y,z = [round(float(n) - c, 5) for n,c in zip(res['CA'].get_coord(), (comx, comy, comz))]
    atomdict['xyz'] = xyzs = [x,y,z]
    
    # if res in charges: atomdict['charge'] = float(charges[res])
    resdict['atoms']['CA'] = atomdict
    atoms.append(resdict)

## Add backbone atoms
for r1, r2 in zip(atoms, atoms[1:]):
    x1,y1,z1 = r1['atoms']['CA']['xyz']
    x2,y2,z2 = r2['atoms']['CA']['xyz']
    x,y,z = (x1+x2)/2., (y1+y2)/2., (z1+z2)/2.
    atomdict = OrderedDict()
    atomdict['mass'] = BBmass
    #atomdict['radius'] = BBradius
    atomdict['LJAttractRepulse'] = OrderedDict(
        [('epsilons', [-1.,-1.]), ('sigma', BBradius),
                                ('indx', 0), ('cut', 1.0)])
    atomdict['xyz'] = ([round(c,5) for c in (x, y, z)])
    # no LJish, we don't do that for the backbone atoms!
    # no charge, not on the backbones
    r1['atoms']['BB'] = atomdict
    

########################################################################
# Atom Bonds
bonds = sorted(
        [[i, 'CA', i, 'BB', BBbondlen, round5(1./BBbondstd/BBbondstd)]
                            for i in range(len(residues)-1)] +
        [[i, 'BB', i+1, 'CA', BBbondlen, round5(1./BBbondstd/BBbondstd)]
                            for i in range(len(residues)-1)],
            key=lambda bd: (bd[0], bd[2], bd[1], bd[3]))

angles = sorted(
        [[i, 'BB', i+1, 'CA', i+1, 'BB', angmean, round5(1./angstd/angstd)]
                for i in range(len(residues)-2)] +
        [[i, 'CA', i, 'BB', i+1, 'CA', BBangmean, round5(1./BBangstd/BBangstd)]
                for i in range(len(residues)-1)],
            key=lambda ang: (ang[0], ang[2], ang[4], ang[1], ang[3], ang[5]))
dihedrals  = [[i, 'CA', i+1, 'CA', i+2, 'CA', i+3, 'CA',
                                            dihedralcos, dihedralsin, False]
                for i in range(len(residues)-3)]

########################################################################
# Basic Parameters
basics = OrderedDict()
esdict = basics['Charges'] = OrderedDict()
esdict['screeninglength'] = 9
esdict['epsilon'] = chargek / 4# = 7.1288713

########################################################################
outdict = OrderedDict()
outdict['basics'] = basics
outdict['structure'] = atoms
outdict['bonds'] = bonds
outdict['angles'] = angles
outdict['dihedrals'] = dihedrals

yaml.safe_dump(outdict, sys.stdout, default_flow_style=None)
