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

Sample code:
parallel --eta python src/PDBtoCGYAML-sixmodes.py -s {1} -r {2} -j {3} pdb/beta-synuclein-phyre2.pdb '2>/dev/null' '>startlocs/cg/bS-six-{1}-{2}-{3}.yml' ::: KD monera hydroavg sharpcorrected aWW ::: None zeroone max minmax ::: arithmetic arithmeticzero geometric
"""

import sys, os, os.path, math
from itertools import count, chain
from optparse import OptionParser
from collections import defaultdict, OrderedDict

mydir = os.path.expanduser('~/idp/')
parser = OptionParser()
parser.add_option('-a', dest='oldangles', action='store_true')
parser.add_option('-b', dest='oldbonds', action='store_true')
parser.add_option('-d', dest='olddihedrals', action='store_true')
parser.add_option('-R', dest='LJradius', type=float, default=4.8)
parser.add_option('-s', '--hset', type='choice', default='KD',
        choices = ['monera', 'sharpcorrected','KD','aWW','hydroavg'])
parser.add_option('-r', '--hrescale', type='choice', default='None',
        choices=['None', 'max', 'minmax', 'zeroone'])
parser.add_option('-j', '--hjoin', type='choice', default='arithmetic',
        choices=['arithmetic', 'arithmeticzero','geometriczero','max','maxzero'])

opts, args = parser.parse_args()

seq = ''.join(args).strip()
seq = seq.replace(' ','')
seq = seq.replace('\n','')

from simpdb import res31, res13
import simpdb
from math import sqrt, pi

bad = False
for n,R in enumerate(seq):
    if R not in res13:
        print("Residue %n (%s) not recognized" % (n,R))
        bad = True 

if bad: exit(1)

from Bio.PDB import PDBParser
pdbparser = PDBParser()
import yaml2 as yaml

def round5(n):
    try:
        return 0. if n == 0 else round(n, int(4 - math.log10(abs(n))))
    except ValueError:
        print(n)
        raise

dihedralcoeffs = [(2.00289+0.00000j), 0.70534-0.17480j, -0.31279-0.09341j,
                                    -0.07877+0.03025j, 0.04106+0.02972j]
#[0.0+0.0j, -0.3824-0.0880j, 0.1178-0.0099j, -0.0304-0.0148j]

if opts.olddihedrals:
    dihedralcoeffs = [2.1701+0.0000j, 0.7294+0.0152j, -0.4723-0.7341j,
    -0.0022+0.0097j, 0.0272+0.1575j, -0.0295+0.0063j, -0.0298-0.0243j]
    
#dihedralcoeffs = [0,0, 0.0 + 0.2j]

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
radius, epsilon = opts.LJradius, 1.
sigcut = 2.9045 / pow(2, 1.0/6.0) #2.5

#bondlen, bondstd = 3.61174679048, 0.242683178422
bondlen, bondstd = 3.87261750303, .04588
if opts.oldbonds:
    bondlen, bondstd = 3.988, 0.024

#angmean, angstd = 2.2820824993, 0.425869565975
angmean, angstd = 2.11948895713, 0.258673706051
if opts.oldangles:
    angmean, angstd = 2.058, 0.189

# ε_ES
# comes from using 293K, angstrom, and e_charge as standard units
# (e_charge^2)/(4pi 80 electric_constant angstrom boltzmann (293K))
#                    ^ 80 comes from relative permittivity of water
chargek = 7.1288713
# note that ε_ES = chargek / 4.8 Å = 1.4852

charges = dict(ARG=1.0, LYS=1.0, HIS=.1, ASP=-1.0, GLU=-1.0)

########################################################################
# Get hydrophobicity values
assert hasattr(simpdb, opts.hset)
hset = getattr(simpdb, opts.hset)
for r in res31: assert r in hset

Hnames, Hcoeffs = zip(*sorted(hset.items()))

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

# Note: now we just use a 0 to 1 (or -1 to 1) rating, and let CG.py
# assign epsilon
# Hcoefflist = [h for h in Hcoefflist]

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
xs = [math.cos((pi-angmean)/2)*bondlen*n for n,R in enumerate(seq)]
ys = [math.sin((pi-angmean)/2)*bondlen*pow(-1,n)/2.0 for n,R in enumerate(seq)]
zs = [0.0 for n,R in enumerate(seq)]

def cmean(coords): return float(sum(coords)) / len(coords)
comx, comy, comz = cmean(xs), cmean(ys), cmean(zs)
xs = [x - comx for x in xs]
ys = [y - comy for y in ys]
zs = [z - comz for z in zs]

#~ print(xs[1]-xs[0], ys[1] - ys[0], file=sys.stderr)
#~ print(xs[2]-xs[1], ys[2] - ys[1], file=sys.stderr)

atoms = []
for n,R in enumerate(seq):
    resdict = OrderedDict()
    resdict['type'] = name = res13[R]
    resdict['atoms'] = OrderedDict()
    atomdict = OrderedDict()
    atomdict['mass'] = m = round5(masses[name])
    #atomdict['radius'] = radius
    Hindx = Hnames.index(name)
    atomdict['LJAttractFixedRepulse'] = OrderedDict(
        [('epsilons', Hcoefflists[Hindx]), ('repeps', epsilon), ('sigma', radius),
        ('indx', Hindx), ('cut', sigcut)])
    if name in charges: atomdict['Charge'] = charges[name]
    x,y,z = [round(float(c), 5) for c in (xs[n], ys[n], zs[n])]
    atomdict['xyz'] = xyzs = [x,y,z]
    
    # if res in charges: atomdict['charge'] = float(charges[res])
    resdict['atoms']['CA'] = atomdict
    atoms.append(resdict)


########################################################################
# Atom Bonds
bonds = sorted(
        [[i, 'CA', i+1, 'CA', bondlen, round5(1./bondstd/bondstd)]
                            for i in range(len(seq)-1)])

angles = sorted(
        [[i, 'CA', i+1, 'CA', i+2, 'CA', angmean, round5(1./angstd/angstd)]
                for i in range(len(seq)-2)])
dihedrals  = [[i, 'CA', i+1, 'CA', i+2, 'CA', i+3, 'CA',
                                            dihedralcos, dihedralsin, False]
                for i in range(len(seq)-3)]

########################################################################
# Basic Parameters
basics = OrderedDict()
esdict = basics['Charges'] = OrderedDict()
esdict['screeninglength'] = 9
esdict['epsilon'] = chargek # = 7.1288713

########################################################################
outdict = OrderedDict()
outdict['basics'] = basics
outdict['structure'] = atoms
outdict['bonds'] = bonds
outdict['angles'] = angles
outdict['dihedrals'] = dihedrals

yaml.safe_dump(outdict, sys.stdout, default_flow_style=None)
