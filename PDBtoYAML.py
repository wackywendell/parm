#!/usr/bin/env python
# encoding: UTF-8
"""
Reads a PDB and converts it to a YAML file suitable for input.

Covers only the geometry:
atoms:
- type: $residue_name
  sizes:
    $atom_name: [$mass, $size]
    ...
  xyz:
    $atom_name: [$x, $y, $z]
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
"""

import sys, os, os.path, math
from itertools import count
from optparse import OptionParser
from collections import defaultdict, OrderedDict
from simpdb import bonddict, backbonds, hydrodict, masses

from Bio.PDB import PDBParser
pdbparser = PDBParser()
import yaml2 as yaml

mydir = os.path.expanduser('~/idp/')
parser = OptionParser()
parser.add_option('--numchains', type=int, default=None)
parser.add_option('-N', '--numres', type=int, default=None)

opts, args = parser.parse_args()

aS = pdbparser.get_structure('Protein', args[0])
chains = list(aS.get_chains())[:opts.numchains]
residues = list(itertools.chain.from_iterable(
                chain.child_list for chain in chains))[:opts.numres]

atoms, bonds, angles = [], [], []

def linetotypes(l, *types, stripout='\n', delim='\t'):
    splitline = l.strip(stripout).split(delim)
    if len(splitline) != len(types):
        raise ValueError("Line %r has length %d, whereas types %r has length %d"
            % (l, len(splitline), types, len(types)))
    return [t(w) for t,w in zip(types, splitline)]


alicesizes = defaultdict(dict)
with open(mydir + 'blengths/alicesizes.tsv') as f:
    for l in f:
        if not l.strip('\n'): continue
        r,aname,s,n = linetotypes(l, str, str, float)
        alicesizes[r][aname] = s

########################################################################
# Atom sizes
for res in residues:
    resdict = OrderedDict()
    resdict['type'] = res.resname
    ressizes = alicesizes[res.resname]
    #myatoms = [a for a in residue.child_list if a.name[0] != 'H']
    resdict['sizes'] = sizes = []
    resdict['xyz'] = xyzs = []
    for a in residue.child_list:
        mass = masses[a.name[0]]
        size = ressizes[a.name] if a.name in ressizes else alicesizes['X'][a.name]
        sizes.append([a.name, mass, size])
        x,y,z = list(map(float, a.get_coord()))
        xyzs.append([a.name, x, y, z])
    atoms.append(res)

########################################################################
# Atom Bonds



