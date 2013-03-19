#!/usr/bin/env python
# encoding: UTF-8

import sys, os, os.path
from itertools import count
from optparse import OptionParser

mydir = os.path.expanduser('~/idp/')

def main():
    parser = OptionParser()
    opts, args = parser.parse_args()
    res, = args
    
    sizelocs = getsizelocs(res)
    bonds = getbonds(res)
    angles = getangles(res, bonds)
    
    import yaml
    yaml.safe_dump(dict(atoms=sizelocs, bonds=bonds, angles=angles,
                            hydrophobicities=fullepsilontable
                        ),
                    sys.stdout)

########################################################################
## Sizes and Locations
masses = dict(C=12.0107, H=1.00794, N=14.0067, O=15.9994, S=32.065)
# ACENOVH
# 0123456
atomindices = dict(A=0,B=0,C=1,E=2,G=0,M=3,N=3,O=4,Q=4,S=1,V=5,H=6)

alicesizes, levittindex = {}, {}
with open(mydir + 'blengths/alicelevittsizes.tsv') as f:
    for l in f:
        if not l.strip(): continue
        r,aname,s,n = l.strip().split('\t')
        rd = alicesizes.get(r, {})
        rd[aname] = float(s)
        alicesizes[r] = rd
        rd = levittindex.get(r, {})
        rd[aname] = atomindices[n]
        levittindex[r] = rd

def getsizelocs(res):
    locs = {}
    with open(mydir + ('pdb/locs-%s.tsv' % res)) as f:
        for l in f:
            a,x,y,z = l.strip().split('\t')
            locs[a] = tuple(map(float, (x,y,z)))
    #~ with open(mydir + ('pdb/locs-BACKBONE-ENDS.tsv')) as f:
        #~ for l in f:
            #~ a,x,y,z = l.strip().split('\t')
            #~ locs[a] = tuple(map(float, (x,y,z)))
    
    finallocs = {}
    sizes = {}
    for a in locs:
        finallocs[a] = [round(n*1e6)/1e6 for n in locs[a]]
        s = alicesizes[res][a] if a in alicesizes[res] else alicesizes['X'][a]
        n = levittindex[res][a] if a in levittindex[res] else levittindex['X'][a]
        m = masses[a[0]]
        sizes[a] = [m,s,n]
    
    return [{'type':res, 'sizes':sizes, 'xyz':finallocs}]

########################################################################
## Hydrophobicity

def gettable(fname):
    with open(fname, 'r') as f:
        return [[float(n) for n in l.strip().split('\t')] for l in f]
Atable = gettable('blengths/levitt-Avalues.tsv')
Btable = gettable('blengths/levitt-Bvalues.tsv')
epsilontable = gettable('blengths/levitt-epsilons.tsv')
sigmatable = gettable('blengths/levitt-sigmas.tsv')

fullepsilontable = (
        [[round(n*1e6)/1e6 if n > 0 else -1.0 for n in r] 
            + [-1.0] for r in epsilontable
        ] + 
        [[-1.0]*(len(epsilontable[0])+1)]
    )



########################################################################
## Bonds

bonddict = {
    'GLY' : [],
    'ALA' : [('CA','CB')],
    'SER' : [('CA','CB'), ('CB','OG')],
    'CYS' : [('CA','CB'), ('CB','SG')],
    'MET' : [('CA','CB'), ('CB','CG'), ('CG','SD'), ('SD','CE')],
    'LYS' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','CE'), ('CE','NZ')],
    'VAL' : [('CA','CB'), ('CB','CG1'), ('CB','CG2')],
    'THR' : [('CA','CB'), ('CB','CG2'), ('CB','OG1')],
    'ILE' : [('CA','CB'), ('CB','CG1'), ('CB','CG2'), ('CG1','CD1')],
    'LEU' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG','CD2')],
    'ASP' : [('CA','CB'), ('CB','CG'), ('CG','OD1'), ('CG','OD2')],
    'ASN' : [('CA','CB'), ('CB','CG'), ('CG','OD1'), ('CG','ND2')],
    'GLU' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','OE1'), ('CD','OE2')],
    'GLN' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','OE1'), ('CD','NE2')],
    'ARG' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD','NE'), 
                          ('NE','CZ'), ('CZ', 'NH1'), ('CZ', 'NH2')],
    'PRO' : [('CA','CB'), ('CB','CG'), ('CG','CD'), ('CD', 'N')],
    'HIS' : [('CA','CB'), ('CB','CG'), ('CG','ND1'), ('CG', 'CD2'),
                          ('ND1', 'CE1'), ('CD2', 'NE2'), ('CE1', 'NE2')],
    'PHE' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), ('CE2', 'CZ')],
    'TYR' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'CE1'), ('CD2', 'CE2'), ('CE1', 'CZ'), 
             ('CE2', 'CZ'), ('CZ', 'OH')],
    'TRP' : [('CA','CB'), ('CB','CG'), ('CG','CD1'), ('CG', 'CD2'),
             ('CD1', 'NE1'), ('CD2', 'CE2'), ('CD2', 'CE3'), ('NE1', 'CE2'), 
             ('CE2', 'CZ2'), ('CE3', 'CZ3'), ('CZ2', 'CH2'), ('CZ3', 'CH2')]
    }

backbonds = [('N','CA'), ('CA','C'),('C','O')]
pairbonds = [('C','N')]
hydrodict = {
'ALA' : [('H', 'N'), ('HA', 'CA'), ('HB1', 'CB'), ('HB2', 'CB'), ('HB3', 'CB')],
'ARG' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD2', 'CD'), ('HD3', 'CD'), ('HE', 'NE'), ('HG2', 'CG'), ('HG3', 'CG'), ('HH11', 'NH1'), ('HH12', 'NH1'), ('HH21', 'NH2'), ('HH22', 'NH2')],
'ASN' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD21', 'ND2'), ('HD22', 'ND2')],
'ASP' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB')],
'CYS' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HG', 'SG')],
'GLN' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HE21', 'NE2'), ('HE22', 'NE2'), ('HG2', 'CG'), ('HG3', 'CG')],
'GLU' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HG2', 'CG'), ('HG3', 'CG')], # ('HE2', 'OE2')
'GLY' : [('H', 'N'), ('HA2', 'CA'), ('HA3', 'CA')],
'HIS' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD2', 'CD2'), ('HE1', 'CE1'), ('HE2', 'NE2')], # ('HD1', 'ND1'), 
'ILE' : [('H', 'N'), ('HA', 'CA'), ('HB', 'CB'), ('HD11', 'CD1'), ('HD12', 'CD1'), ('HD13', 'CD1'), ('HG12', 'CG1'), ('HG13', 'CG1'), ('HG21', 'CG2'), ('HG22', 'CG2'), ('HG23', 'CG2')],
'LEU' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD11', 'CD1'), ('HD12', 'CD1'), ('HD13', 'CD1'), ('HD21', 'CD2'), ('HD22', 'CD2'), ('HD23', 'CD2'), ('HG', 'CG')],
'LYS' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD2', 'CD'), ('HD3', 'CD'), ('HE2', 'CE'), ('HE3', 'CE'), ('HG2', 'CG'), ('HG3', 'CG'), ('HZ1', 'NZ'), ('HZ2', 'NZ'), ('HZ3', 'NZ')],
'MET' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HE1', 'CE'), ('HE2', 'CE'), ('HE3', 'CE'), ('HG2', 'CG'), ('HG3', 'CG')],
'PHE' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD1', 'CD1'), ('HD2', 'CD2'), ('HE1', 'CE1'), ('HE2', 'CE2'), ('HZ', 'CZ')],
'PRO' : [('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD2', 'CD'), ('HD3', 'CD'), ('HG2', 'CG'), ('HG3', 'CG')],
'SER' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HG', 'OG')],
'THR' : [('H', 'N'), ('HA', 'CA'), ('HB', 'CB'), ('HG1', 'OG1'), ('HG21', 'CG2'), ('HG22', 'CG2'), ('HG23', 'CG2')],
'TRP' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD1', 'CD1'), ('HE1', 'NE1'), ('HE3', 'CE3'), ('HH2', 'CH2'), ('HZ2', 'CZ2'), ('HZ3', 'CZ3')],
'TYR' : [('H', 'N'), ('HA', 'CA'), ('HB2', 'CB'), ('HB3', 'CB'), ('HD1', 'CD1'), ('HD2', 'CD2'), ('HE1', 'CE1'), ('HE2', 'CE2'), ('HH', 'OH')],
'VAL' : [('H', 'N'), ('HA', 'CA'), ('HB', 'CB'), ('HG11', 'CG1'), ('HG12', 'CG1'), ('HG13', 'CG1'), ('HG21', 'CG2'), ('HG22', 'CG2'), ('HG23', 'CG2')]
}

def getbonds(res):
    bonds = backbonds + bonddict[res] + hydrodict[res]
    #~ bonds = [b for b in bonds if 'H' not in b]
    #~ bonds.extend([('H1','N'), ('H2','N'), ('H3','N'), ('C','OXT')])
    
    
    with open(mydir + 'blengths/pdbbonds.tsv', 'r') as f:
        blist = []
        for l in f:
            r1,a1,r2,a2,x0,k = l.strip().split('\t')
        
            if (a1,a2) in bonds and r1 in (res, 'X1') and r2 in (res, 'X1'):
                blist.append([0, a1, 0, a2, float(x0), float(k)])
            
        return blist

########################################################################
# Angles

def getangles(res, bonds):
    alist = []
    minbonds = [{a1,a2} for (r1,a1,r2,a2,x0,k) in bonds]
    with open(mydir + 'blengths/pdbangles.tsv', 'r') as f:
        for l in f:
            r1,a1,r2,a2,r3,a3,x0,k = l.strip().split('\t')
            if r1 not in (res, 'X1'): continue
            if r2 not in (res, 'X1'): continue
            if r3 not in (res, 'X1'): continue
            if {a1,a2} not in minbonds: continue
            if {a2,a3} not in minbonds: continue
            alist.append([0, a1, 0, a2, 0, a3, float(x0), float(k)])
            
            #~ else:
                #~ print(r1,a1,r2,a2,r3,a3,
                    #~ {a1,a2} in minbonds, {a2,a3} in minbonds,
                    #~ r1 in (res, 'X1'), r2 in (res, 'X1'), r3 in (res, 'X1'),
                    #~ file=sys.stderr)
        
    return alist

if __name__ == "__main__":
    main()
