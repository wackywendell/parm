# encoding: UTF8

from Bio.PDB import PDBParser, Vector
from Bio.PDB.Residue import Residue as _Residue
from Bio.PDB.Atom import Atom as _Atom
from simw import *
import itertools
from collections import defaultdict
import math

# adapted from IUPAC standards, available at
# http://www.chem.qmul.ac.uk/iupac/misc/ppep4.html#400

from os.path import expanduser
mydir = expanduser('~/idp/')
defaultpdbfile  = mydir + 'pdb/aS.pdb'
defaultloadfile = mydir + 'blengths/stats-H.pkl'

res31 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
    HIS='H', ILE='I', LYS='K', LEU='L', MET='M', ASN='N', PRO='P',
    GLN='Q', ARG='R', SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
res13 = {v:k for k,v in res31.items()}

masses = dict(C=12.0107, H=1.00794, N=14.0067, O=15.9994, S=32.065)

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

#first is ph 2.0, second is ph 7.0
hydroindex2={'ALA': 0.735, 'ARG': 0.37, 'ASN': 0.295, 'ASP': 0.41, 
                'CYS': 0.76, 'GLN': 0.41, 'GLU': 0.54, 'GLY': 0.5, 
                'HIS': 0.29, 'ILE': 1, 'LEU': 1, 'LYS': 0.315, 
                'MET': 0.87, 'PHE': 0.96, 'PRO': 0.27, 'SER': 0.465, 
                'THR': 0.565, 'TRP': 0.92, 'TYR': 0.745, 'VAL': 0.895}
hydroindex7={'ALA': 0.705, 'ARG': 0.43, 'ASN': 0.36, 'ASP': 0.225, 
                'CYS': 0.745, 'GLN': 0.45, 'GLU': 0.345, 'GLY': 0.5, 
                'HIS': 0.54, 'ILE': 0.995, 'LEU': 0.985, 'LYS': 0.385, 
                'MET': 0.87, 'PHE': 1, 'PRO': 0.27, 'SER': 0.475, 
                'THR': 0.565, 'TRP': 0.985, 'TYR': 0.815, 'VAL': 0.88}

monera = hydroindex7
# Original Monera scale
moneraneg = {k : (v * 2 - 1) for k,v in hydroindex7.items()}
monerascale = {k : (v - min(monera.values())) / 2 for k,v in monera.items()}
moneranorm = {k : v / max(monerascale.values()) for k,v in monerascale.items()}
                
                
hydroindex3 = {'ALA': 0.729, 'ARG': 0.382, 'ASN': 0.308, 'ASP': 0.373,
                'CYS': 0.757, 'GLN': 0.418, 'GLU': 0.501, 'GLY': 0.5,
                'HIS': 0.34,  'ILE': 0.999, 'LEU': 0.997, 'LYS': 0.329,
                'MET': 0.87, 'PHE': 0.968, 'PRO': 0.27, 'SER': 0.467,
                'THR': 0.565, 'TRP': 0.933, 'TYR': 0.759, 'VAL': 0.892}

#~ hydroindex={'CYS': (52, 49), 'ILE': (100, 99), 'GLN': (-18, -10), 
#~ 'VAL': (79, 76), 'LYS': (-37, -23), 'GLY': (0, 0), 'PRO': (-46, -46),
#~ 'ASP': (-18, -55), 'THR': (13, 13), 'PHE': (92, 100), 'ALA': (47, 41),
#~ 'MET': (74, 74), 'HIS': (-42, 8), 'GLU': (8, -31), 'LEU': (100, 97),
#~ 'ARG': (-26, -14), 'TRP': (84, 97), 'ASN': (-41, -28), 'TYR': (49,63),
#~ 'SER': (-7, -5)}
# see Carl's email; I don't know where these came from.
# Looks like the image from 
# http://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html

# Compare to Nozack, and got that on a gly=0, top=100, it should be 
# shifted to gly=-14, top = 100-14

# K. A. Sharp, A. Nicholls, R. Friedman, and B. Honig, Biochemistry 30, 9686 (1991).
# from Sharp et al., unnormalized, no proline (Table II, uncorrected c-> w):
{'ALA': 1.81, 'ARG': -14.92, 'ASN': -6.64, 'ASP': -8.72, 'CYS': 1.28,
 'GLN': -5.54, 'GLU': -6.81, 'GLY': 0.94, 'HIS': -4.66, 'ILE': 4.92,
 'LEU': 4.92, 'LYS': -5.55, 'MET': 2.35, 'PHE': 2.98, 'SER': -3.4,
 'THR': -2.57, 'TRP': 2.33, 'TYR': -0.14, 'VAL': 4.04}
# fit to Monera: m=0.04163418,  b=0.72607708
# that gives proline a value of -10.95
# or without asp, gly, which were clear outliers:
# m = 0.057, b=0.743
# gives Proline -8.30
# gives Asp -9.09
# gives Gly -4.27

# K. A. Sharp, A. Nicholls, R. Friedman, and B. Honig, Biochemistry 30, 9686 (1991).
# std. dev. is 1. Pro averaged from Kyte and Doolittle and Monera.
# from Sharp et al., unnormalized, no proline (Table II, uncorrected c-> w):
sharpuncorrected={
 'ALA': 0.336, 'ARG': -2.77, 'ASN': -1.233, 'ASP': -1.619, 'CYS': 0.238,
 'GLN': -1.029, 'GLU': -1.264, 'GLY': 0.175, 'HIS': -0.865, 'ILE': 0.913,
 'LEU': 0.913, 'LYS': -1.03, 'MET': 0.436, 'PHE': 0.553, 'PRO': -0.703,
 'SER': -0.631, 'THR': -0.477, 'TRP': 0.433, 'TYR': -0.026, 'VAL': 0.75}

# K. A. Sharp, A. Nicholls, R. Friedman, and B. Honig, Biochemistry 30, 9686 (1991).
# std. dev. is 1. Pro averaged from Kyte and Doolittle and Monera.
# from Sharp et al., normalized, no proline (Table II, corrected c-> w):
 
sharpcorrected={
 'ALA': 0.519, 'ARG': -2.221, 'ASN': -0.898, 'ASP': -1.28,
 'CYS': 0.506, 'GLN': -0.601, 'GLU': -0.851, 'GLY': 0.264,
 'HIS': -0.443, 'ILE': 1.382, 'LEU': 1.409, 'LYS': -0.506,
 'MET': 0.873, 'PHE': 1.065, 'PRO': -0.703, 'SER': -0.409,
 'THR': -0.171, 'TRP': 1.029, 'TYR': 0.505, 'VAL': 1.129}


KD = {'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
 'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2}

KD1 = dict([(k, v / max(KD.values())) for k,v in list(KD.items())])
KDscaled = dict([(k, (v - min(KD.values())) / 2) for k,v in list(KD.items())])
KDnorm = dict([(k, v / max(KDscaled.values())) for k,v in list(KDscaled.items())])

aWW = {'ALA': 0.67, 'ARG': -0.5, 'ASN': -3.64, 'ASP': -1.81,
 'CYS': 1.71, 'GLN': -0.77, 'GLU': -0.85, 'GLY': 1.25,
 'HIS': -1.15, 'ILE': 0.46, 'LEU': 0.71, 'LYS': 0.02,
 'MET': -0.46, 'PHE': -0.25, 'PRO': -3.63, 'SER': -0.11,
 'THR': 1.12, 'TRP': -0.14, 'TYR': -2.8, 'VAL': 2.09}

aWWscaled = dict([(k, (v - min(aWW.values())) / 2) for k,v in list(aWW.items())])
aWWnorm = dict([(k, v / max(aWWscaled.values())) for k,v in list(aWWscaled.items())])

# averages over several different scales.
hydroavg = {'ALA': 0.328, 'ARG': -1.601, 'ASN': -0.878,
             'ASP': -1.335, 'CYS': 0.54, 'GLN': -0.75,
             'GLU': -1.13, 'GLY': -0.082, 'HIS': -0.436,
             'ILE': 1.307, 'LEU': 1.213, 'LYS': -1.126,
             'MET': 0.836, 'PHE': 1.17, 'PRO': -0.46,
             'SER': -0.385, 'THR': -0.211, 'TRP': 0.851,
             'TYR': 0.279, 'VAL': 0.989}

hydroavgstd = dict(
ARG=0.952,ASP=0.471,GLU=0.552,LYS=0.500,ASN=0.278,
GLN=0.342,PRO=0.444,HIS=0.425,SER=0.244,THR=0.274,
GLY=0.437,TYR=0.511,ALA=0.422,CYS=0.390,MET=0.442,
TRP=0.673,VAL=0.400,PHE=0.429,LEU=0.352,ILE=0.388)

hydroavgscale = dict([(k, (v - min(hydroavg.values())) / 2) for k,v in list(hydroavg.items())])
hydroavgnorm = dict([(k, v / max(hydroavgscale.values())) for k,v in list(hydroavgscale.items())])


# For Ph 2:
# Sereda, Terrance J., Colin T. Mant, Frank D. Sönnichsen, and Robert S. Hodges. “Reversed-phase Chromatography of Synthetic Amphipathic Α-helical Peptides as a Model for Ligand/receptor Interactions Effect of Changing Hydrophobic Environment on the Relative Hydrophilicity/hydrophobicity of Amino Acid Side-chains.” Journal of Chromatography A 676, no. 1 (July 29, 1994): 139–153.
# For Ph 7:
# Monera, Oscar D., Terrance J. Sereda, Nian E. Zhou, Cyril M. Kay, and Robert S. Hodges. “Relationship of Sidechain Hydrophobicity and Α-helical Propensity on the Stability of the Single-stranded Amphipathic Α-helix.” Journal of Peptide Science 1, no. 5 (1995): 319–329.



def add_first_bonds(bonds_so_far, H=True):
    if ('H','N') in bonds_so_far:
        bonds_so_far.remove(('H','N'))
    if H: bonds_so_far.extend([('H1','N'), ('H2','N')])
    if ('CD', 'N') not in bonds_so_far and H:
        bonds_so_far.extend([('H3','N')])
    return bonds_so_far
    
def add_last_bonds(bonds_so_far, H=True):
    bonds_so_far.extend([('C','OXT')])
    return bonds_so_far
    
def get_bonds(resname, lastres=None, nextres=None, H=True):
    """Yields quadruplets of (r1 idx, atom1 name, r2 idx, atom2 name).
    
    r1 idx = 0 for prev, 1 for this."""
    bonds = bonddict[resname] + backbonds
    if H:
        bonds.extend(hydrodict[self.resname])
    if lastres is None: bonds = add_first_bonds(bonds, H=H)
    if nextres is None: bonds = add_last_bonds(bonds, H=H)
    bondsnamed = [(1, a1, 1, a2) for a1,a2 in bonds]
    if nextres is not None: bondsnamed = (
                bondsnamed + [(0, a1, 1, a2) for a1,a2 in bonds])
    return bondsnamed

def get_angles(resname, lastres=None, nextres=None, H=True):
    """Yields sextuplets of (r1 idx, atom1 name, r2 idx, atom2 name, r3 idx, atom3 name).
    
    r1 idx = 0 for prev, 1 for this."""
    self.load_data(f)
    myangles = self.resangles[self.resname]
    backangles = self.backangles
    dlist = defaultdict(list)
    for a1,a2,d,s in self.get_bonds(lastres, nextres, f):
        dlist[(a1.group, a1.name)].append(a2)
        dlist[(a2.group, a2.name)].append(a1)
    for center, pairlist in dlist.items():
        cgroup, cname = center
        center = cgroup[cname]
        #print('get_angles:',self.resname, center.name, pairlist)
        #print(backangles.keys(), myangles.keys())
        for l,r in itertools.combinations(pairlist, 2):
            aname = (l.name, center.name, r.name)
            aname2 = (r.name, center.name, l.name)
            if aname in myangles: stats = myangles[aname]
            elif aname2 in myangles: stats = myangles[aname2]
            elif aname in backangles: stats = backangles[aname]
            elif aname2 in backangles: stats = backangles[aname2]
            else:
                n1, n2, n3 = aname
                raise KeyError("Could not find angle %s-%s-%s for residue %s" % (n1, n2, n3, self.resname))
            yield (l,center,r, stats['mean'], stats['std'])

# other effects:
# ASP, GLU have extra Hs at some pH values
# Histidine has proton exchange issues
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

# to_attach = {'GLU':[('H','HE2','OE2')], 'PRO':[('H','H','N')]}

# from richards
RichardsSizes = {
    'CH3' : 2.0, 'CH2' : 2.0, 'CH' : 2.0, 'C' : 1.7, 'O' : 1.4,
    'OH' : 1.6, 'NH3' : 2.0, 'NH2' : 1.75, 'NH' : 1.7, 'N' : 1.7, # 'N' not in Richards
     'S' : 1.8
    }

AliceSizes = {'H':1.05, 'Csp3':1.5, 'Csp2':1.4, 'N':1.4, 'O':1.45,
                    'S' : 1.6}
PorterRoseSizedict = {'H':1.0, 'Csp3':1.64, 'Csp2':1.5, 'N':1.35, 'O':1.35,
                    'S' : 1.8} # Sulfur taken from Richards

# following RasMol
#~ bonddists = {
    #~ 'CC':2,'CN':1.96,'CO':1.96,'PP':2.632,'OP':2.276,'SS':2.6,'OS':2.26,
    #~ 'CaO':2.232,'SZn':3.028
    #~ }

# Charges at Ph 3, from Abhi
Charges3 = { ('LYS', 'NZ') : 1,
    ('ARG', 'NH1') : .39, ('ARG', 'NH2') : .39, ('ARG', 'NE') : .22,
    ('HIS', 'ND1') : .5, ('HIS', 'NE2') : .5,
    ('ASP', 'OD1') : -.025, ('ASP', 'OD2') : -.025, 
    ('GLU', 'OE1') : -.025, ('GLU', 'OE2') : -.025
    }

# Charges at Ph 7.4, from Abhi
Charges74 = { ('LYS', 'NZ') : 1,
    ('ARG', 'NH1') : .39, ('ARG', 'NH2') : .39, ('ARG', 'NE') : .22,
    ('HIS', 'ND1') : .05, ('HIS', 'NE2') : .05,
    ('ASP', 'OD1') : -.5, ('ASP', 'OD2') : -.5, 
    ('GLU', 'OE1') : -.5, ('GLU', 'OE2') : -.5
    }

# charges per residue
rescharges74 = dict()
for r,a in Charges74: rescharges74[r] = rescharges74.get(r,0) + Charges74[(r,a)]

class Resvec(atomvec):
    def __init__(self, residue, H=False, amu=1, loadfile=None):
        residue.sort()
        if not H:
            self.atoms = [a for a in residue.child_list if a.name[0] != 'H']
            self.hydrogens = [self._getHs(residue, a) for a in self.atoms]
            #~ self.hydrogens = None
        else:
            self.atoms = residue.child_list
            self.hydrogens = None
        self.masses = [self._get_mass(residue, a, H)*amu for a in self.atoms]
            
        atomvec.__init__(self, self.masses)
        self.set_locs(residue)
        self.loadfile = loadfile
        self.resname = residue.resname
        
    # dictionaries with means, to be loaded as necessary
    resbonds = None
    backbonds = None
    resangles = None
    backangles = None
    
    def set_locs(self, residue):
        atoms = self.atoms
        if(len(atoms) != len(self)):
            raise TypeError("Resvec.setLocs: lengths do not match")
        #locs = [a.get_coord().tolist() for a in atoms]
        locs = [tuple(map(float, a.get_coord())) for a in atoms]
        for loc, atom in zip(locs, self):
            atom.x = Vec(*loc)
    
    def _getHs(self, res, atom):
        """res is the Bio.PDB Residue class"""
        if atom.name == 'C' or atom.name == 'O':
            return []
        if atom.name == 'N':
            return [a for a in res.child_list if a.name == 'H']
        name = 'H' + atom.name[1:]
        names = (name, name + '1', name + '2', name + '3', name + '4')
        return [a for a in res.child_list if a.name in names]
    
    def getHs(self, atom):
        if atom.name == 'C' or atom.name == 'O':
            return []
        if atom.name == 'N':
            return [a for a in self.atoms if a.name == 'H']
        name = 'H' + atom.name[1:]
        names = (name, name + '1', name + '2', name + '3', name + '4')
        return [a for a in self.atoms if a.name in names]
    
    def _get_mass(self, res, atom, H = True):
        """H means treat Hydrogens as separate"""
        if H:
            return atom.mass
        return sum(h.mass for h in self._getHs(res, atom)) + atom.mass
    
    def _formula(self, indx):
        atom = self.atoms[indx]
        if not self.hydrogens:
            return atom.element
        Hs = len(self.hydrogens[indx])
        if Hs == 0:
            return atom.element
        elif Hs == 1:
            return atom.element + 'H'
        return atom.element + 'H' + str(Hs)
    
    def __contains__(self, atom):
        if isinstance(atom, str):
            return any(a.name == atom for a in self.atoms)
        return atom in self.atoms
    
    @classmethod
    def load_full(cls, f=defaultloadfile):
        if None not in (cls.resbonds, cls.backbonds, cls.resangles, cls.backangles):
            return
        opened = False
        if not (hasattr(f, 'read') and hasattr(f, 'readline')):
            f = open(f, 'rb')
            opened = True
        
        import pickle
        cls.resbonds, cls.backbonds, cls.resangles, cls.backangles = pickle.load(f)
        if opened: f.close()
    
    def load_data(self, f=None):
        if None not in (self.resbonds, self.backbonds, self.resangles, self.backangles):
            return
        if f is None:
            f = self.loadfile
            if self.loadfile is None:
                raise KeyError("No file or filename provided")
        self.load_full(f)
    
    def getCOM(self):
        vs,ms = list(zip(*[(a.x*a.mass, a.mass) for a in self]))
        return sum(vs, Vec(0,0,0)) / sum(ms)
    
    def getCORg(self):
        vs = [a.x for a in self]
        return sum(vs, Vec(0,0,0)) / len(vs)
    
    def get_orbital(self, atom):
        if atom.element != 'C':
            raise NotImplementedError
        self.load_data()
        rbonds = self.resbonds[self.resname]
        mybonds = ([b for b in rbonds if atom.name in b] +
                    [b for b in self.backbonds if atom.name in b])
        #~ mybonds += [(atom.name, H.name) for H in self.getHs(atom)]
        #~ if 2 <= len(mybonds) <= 3: return 'sp2'
        if len(mybonds) == 3: return 'sp2'
        elif len(mybonds) == 4: return 'sp3'
        else:
            raise NotImplementedError(
                "Can't deal with atom %s (%s) with %d bonds" % (
                    atom.name, atom.element, len(mybonds)), mybonds)
    
    def get_angles(self, lastres=None, nextres=None, f=None):
        """Yields quintuplets of (atom1 ptr, atom2 ptr, atom3 ptr, bond angle in radians, std in radians)"""
        self.load_data(f)
        myangles = self.resangles[self.resname]
        backangles = self.backangles
        dlist = defaultdict(list)
        for a1,a2,d,s in self.get_bonds(lastres, nextres, f):
            dlist[(a1.group, a1.name)].append(a2)
            dlist[(a2.group, a2.name)].append(a1)
        for center, pairlist in dlist.items():
            cgroup, cname = center
            center = cgroup[cname]
            #print('get_angles:',self.resname, center.name, pairlist)
            #print(backangles.keys(), myangles.keys())
            for l,r in itertools.combinations(pairlist, 2):
                aname = (l.name, center.name, r.name)
                aname2 = (r.name, center.name, l.name)
                if aname in myangles: stats = myangles[aname]
                elif aname2 in myangles: stats = myangles[aname2]
                elif aname in backangles: stats = backangles[aname]
                elif aname2 in backangles: stats = backangles[aname2]
                else:
                    n1, n2, n3 = aname
                    raise KeyError("Could not find angle %s-%s-%s for residue %s" % (n1, n2, n3, self.resname))
                yield (l,center,r, stats['mean'], stats['std'])
        
    def old_get_angles(self, lastres=None, nextres=None, f=None):
        """Yields triplets of (atom1 ptr, atom2 ptr, atom3 ptr, bond angle in radians)"""
        self.load_data(f)
        myangles = self.resangles[self.resname]
        for a1, a2, a3 in myangles:
            if (any(a.startswith('H') for a in (a1,a2,a3)) 
                    and self.hydrogens is not None):
                # we have NO hydrogens
                continue
            if (any(a.startswith('H') for a in (a1,a2,a3)) 
                    and self.hydrogens is None
                    and any(a not in self for a in (a1,a2,a3))):
                        # angle with H in it, and we have H1, H2, H3
                continue
            if 'OXT' in (a1,a2,a3):
                try:
                    yield (self[a1], self[a2],self[a3], myangles[a1, a2, a3]['mean'])
                except KeyError:
                    continue
            yield (self[a1], self[a2],self[a3], myangles[a1, a2, a3]['mean'])
        
        for trip, val in self.backangles.items():
            a1, a2, a3 = trip
            if (any([a in trip for a in ['H1','H2','H3','OXT']]) 
                    or ('H' in trip and self.hydrogens is not None)
                    or ('H' in trip and any([a in self for a in ['H1','H2','H3']]))):
                try:
                    a1 = self[a1]
                    a2 = self[a2]
                    a3 = self[a3]
                    yield (a1, a2, a3, val['mean'])
                except KeyError:
                    continue
            elif a2 == 'C' and 'N' in trip:
                #this is the neighbor-neighbor peptide bond
                if nextres:
                    yield (self[a1], self['C'], nextres[a3], val['mean'])
            elif a2 == 'N' and 'C' in trip:
                #this is the neighbor-neighbor peptide bond
                if lastres:
                    yield (lastres[a1], self['N'], self[a3], val['mean'])
            else:
                yield (self[a1], self[a2], self[a3], val['mean'])
    
    def get_bonds(self, lastres=None, nextres=None, f=None):
        """Yields quadruplets of (atom1 ptr, atom2 ptr, dist, std)"""
        self.load_data(f)
        bonds = bonddict[self.resname] + backbonds
        if self.hydrogens is None:
            #~ print('Extending...')
            bonds.extend(hydrodict[self.resname])
        if lastres is None: bonds = add_first_bonds(bonds, H=(self.hydrogens is None))
        else: bonds.extend(pairbonds)
        if nextres is None: bonds = add_last_bonds(bonds, H=(self.hydrogens is None))
        
        mybonds = self.resbonds[self.resname]
        for pair in bonds:
            a1,a2 = pair
            stats = (mybonds[pair] 
                    if pair in mybonds else self.backbonds[pair])
            dist = stats['mean']
            std = stats['std']
            if pair in pairbonds:
                yield (lastres[a1], self[a2], dist, std)
                continue
            yield (self[a1], self[a2], dist, std)
        
    def old_get_bonds(self, lastres=None, f=None):
        """Yields triplets of (atom1 ptr, atom2 ptr, bond length in angstrom)"""
        self.load_data(f)
        mybonds = self.resbonds[self.resname]
        madebonds = set()
        for a1, a2 in mybonds:
            if (((a1.startswith('H') and a1 not in self)
                 or (a2.startswith('H') and a2 not in self))
                    and self.hydrogens is not None):
                        # angle with H in it, and we don't have Hs
                        #~ print('ignoring', a1, a2)
                        continue
            madebonds.add(a1)
            madebonds.add(a2)
            yield (self[a1], self[a2], mybonds[a1, a2]['mean'])
        for pair, val in self.backbonds.items():
            a1, a2 = pair
            if 'OXT' in pair:
                if 'OXT' not in self: continue
                a1 = self[a1]
                a2 = self[a2]
                madebonds.add(a1.name)
                madebonds.add(a2.name)
                yield (a1, a2, val['mean'])   
                madebonds.add(a1)
                madebonds.add(a2)
            elif 'C' in pair and 'N' in pair:
                #this is the neighbor-neighbor peptide bond
                if lastres:
                    yield (lastres['C'], self['N'], val['mean'])
                    madebonds.add('N')
            elif 'H' in pair and 'H' not in self and self.hydrogens is not None:
                continue
            elif 'H' in pair  and 'H' not in self and 'H1' in self:
                continue
            elif a1 in ('H1', 'H2', 'H3') and a1 not in self:
                continue
            elif a2 in ('H1', 'H2', 'H3') and a2 not in self:
                continue
            else:
                try:
                    yield (self[a1], self[a2], val['mean'])
                    madebonds.add(a1)
                    madebonds.add(a2)
                except KeyError:
                    atm = a1 if a1 not in self else a2
                    raise KeyError("Could not find atom %s in this %s residue %s" %
                     (atm, self.resname, [a.name for a in self] ))
        allnames = set([a.name for a in self])
        missed = allnames.difference(madebonds)
        if len(missed) > 0:
            raise RuntimeError('Missing bonds for %s; it is not attached to anything.' % sorted(missed))
    
    def __getitem__(self, obj):
        indx = obj
        if not isinstance(obj, int):
            try:
                atom, = [a for a in self.atoms if a.name == obj]
            except ValueError:
                raise KeyError(str(obj) + " not found, or multiple found (" +
                    self.resname + ") " + str([a.name for a in self.atoms]))
            indx = self.atoms.index(atom)
        
        a = atomvec.__getitem__(self, indx)
        a.name = self.atoms[indx].name
        a.mass = self.masses[indx]
        a.pdbatom = self.atoms[indx]
        a.element = a.pdbatom.element
        if self.hydrogens:
            a.hydrogens = self.hydrogens[indx]
            a.formula = self._formula(indx)
        else:
            a.formula = a.element
        a.group = self
        return a
    
    def get_name(self, atom):
        for avec_atom, res_atom in zip(self, self.atoms):
            #~ print(atom, avec_atom)
            if atom.this == avec_atom.this:
                return res_atom.name
    
    @classmethod
    def all_bonds(cls, reslist):
        """For a list of residues, yields (atom1, atom2, bond length) triplets.
        
        The bond length comes from the statistics file."""
        reslist = list(reslist)
        last_residues = [None] + reslist
        next_residues = reslist[1:] + [None]
        lists_of_bonds = (res.get_bonds(lres,nres) for res, lres, nres
                                        in zip(reslist, last_residues, next_residues))
        
        return itertools.chain.from_iterable(lists_of_bonds)
        
    @classmethod
    def stds(cls, reslist):
        """returns the standard deviation"""
        diffs = np.array([(a1.x - a2.x).mag() - blength
            for a1, a2, blength in cls.all_bonds(reslist)])
        return diffs.std()
    
    @classmethod
    def all_angles(cls, reslist):
        reslist = list(reslist)
        last_residues = [None] + reslist
        next_residues = reslist[1:] + [None]
        all_angs = (res.get_angles(lres,nres) for res, lres, nres
                                        in zip(reslist, last_residues, next_residues))
        
        return itertools.chain.from_iterable(all_angs)
    
    @classmethod
    def to_xyz(cls, reslist, f):
        """Print out coordinates into a file"""
        reslist = list(reslist)
        print(sum(len(r) for r in reslist), file=f)
        for res in reslist:
            for a in res:
                print("%s %.3f %.3f %.3f" % tuple(a.name[:1], *(a.x)), file=f)
    
    @classmethod
    def from_pdb(cls, strucname='aS', pdbfilename=defaultpdbfile, 
                        loadfile=defaultloadfile, fix=True,
                        numchains=None, numres=None, amu=1, H=False):
        pdbp = PDBParser(QUIET=True)
        aS = pdbp.get_structure(strucname, pdbfilename)
        chains = list(aS.get_chains())[:numchains]
        residues = list(itertools.chain.from_iterable(
                            chain.child_list for chain in chains))[:numres]
        if fix:
            #~ print('fixing...')
            for r in residues:
                bs = bonddict[r.resname] + backbonds + hydrodict[r.resname]
                if r is residues[0]: bs = add_first_bonds(bs)
                if r is residues[-1]: bs = add_last_bonds(bs)
                if not H:
                    bs = [(a1, a2) for (a1, a2) in bs if 'H' not in (a1[0], a2[0])]
                added, removed = fixResidue(r, bs)
                if not H:
                    removed = [a for a in removed if a[0] != 'H']
                if added or removed:
                    print('Fixed',r.resname, 
                            'added:', ', '.join(added),
                            'removed:', ', '.join(removed))
        
        rvecs = [Resvec(res, amu=amu, H=H, loadfile=loadfile) for res in residues]
        return rvecs
    
    def dihedral(self, ang, prev, nxt):
        if ang=='omega': angs, sgn = ([], ['CA', 'C'],['N','CA']), -1
        elif ang=='psi': angs, sgn = ([], ['N','CA','C'],['N']), 1
        elif ang=='phi': angs, sgn = (['C'],['N','CA','C'], []), 1
        else: raise ValueError("Angle %r not recognized." % ang)
        
        atoms = [r[a] 
                    for r,anames in zip([prev, self, nxt], angs)
                    for a in anames]
        assert len(atoms) == 4
        locs = [a.x for a in atoms]
        r1,r2,r3 = [x2-x1 for x1,x2 in zip(locs[:-1], locs[1:])]
        anames = "-".join(["%s%s:%s" % 
                (r.resname if r else '#', 
                            str([prev,self,nxt].index(r)) if r else '', a)
                    for r,anames in zip([prev, self, nxt], angs)
                    for a in anames])
        
        lmin = 3
        assert r1.mag() < lmin, (
            "r1 mag was %.3f for angle \'%s\', residues %s,%s,%s" % (
            r1.mag(), ang, prev.resname if prev else '-', 
            self.resname if self else '-', nxt.resname if nxt else '-') + " " + anames)
        assert r2.mag() < lmin, (
            "r2 mag was %.3f for angle \'%s\', residues %s,%s,%s" % (
            r2.mag(), ang, prev.resname if prev else '-', 
            self.resname if self else '-', nxt.resname if nxt else '-') + " " + anames)
        assert r3.mag() < lmin, (
            "r3 mag was %.3f for angle \'%s\', residues %s,%s,%s" % (
            r3.mag(), ang, prev.resname if prev else '-', 
            self.resname if self else '-', nxt.resname if nxt else '-') + " " + anames)
    
        return sgn * dihedral.getang(r1,r2,r3)
        #cos = dihedral.getcos(r1,r2,r3)
        #return math.acos(cos)
    
    @classmethod
    def phipsis(cls, rlist):
        """Note that psis are for 0:N-1, and phis for 1:N. You need
        to do phis[:-1], psis[1:] to get matching pairs..."""
        restriplets = list(zip([None] + rlist, rlist, rlist[1:]))
        psis = [res.dihedral('psi', p, n) for p,res,n in restriplets]
        restriplets = list(zip(rlist, rlist[1:], rlist[2:] + [None]))
        phis = [res.dihedral('phi', p, n) for p,res,n in restriplets]
        assert len(phis) == len(rlist) - 1
        assert len(psis) == len(rlist) - 1
        return phis,psis

class Residue(_Residue):
    """I don't think I'm using this anymore."""
    def __init__(self, res):
        self.__dict__.update(res.__dict__)
        
    def atoms(self, H=False):
        self.sort()
        if not H:
            return [a for a in self.child_list if a.name[0] != 'H']
        return self.child_list
    
    def masses(self, H=False, amu=1):
        if H:
            return [a.mass/amu for a in self.atoms(H)]
        return [sum(h.mass for h in self.get_Hs(a)) + a.mass 
                                for a in self.atoms(H)]
    
    def set_locs(self, avec, H=False):
        atoms = self.atoms(H)
        if(len(atoms) != len(avec)):
            raise TypeError("Residue.setLocs: lengths do not match")
        locs = [a.get_coord().tolist() for a in atoms]
        for loc, atom in zip(locs, avec):
            atom.x = Vec(*loc)
    
    def make_atomvec(self, H=False, amu=1):
        avec = atomvec([m/amu for m in self.masses(H,amu)])
        self.set_locs(avec, H)
        return avec
    
    def get_Hs(self, atom):
        name = atom.name[:]
        names = (name, name + '1', name + '2', name + '3', name + '4')
        return [a for a in self.child_list if a.name in names]
        
    
    def get_bonds(self, avec):
        """Returns a list of tuples (atom1, atom2, length) indicating bonds.
        
        Note: does not include hydrogens."""
        atoms = self.atoms(False)
        lst = bonddict[self.resname]
        for a in atoms:
            aindx = atoms.index(a)
            for b in atoms[aindx+1:]:
                bindx = atoms.index(b)
                if (a.name, b.name) in lst:
                    yield (avec.get(aindx), avec.get(bindx), 1)
    
    def old_get_bonds(self, avec):
        """Returns a list of tuples (atom1, atom2) indicating bonds.
        
        Note: does not include hydrogens.
        This version just looks for minimum distances."""
        atoms = self.atoms(False)
        for a in atoms:
            aindx = atoms.index(a)
            for b in atoms[aindx+1:]:
                bindx = atoms.index(b)
                aloc = Vec(*a.get_coord().tolist())
                bloc = Vec(*b.get_coord().tolist())
                dist = (aloc - bloc).mag()
                if dist > 3:
                    continue
                mindist = bonddists["".join(sorted(a.element + b.element))]
                if dist <= mindist:
                    print(a.name, b.name)
                    yield (avec.get(aindx), avec.get(bindx), a, b)

def addAtom(residue, attachto, newname, element='H', coord=None, **kwargs):
    """This takes a Bio.PDB.Residue, the atom in that residue to attach
    to, and the name of the new atom to make, and it makes a new atom."""
    rname = residue.resname
    if isinstance(attachto, str): attachto = residue[attachto]
    if coord is None:
        for a in residue:
            print(attachto, a)
            print(attachto - a)
        bonds = [a for a in residue if attachto - a < 4 and a is not attachto]
        #~ def getpaired(bond):
            #~ a1, a2 = bond
            #~ return a2 if a1 == attachto.name else a1
        #~ # get just the neighbors, not the tuple pairs
        #~ bonds = [getpaired(b) for b in bonds]
        #~ # get the atoms that exist
        #~ bonds = [residue[a] for a in bonds if a in residue.child_dict]
        # vectors from bonded atoms to attachto location
        attachloc = attachto.get_vector()
        dvecs = [a.get_vector() - attachloc for a in bonds]
        # average them to get a new location for the next atom
        # this will hopefully put it as far away from the atoms attached to the
        # 'attachto' atom as possible
        coordoffset =  avg(dvecs)
        # make sure our newloc has length 2. 2 angstroms sounds good; the simulation can handle that.
        coord = attachloc + (coordoffset.normalized() ** 2)
    
    # make our atom!
    # these are my 'default' arguments for Bio.PDB.Atom
    kw = dict(coord = coord, bfactor=0,
            occupancy=0, altloc = 0, fullname=newname,
            serial_number=0, element=element)
    # update with anything given
    kw.update(kwargs)
    newatom = _Atom(newname, **kw)
    # add it to the residue, which will also give it everything we need
    residue.add(newatom)

def fixResidue(res, bonds):
    """Adds and deletes atoms from a residue to match what is in bonds.
    
    Returns a tuple pair (added, removed)"""
    reshas = set([a.name for a in res])
    atoms = set([a for p in bonds for a in p])
    for a in reshas.difference(atoms):
        del res[a]
    for atm in atoms.difference(reshas):
        bs = [a for p in bonds if atm in p for a in p if a != atm]
        addAtom(res, bs[0], atm, atm[0])
    return (atoms.difference(reshas), reshas.difference(atoms))

def addAtoms(residues, atomdict):
    for res, lst in atomdict.items():
        for element,aname, attach in lst:
            allres = [r for r in residues if r.resname == res and aname not in r]
            #~ print('%d residues of type %s' % (len(allres), res))
            for r in allres:
                addAtom(r, attach, aname, element=element)

def make_bonds(resvecs, bond_k, angstrom=1, usestd=False):
    bond_pairs = bondpairs()
    allks = []
    for a1,a2,l,s in Resvec.all_bonds(resvecs):
        k = bond_k / s /s if usestd else bond_k
        allks.append(k)
        length = l * angstrom
        bond_pairs.add(k, length, a1, a2)
    print('Made %d bonds from %d residues, with k=%.2f +- %.2f [%.2f - %.2f]' % (
        bond_pairs.size(), len(resvecs), 
        np.mean(allks), np.std(allks), np.min(allks), np.max(allks)))
    return bond_pairs

def make_angles(resvecs, angle_k, usestd=False):
    bond_angles = angletriples()
    allks = []
    for a1,a2,a3,angle,s in Resvec.all_angles(resvecs):
        #~ print("angle", "-".join(a.name for a in (a1,a2,a3)), angle)
        k = angle_k / s / s if usestd else angle_k
        allks.append(k)
        bond_angles.add(k, angle, a1, a2, a3)
    print('Made %d angles from %d residues, with k=%.2f +- %.2f [%.2f - %.2f]' % (
            bond_angles.size(), len(resvecs), 
        np.mean(allks), np.std(allks), np.min(allks), np.max(allks)))

    return bond_angles

def make_dihedrals(resvecs, k):
    d = dihedrals()
    k = float(k)
    piang = fvector([k/2,k,k/2])
    zeroang = fvector([k/2,-k,k/2])
    for r1,r2 in zip(resvecs, resvecs[1:]):
        #~ ang = piang if r1.resname is 'PRO' or r2.resname is 'PRO' else zeroang
        # ang = zeroang # corresponds to planar zigzag, not planar C <-- OLD!!
        ang = piang
        a1, a2, a3, a4 = r1['CA'],r1['C'],r2['N'],r2['CA']
        r1,r2,r3 = (a2.x-a1.x), (a3.x-a2.x), (a4.x-a3.x)
        mags = r1.mag(), r2.mag(), r3.mag()
        #~ if ang == piang:
            #~ print('added piang', mags)
        #~ if ang == zeroang:
            #~ print('added zeroang', mags)
        for m in mags:
            if m >= 2.5: raise RuntimeError("Found bad angle")
        
        curang = dihedral.getang(r1,r2,r3)
        #print('angle:', curang)
        #if abs(curang) < 2.8:
        #    raise RuntimeError, "Angle %.2f is not near pi" % curang
        d.add(ang, a1, a2, a3, a4)
    return d

def make_LJ(resvecs, epsilon=1, sizefactor=1, neighborcutoff=1.4, LJdict=AliceSizes):
    """Neighborcutoff relative to sigma.
    
    returns (LJ, neighbors)"""
    factor = 2 * sizefactor * (2**(1.0/6.0))
    # 2 because LJ-sigma is the sum of two radii, (2**(1.0/6.0)) for the shape
    
    maxsigma = max(LJdict.values()) * factor
    innerradius = maxsigma
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(innerradius, outerradius)
    LJ = LJgroup(neighbors)
    
    keygens = (lambda a,r: a.formula, lambda a,r: atom.element, 
                lambda a,r: a.element + r.get_orbital(atom))
    def getsize(atom, res):
        for keyfunc in keygens:
            key = keyfunc(atom, res)
            if key in LJdict: return LJdict[key]
        #~ print(*((atom.name, atom.hydrogens, atom.formula) for atom in atom.group))
        raise KeyError('Atom ' + atom.name + ' (' + atom.formula
                + ') in ' + rname + ',' + str(resvecs.index(atom.group))
                + ' not found in LJdict.')
    
    for atom,res in ((atom,r) for r in resvecs for atom in r):
        aref = res.get_id(atom)
        core = getsize(atom, res)
        sigma = core * factor
        LJ.add(LJatom(epsilon, sigma, aref))
    
    for a1,a2,l,bondstd in Resvec.all_bonds(resvecs):
        neighbors.ignore(a1,a2)
    
    for a1,a2,a3,angle,anglestd in Resvec.all_angles(resvecs):
        #~ neighbors.ignore(a1,a2)
        #~ neighbors.ignore(a2,a3)
        neighbors.ignore(a1,a3)
    
    return LJ, neighbors

def make_LJ_simple(resvecs, LJcutoff=2.5, epsilon=1):
    """LJcutoff in sigma units. Useful only for testing the neighbor list."""
    LJ = LJsimple(LJcutoff)
    
    fullgroup = metagroup(resvecs)
    
    for i, atom in zip(range(fullgroup.size()), (atom for r in resvecs for atom in r)):
        assert fullgroup.get(i) == atom
        try:
            sigma = LJdict[atom.formula]
        except KeyError:
            print(*((atom.name, atom.hydrogens, atom.formula) for atom in atom.group))
            raise KeyError('Atom ' + atom.name + ' (' + atom.formula + ') in ' 
                    + atom.group.resname + ',' + str(resvecs.index(atom.group))
                    + ' not found in LJdict.')
        LJ.add(fullgroup.get_id(i), epsilon, sigma)
    
    for a1,a2,l in Resvec.all_bonds(resvecs):
        LJ.ignore(a1,a2)
    
    for a1,a2,a3,angle in Resvec.all_angles(resvecs):
        #~ LJ.ignore(a1,a2)
        #~ LJ.ignore(a2,a3)
        LJ.ignore(a1,a3)
    
    return LJ

def make_charges(resvecs, screen, k, ph=7.4, expectOXT=True, subtract=True):
    """Screening distance in base units - 0 or less for no screening.
    
    Epsilon for an extra charge."""
    charges = Charges(screen, k)
    if ph == 7.4:
        chargedict = Charges74
    elif ph == 3:
        chargedict = Charges3
    else:
        raise NotImplementedError
    
    fullgroup = metagroup(resvecs)
    atoms = [(atom, chargedict[(r.resname, atom.name)]) for r in resvecs for atom in r
                        if (r.resname, atom.name) in chargedict]
    # Add extra pieces
    # N-terminus
    #~ aref = fullgroup.get_id(resvecs[0]['N'])
    #~ charges.add(aref, 1)
    atoms.append((resvecs[0]['N'], 1))
    # C-terminus
    if expectOXT or 'OXT' in resvecs[-1]:
        atoms.append((resvecs[-1]['OXT'], -.5))
        atoms.append((resvecs[-1]['O'], -.5))
    else:
        # no OXT, put all on the 'O'
        #~ aref = fullgroup.get_id(resvecs[-1]['O'])
        #~ charges.add(aref, -1)
        atoms.append((resvecs[-1]['O'], -1))
    
    totcharge = sum([float(q) for a,q in atoms])
    tosubtract = totcharge / len(atoms)
    if subtract:
        atoms = [(a,q-tosubtract) for a,q in atoms]
        #~ oldtot = totcharge
        #~ print("new sum:", sum([float(q) for a,q in atoms]))
    
    if subtract:
        print(len(atoms), "charged atoms, total charge was %.2f" % totcharge)
    else: print(len(atoms), "charged atoms, total charge %.2f" % totcharge)
    
    #~ residues = set(a.group for a,q in atoms)
    #~ import collections
    #~ d = collections.defaultdict(int)
    #~ for r in residues: d[r.resname] += 1
    #~ for r, n in sorted(d.items()):
        #~ print(n, r)
    
    for atom,charge in atoms:
        #~ print('ref for', atom)
        aref = fullgroup.get_id(atom)
        #~ charge = chargedict[(atom.group.resname, atom.name)]
        charges.add(aref, charge)
    
    for a1,a2,l,bondstd in Resvec.all_bonds(resvecs):
        key1, key2 = (a1.group.resname, a1.name), (a2.group.resname, a2.name)
        if key1 not in chargedict or key2 not in chargedict:
            continue
        #~ print('ignoring', a1.group.resname, a1.name, a2.name)
        charges.ignore(a1,a2)
    
    for a1,a2,a3,angle,anglestd in Resvec.all_angles(resvecs):
        key1, key2 = (a1.group.resname, a1.name), (a3.group.resname, a3.name)
        if key1 not in chargedict: continue
        if key2 not in chargedict: continue
        charges.ignore(a1,a3)
    
    return charges

def make_hydrophobicity_old(resvecs, epsilons, LJcutoff=2.5, neighborcutoff=1.4, sigma=5.38):
    """LJcutoff in sigma units, neighborcutoff relative to LJcutoff.
    Makes a CA hydrophobicity attractive force.
    
    returns (interaction, neighbors)"""
    
    # 2 because LJ-sigma is the sum of two radii, (2**(1.0/6.0)) for the shape
    
    innerradius = sigma*LJcutoff
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(innerradius, outerradius)
    Hphob = Hydrophobicity(neighbors)
    
    def makeatom(res):
        atom = res['CA']
        Hindex = hydroindex[res.resname][1]
        index = 0 if Hindex < 0 else 1
        epsvec = epsilons[index]
        return HydroAtom(epsvec, index, sigma, res.get_id(atom), LJcutoff)
    
    atoms = [makeatom(r) for r in resvecs]
    for atom in atoms:
        Hphob.add(atom)
    for a1, a2 in zip(atoms, atoms[1:]):
        neighbors.ignore(a1,a2)
    
    return Hphob, neighbors

def make_hydrophobicity(resvecs, epsilon, LJcutoff=2.5, neighborcutoff=1.4,
                    sigma=5.38, hydroindex=hydroindex7, method='geometric'):
    """LJcutoff in sigma units, neighborcutoff relative to LJcutoff.
    Makes a CA hydrophobicity attractive force.
    
    returns (interaction, neighbors)"""
    
    # 2 because LJ-sigma is the sum of two radii, (2**(1.0/6.0)) for the shape
    
    innerradius = sigma*LJcutoff
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(innerradius, outerradius)
    Hphob = Hydrophobicity(neighbors)
    
    #~ print("Hlist:", *Hlist)
    Rnames, Hlist = list(zip(*sorted(hydroindex.items())))
    
    if method == 'geometric':
        assert all(0 <= H <= 1 for H in Hlist)
        def combine(H1, H2): return math.sqrt(H1*H2)
    elif method == 'geometriczero':
        def combine(H1, H2): return math.sqrt(H1*H2) if (H1 > 0 and H2 > 0) else 0
    elif method == 'arithmetic':
        def combine(H1, H2): return (H1 + H2)/2
    elif method == 'arithmeticzero':
        def combine(H1, H2): return (H1 + H2)/2 if (H1 + H2)/2 > 0 else 0
    else:
        raise NotImplementedError("Method %s not recognized" % method)
    
    printedres = set()
    def makeatom(res):
        atom = res['CA']
        #~ if res.resname not in Hlist: return None
        idx = Rnames.index(res.resname)
        myH = Hlist[idx]
        epsvec = [epsilon*combine(H,myH) for H in Hlist]
        #~ print(res.resname, ' '.join([str(int(H/epsilon*100)) for H in epsvec]))
        #~ if res.resname not in printedres:
            #~ print(res.resname, epsilon, '%.3f' % epsvec[idx])
            #~ printedres.add(res.resname)
        return HydroAtom(epsvec, idx, sigma, res.get_id(atom), LJcutoff)
    
    atoms = [makeatom(r) for r in resvecs]
    atoms = [a for a in atoms if a is not None]
    
    for atom in atoms:
        Hphob.add(atom)
    for a1, a2 in zip(atoms, atoms[1:]):
        neighbors.ignore(a1,a2)
    
    return Hphob, neighbors
    
def correct_phi(residues):
    restrip = [(prev, res, nxt) for prev, res, nxt in 
            zip(residues[:-2], residues[1:-1], residues[2:]) 
            if res.resname not in ('PRO', 'GLY')]
    phis = [ res.dihedral('phi', prev, nxt) for prev, res, nxt in restrip]
    phiscorrect = [(phi > 2*math.pi/3) or (phi < 0) for phi in phis]
    return phiscorrect.count(True), len(restrip)

########################################################################
def printlist(lst, name, **kw):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    stdstr = '%8.3f' % std if std < 10000 else '%8.5g' % std
    sovm = 100*std / mean if mean > 0 else float('nan')
    last = float(lst[-1])
    lstd = 100*(last - mean) / std if std > 0 else float('nan')
    print("{name:10s}={mean:10.5g}, σ={std:10.4g} ({sovm:9.3g}%) [{last:9.2f} ({lstd:7.3g}%)]".format(**locals()), **kw)

def to_dhms(tdelta):
    days = tdelta.days
    hr, sec = divmod(tdelta.seconds, 3600)
    mn, sec = divmod(sec, 60)
    return (days, hr, mn, sec)

def get_eta(cur, tot, starttime):
    from datetime import datetime, timedelta
    now = datetime.now()
    tused = now - starttime
    
    # timedeltas cannot be multiplied by floats (?), so we get out the 
    #seconds, use that, and make it back into a dt
    usedsecs = 24*60*60 * tused.days + tused.seconds
    tleft = timedelta(0, (usedsecs) * ((tot - cur) / cur))
    
    endtime = now + tleft
    return endtime, to_dhms(tleft)

def interval_str(days, hr, mn, sec):
    return (('%dd ' % days) if days else '') + ('%d:%02d:%02d'% (hr, mn, sec))
