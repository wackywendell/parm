from os.path import expanduser
mydir = expanduser('~/idp/')
defaultpdbfile  = mydir + 'pdb/aS.pdb'
defaultloadfile = mydir + 'blengths/stats-H.pkl'

res31 = dict(ALA='A', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G',
    HIS='H', ILE='I', LYS='K', LEU='L', MET='M', ASN='N', PRO='P',
    GLN='Q', ARG='R', SER='S', THR='T', VAL='V', TRP='W', TYR='Y')
res13 = {v:k for k,v in res31.items()}

masses = dict(C=12.0107, H=1.00794, N=14.0067, O=15.9994, S=32.065)

# adapted from IUPAC standards, available at
# http://www.chem.qmul.ac.uk/iupac/misc/ppep4.html#400
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


# Original Monera scale
monera = moneraneg = {k : (v * 2 - 1) for k,v in hydroindex7.items()}
monerascale = {k : (v - min(monera.values())) / 2 for k,v in hydroindex7.items()}
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

# the double-bonded oxygen (Od) needed to be redone; 
# Diego thinks 1.35, Alice 1.4, jury is still out
AliceSizes = {'H':1.05, 'Csp3':1.5, 'Csp2':1.4, 'N':1.4, 'O2':1.45, 'O1':1.35,
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
