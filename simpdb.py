from __future__ import print_function
from Bio.PDB import PDBParser
from Bio.PDB.Residue import Residue as _Residue
from simw import *
import itertools

# adapted from IUPAC standards, available at
# http://www.chem.qmul.ac.uk/iupac/misc/ppep4.html#400

from os.path import expanduser
mydir = expanduser('~/idp/')
defaultpdbfile  = mydir + 'pdb/aS.pdb'
defaultloadfile = mydir + 'blengths/stats.pkl'

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

# from richards
LJdict = {
    'CH3' : 2.0, 'CH2' : 2.0, 'CH' : 2.0, 'C' : 1.7, 'O' : 1.4,
    'OH' : 1.6, 'NH3' : 2.0, 'NH2' : 1.75, 'NH' : 1.7, 'N' : 1.7, # 'N' not in Richards
     'S' : 1.8
    }

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
        locs = [a.get_coord().tolist() for a in atoms]
        for loc, atom in zip(locs, self):
            atom.x = Vec(*loc)
    
    def _getHs(self, res, atom):
        if atom.name == 'C' or atom.name == 'O':
            return []
        name = 'H' + atom.name[1:]
        names = (name, name + '1', name + '2', name + '3', name + '4')
        return [a for a in res.child_list if a.name in names]
    
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
    
    @classmethod
    def load_full(cls, f=defaultloadfile):
        if None not in (cls.resbonds, cls.backbonds, cls.resangles, cls.backangles):
            return
        if not (hasattr(f, 'read') and hasattr(f, 'readline')):
            f = open(f, 'rb')
        
        import cPickle as pickle
        cls.resbonds, cls.backbonds, cls.resangles, cls.backangles = pickle.load(f)
    
    def load_data(self, f=None):
        if None not in (self.resbonds, self.backbonds, self.resangles, self.backangles):
            return
        if f is None:
            f = self.loadfile
            if self.loadfile is None:
                raise KeyError("No file or filename provided")
        self.load_full(f)
    
    def get_angles(self, lastres=None, nextres=None, f=None):
        """Yields triplets of (atom1 ptr, atom2 ptr, atom3 ptr, bond angle in radians)"""
        self.load_data(f)
        myangles = self.resangles[self.resname]
        for a1, a2, a3 in myangles:
            #~ print("bonds:",a1,a2)
            yield (self[a1], self[a2],self[a3], myangles[a1, a2, a3]['mean'])
        
        for trip, val in self.backangles.items():
            #~ print("backbonds:",*pair)
            a1, a2, a3 = trip
            if 'OXT' in trip:
                try:
                    a1 = self[a1]
                    a2 = self[a2]
                    a3 = self[a3]
                    #~ print(a1,a2,*pair)
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
    
    def get_bonds(self, lastres=None, f=None):
        """Yields triplets of (atom1 ptr, atom2 ptr, bond length in angstrom)"""
        self.load_data(f)
        mybonds = self.resbonds[self.resname]
        for a1, a2 in mybonds:
            #~ print("bonds:",a1,a2)
            yield (self[a1], self[a2], mybonds[a1, a2]['mean'])
        for pair, val in self.backbonds.items():
            #~ print("backbonds:",*pair)
            a1, a2 = pair
            if 'OXT' in pair:
                try:
                    a1 = self[a1]
                    a2 = self[a2]
                    #~ print(a1,a2,*pair)
                    yield (a1, a2, val['mean'])
                except KeyError:
                    continue
            elif 'C' in pair and 'N' in pair:
                #this is the neighbor-neighbor peptide bond
                if lastres:
                    yield (lastres['C'], self['N'], val['mean'])
            else:
                yield (self[a1], self[a2], val['mean'])
    
    def __getitem__(self, obj):
        indx = obj
        if not isinstance(obj, int):
            try:
                atom, = [a for a in self.atoms if a.name == obj]
            except ValueError:
                raise KeyError(str(obj) + " not found, or multiple found")
            indx = self.atoms.index(atom)
        
        a = atomvec.__getitem__(self, indx)
        a.name = self.atoms[indx].name
        a.mass = self.masses[indx]
        a.pdbatom = self.atoms[indx]
        a.element = a.pdbatom.element
        if self.hydrogens:
            a.hydrogens = self.hydrogens[indx]
            a.formula = self._formula(indx)
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
        lastvecs = [None] + reslist
        lists_of_bonds = (res.get_bonds(lvec) for res, lvec
                                        in zip(reslist, lastvecs))
        
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
        lists_of_bonds = (res.get_angles(lres,nres) for res, lres, nres
                                        in zip(reslist, last_residues, next_residues))
        
        return itertools.chain.from_iterable(lists_of_bonds)
    
    @classmethod
    def to_xyz(cls, reslist, f):
        reslist = list(reslist)
        print(sum(len(r) for r in reslist), file=f)
        for res in reslist:
            for a in res:
                print("%s %.3f %.3f %.3f" % tuple(a.name[:1], *(a.x)), file=f)
    
    @classmethod
    def from_pdb(cls, strucname='aS', pdbfilename=defaultpdbfile, 
                        loadfile=defaultloadfile, numchains=None, amu=1):
        pdbp = PDBParser()
        aS = pdbp.get_structure(strucname, pdbfilename)
        chains = list(aS.get_chains())[:numchains]
        residues = itertools.chain.from_iterable(chain.child_list for chain in chains)
        rvecs = [Resvec(res, amu=amu, loadfile=loadfile) for res in residues]
        return rvecs

class Residue(_Residue):
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

def make_bonds(resvecs, bond_k, angstrom=1):
    bond_pairs = bondpairs()
    for a1,a2,l in Resvec.all_bonds(resvecs):
        length = l * angstrom
        bond_pairs.add(bond_k, length, a1, a2)
    return bond_pairs

def make_angles(resvecs, angle_k):
    bond_angles = angletriples()
    for a1,a2,a3,angle in Resvec.all_angles(resvecs):
        bond_angles.add(angle_k, angle, a1, a2, a3)
    return bond_angles

def make_dihedrals(resvecs, k):
    d = dihedrals()
    k = float(k)
    piang = fvector([k/2,k,k/2])
    zeroang = fvector([k/2,-k,k/2])
    for r1,r2 in zip(resvecs, resvecs[1:]):
        #~ ang = piang if r1.resname is 'PRO' or r2.resname is 'PRO' else zeroang
        ang = zeroang # corresponds to planar zigzag, not planar C
        a1, a2, a3, a4 = r1['CA'],r1['C'],r2['N'],r2['CA']
        r1,r2,r3 = (a2.x-a1.x), (a3.x-a2.x), (a4.x-a3.x)
        mags = r1.mag(), r2.mag(), r3.mag()
        #~ if ang == piang:
            #~ print('added piang', mags)
        #~ if ang == zeroang:
            #~ print('added zeroang', mags)
        for m in mags:
            if m >= 2.1: raise RuntimeError, "Found bad angle"
        
        cosine = dihedral.getcos(r1,r2,r3)
        #~ print('cosine:', cosine)
        #~ if abs(cosine) < 3:
            #~ raise RuntimeError, "Angle %.2f is near pi" % cosine
        d.add(ang, a1, a2, a3, a4)
    return d

def make_LJ(resvecs, epsilon=1, neighborcutoff=1.4):
    """LJcutoff in sigma units, neighborcutoff relative to LJcutoff.
    
    returns (LJ, neighbors)"""
    maxsigma = max(LJdict.values())*2
    fullgroup = metagroup(resvecs)
    innerradius = maxsigma
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(fullgroup, innerradius, outerradius)
    LJ = LJgroup(neighbors)
    
    for atom in (atom for r in resvecs for atom in r):
        aref = fullgroup.get_id(atom)
        try:
            sigma = LJdict[atom.formula] * 2
        except KeyError:
            print(*((atom.name, atom.hydrogens, atom.formula) for atom in atom.group))
            raise KeyError('Atom ' + atom.name + ' (' + atom.formula + ') in ' 
                    + atom.group.resname + ',' + str(resvecs.index(atom.group))
                    + ' not found in LJdict.')
        LJ.add(LJatom(sigma, epsilon, aref))
    
    for a1,a2,l in Resvec.all_bonds(resvecs):
        neighbors.ignore(a1,a2)
    
    for a1,a2,a3,angle in Resvec.all_angles(resvecs):
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
        LJ.add(fullgroup.get_id(i), sigma, epsilon)
    
    for a1,a2,l in Resvec.all_bonds(resvecs):
        LJ.ignore(a1,a2)
    
    for a1,a2,a3,angle in Resvec.all_angles(resvecs):
        #~ LJ.ignore(a1,a2)
        #~ LJ.ignore(a2,a3)
        LJ.ignore(a1,a3)
    
    return LJ

def make_charges(resvecs, screen, ph=7.4, k=1):
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
    atoms = [atom for r in resvecs for atom in r
                        if (r.resname, atom.name) in chargedict]
    
    residues = set(a.group for a in atoms)
    import collections
    d = collections.defaultdict(int)
    for r in residues: d[r.resname] += 1
    print(len(atoms), "charged atoms")
    for r, n in sorted(d.items()):
        print(n, r)
    
    for atom in atoms:
        #~ print('ref for', atom)
        aref = fullgroup.get_id(atom)
        charge = chargedict[(atom.group.resname, atom.name)]
        charges.add(aref, charge)
    
    for a1,a2,l in Resvec.all_bonds(resvecs):
        key1, key2 = (a1.group.resname, a1.name), (a2.group.resname, a2.name)
        if key1 not in chargedict or key2 not in chargedict:
            continue
        #~ print('ignoring', a1.group.resname, a1.name, a2.name)
        charges.ignore(a1,a2)
    
    for a1,a2,a3,angle in Resvec.all_angles(resvecs):
        key1, key2 = (a1.group.resname, a1.name), (a3.group.resname, a3.name)
        if key1 not in chargedict:
            continue
        #~ print("looking at", a1.group.resname, a1.name, a3.name)
        if key2 not in chargedict:
            #~ print("not ignoring", a1.group.resname, a1.name, a3.name)
            continue
        #~ print('ignoring', a1.group.resname, a1.name, a3.name)
        charges.ignore(a1,a3)
    
    return charges
