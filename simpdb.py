from __future__ import print_function
from Bio.PDB import PDBParser
from Bio.PDB.Residue import Residue as _Residue
from simw import *
import itertools


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

# from richards
LJdict = {
    'CH3' : 2.0, 'CH2' : 2.0, 'CH' : 2.0, 'C' : 1.7, 'O' : 1.4,
    'OH' : 1.6, 'NH3' : 2.0, 'NH2' : 1.75, 'NH' : 1.7, 'S' : 1.8
    }

# following RasMol
bonddists = {
    'CC':2,'CN':1.96,'CO':1.96,'PP':2.632,'OP':2.276,'SS':2.6,'OS':2.26,
    'CaO':2.232,'SZn':3.028
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
    
    def load_data(self, f=None):
        if None not in (self.resbonds, self.backbonds, self.resangles, self.backangles):
            return
        if f is None:
            f = self.loadfile
            if self.loadfile is None:
                raise KeyError("No file or filename provided")
        
        if not (hasattr(f, 'read') and hasattr(f, 'readline')):
            f = open(f, 'rb')
        
        import cPickle as pickle
        self.resbonds, self.backbonds, self.resangles, self.backangles = pickle.load(f)
    
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
        reslist = list(reslist)
        lastvecs = [None] + reslist
        lists_of_bonds = (res.get_bonds(lvec) for res, lvec
                                        in zip(reslist, lastvecs))
        
        return itertools.chain.from_iterable(lists_of_bonds)
    
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
    def from_pdb(cls, strucname, pdbfilename, loadfile, amu=1):
        pdbp = PDBParser()
        aS = pdbp.get_structure(strucname, pdbfilename)
        chains = list(aS.get_chains())
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

def make_structure(residue_list, loadfile, bond_k, angle_k, epsilon=1, amu=1, angstrom=1):
    """Returns a tuple (Resvecs, list of interactions, list of trackers)
    
    Resvecs is a list of Resvecs (atomvecs) corresponding to residues.
    bondgroup is a bondgrouping interaction corresponding to all the bonds,
    backbone included."""
    #so that it works on a chain object as well
    if hasattr(residue_list, 'child_list'):
        residue_list = residue_list.child_list
    lastres = None
    
    avecs = [Resvec(res, amu=amu, loadfile=loadfile) for res in residue_list]
    fullgroup = metagroup(avecs)
    
    bond_pairs = bondpairs()
    bond_angles = angletriples()
    maxsigma = max(LJdict.values())
    neighbors = neighborlist(fullgroup, 2.5*maxsigma, 4*maxsigma)
    LJ = LJgroup(neighbors, 2.5)
    for i, atom in zip(range(fullgroup.size()), (atom for r in avecs for atom in r)):
        assert fullgroup.get_id(i) == atom
        sigma = LJdict[atom.formula]
        LJ.add(fullgroup.get_id(i), sigma, epsilon)
    
    for a1,a2,l in Resvec.all_bonds(avecs):
        length = l * angstrom
        bond_pairs.add(bond_k, length, a1, a2)
        neighbors.ignore(a1,a2)
    
    for group in Resvec.all_angles(avecs):
        a1,a2,a3,angle = group
        bond_angles.add(angle_k, angle, a1, a2, a3)
        neighbors.ignore(a1,a3)
    
    neighbors.update_list(True)
    
    #~ return avecs, [bond_pairs, bond_angles], [neighbors]
    return avecs, [bond_pairs, bond_angles, LJ], [neighbors]
