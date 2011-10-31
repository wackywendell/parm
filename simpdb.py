from __future__ import print_function
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
        else:
            self.atoms = residue.child_list
        self.masses = [self._get_mass(residue, a, H)*amu for a in self.atoms]
        
        atomvec.__init__(self, self.masses)
        self.set_locs(residue)
        self.loadfile = loadfile
        self.resname = residue.resname
    
    # dictionaries with means, to be loaded as necessary
    backbonds = None
    resbonds = None
    angles = None
    
    def set_locs(self, residue):
        atoms = self.atoms
        if(len(atoms) != len(self)):
            raise TypeError("Resvec.setLocs: lengths do not match")
        locs = [a.get_coord().tolist() for a in atoms]
        for loc, atom in zip(locs, self):
            atom.x = Vec(*loc)
    
    def _get_mass(self, res, atom, H = True):
        if not H:
            return atom.mass
        name = atom.name[:]
        names = (name, name + '1', name + '2', name + '3', name + '4')
        Hs = [a.mass for a in res.child_list if a.name in names]
        return sum(Hs) + a.mass
    
    def load_data(self, f=None):
        if None not in (self.resbonds, self.backbonds, self.angles):
            return
        if f is None:
            f = self.loadfile
            if self.loadfile is None:
                raise KeyError("No file or filename provided")
        
        if not (hasattr(f, 'read') and hasattr(f, 'readline')):
            f = open(f, 'rb')
        
        import cPickle as pickle
        self.resbonds, self.backbonds, self.angles = pickle.load(f)
    
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
        a.group = self
        return a
    
    #~ def __iter__(self):    
        #~ for i in range(self.N()):
            #~ yield self[i]
    
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

def make_structure(residue_list, bond_k, loadfile, amu=1, angstrom=1):
    """Returns a tuple (atomvecs, bondgroup)
    
    atomvecs is a list of atomvecs corresponding to residues.
    bondgroup is a bondgrouping interaction corresponding to all the bonds,
    backbone included."""
    #so that it works on a chain object as well
    if hasattr(residue_list, 'child_list'):
        residue_list = residue_list.child_list
    lastres = None
    
    avecs = [Resvec(res, amu=amu, loadfile=loadfile) for res in residue_list]
    
    bond_pairs = bondpairs()
    for a1,a2,l in Resvec.all_bonds(avecs):
        length = l * angstrom
        bond_pairs.add(bond_k, length, a1, a2)
    
    return avecs, bond_pairs
