# encoding: UTF8

from Bio.PDB import PDBParser, Vector
from Bio.PDB.Residue import Residue as _Residue
from Bio.PDB.Atom import Atom as _Atom
from simw import *
from proteindata import *
import itertools
from collections import defaultdict
import math

# adapted from IUPAC standards, available at
# http://www.chem.qmul.ac.uk/iupac/misc/ppep4.html#400

from os.path import expanduser
mydir = expanduser('~/idp/')
defaultpdbfile  = mydir + 'pdb/aS.pdb'
defaultloadfile = mydir + 'blengths/stats-H.pkl'


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

class Resvec(atomvec):
    def __init__(self, residue, H=False, amu=1, loadfile=None):
        #residue.sort()
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
    
    def get_n_bonds(self, atom):
        self.load_data()
        rbonds = self.resbonds[self.resname]
        mybonds = ([b for b in rbonds if atom.name in b] +
                    [b for b in self.backbonds if atom.name in b])
        return len(mybonds)
        
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
        
        #coordoffset =  avg(dvecs) # can't divide vector by float? what?
        dvecsum = dvecs[0]
        for v in dvecs[1:]:
            dvecsum += v
        coordoffset = dvecsum.__div__(len(dvecs)) # __div__ doesn't work for '/' in py3?
        
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

def make_LJ(box, resvecs, epsilon=1, sizefactor=1, neighborcutoff=1.4, LJdict=AliceSizes):
    """Neighborcutoff relative to sigma.
    
    returns (LJ, neighbors)"""
    factor = 2 * sizefactor * (2**(1.0/6.0))
    # 2 because LJ-sigma is the sum of two radii, (2**(1.0/6.0)) for the shape
    
    maxsigma = max(LJdict.values()) * factor
    innerradius = maxsigma
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(box, innerradius, outerradius)
    LJ = LJgroup(neighbors)
    
    keygens = (lambda a,r: a.formula, lambda a,r: atom.element, 
                lambda a,r: atom.element + str(r.get_n_bonds(atom)),
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
        core = getsize(atom, res)
        sigma = core * factor
        LJ.add(LJatom(epsilon, sigma, atom))
    
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
        #~ charge = chargedict[(atom.group.resname, atom.name)]
        charges.add(atom, charge)
    
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

def make_hydrophobicity(box, resvecs, epsilon, LJcutoff=2.5, neighborcutoff=1.4,
                    sigma=5.38, hydroindex=hydroindex7, method='geometric', scale='max'):
    """LJcutoff in sigma units, neighborcutoff relative to LJcutoff.
    Makes a CA hydrophobicity attractive force.
    
    returns (interaction, neighbors)"""
    
    # 2 because LJ-sigma is the sum of two radii, (2**(1.0/6.0)) for the shape
    
    innerradius = sigma*LJcutoff
    outerradius = innerradius * neighborcutoff
    neighbors = neighborlist(box, innerradius, outerradius)
    Hphob = Hydrophobicity(neighbors)
    
    #~ print("Hlist:", *Hlist)
    Rnames, Hlist = list(zip(*sorted(hydroindex.items())))
    
    if scale.lower() == 'none':
        pass
    elif scale.lower() == 'max':
        mx = max([abs(h) for h in Hlist])
        Hlist = [h / mx for h in Hlist]
    elif scale.lower() == 'minmax':
        mn, mx = min(Hlist), max(Hlist)
        Hlist = [(h - mn) / (mx - mn) for h in Hlist]
    elif scale.lower() == 'zeroone':
        mx = max([abs(h) for h in Hlist])
        Hlist = [((h/mx) + 1) / 2 for h in Hlist]
        assert min(Hlist) >= 0
    
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
    
    allHs = [combine(H1, H2) for H1 in Hlist for H2 in Hlist]
    print("Hydrophobicity made with %s-%s: %.2f - %.2f (x %.2f)" % 
            (scale, method, min(allHs), max(allHs), epsilon))
    
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
        return HydroAtom(epsvec, idx, sigma, atom, LJcutoff)
    
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
