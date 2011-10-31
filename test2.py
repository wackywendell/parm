from __future__ import print_function
from Bio.PDB import PDBParser
from simw import *
from simpdb import Resvec, Residue, make_structure





pdbp = PDBParser()
aS = pdbp.get_structure('AS','/home/wendell/idp/pdb/aS.pdb')
loadfile='../blengths/stats.pkl'
chain, = aS.get_chains()
#~ 
#~ r = Resvec(chain.child_list[-1], loadfile=loadfile)
#~ print(r[0].x)
#~ print('C', r['C'].x)
#~ print(r.resname, [a.name for a in r.atoms])
#~ print(len(list(r.get_bonds())), len(r.resbonds[r.resname]))
#~ 
#~ newbs = [(r.get_name(a1), r.get_name(a2), l) for a1, a2, l in r.get_bonds()]
#~ print(newbs)
#~ 
#~ r2 = Resvec(chain.child_list[-2], loadfile=loadfile)

print("Making structure...")
avecs, bondpair = make_structure(chain, 100, loadfile)
print("Sizes:",len(avecs), bondpair.size())

a = avecs[0][0]
import pickle
s = pickle.dumps(a)
b = pickle.loads(s)

if False:
    r0 = Residue(chain.child_list[-1])

    avec = r0.make_atomvec()

    print(sum(r0.masses(False)), sum(r0.masses(True)))

    print([a.name for a in r0.atoms(True)])

    for a, atom, mass in zip(avec, r0.atoms(False), r0.masses(False)):
        print(atom.name, a.x, mass, [h.name for h in r0.get_Hs(atom)])
