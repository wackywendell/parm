from __future__ import print_function
from fpath import *
from Bio.PDB import PDBList
from urllib2 import URLError

mydir = Dir('~/idp/').norm()

pdbcsvfile = File(mydir + 'blengths/Dunbrack_bbdep02.May.cmpdlist').norm()

with pdbcsvfile.open('r') as f:
    names = [line[:4] for line in f if line[:4].strip()]

print('Found', len(names), 'Names')

pdbl = PDBList()

broken = []

count=0
for i, n in enumerate(names):
    print(i, end=': ')
    try:
        pdbl.retrieve_pdb_file(n, pdir = str(mydir + 'allpdbs'))
        count += 1
    except URLError:
        broken.append(n)
        print(n, 'BROKEN')

print(count, "Done!", len(broken), 'Broken:', *broken)
