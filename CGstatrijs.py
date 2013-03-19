from simw import *
from fpath import File, Path, Dir
import xyzfile
import yaml
import numpy
import sys
from hashlib import sha256

def getRijs(fname, cutoff, ijs, template=None, outf=None, force=False, printevery=None):
    fname = File(fname)
    if not fname.exists():
        raise NotImplementedError
        
    if outf is None:
        extlen = len(fname.extension) + 1
        try:
            hash = sha256(str(list(map(tuple, ijs))).encode('utf8')).hexdigest()
        except TypeError:
            print('ijs:', type(ijs), ijs)
            raise
        extra = '-' + 'Rij' + ('-%s-%d' % (hash[-8:], cutoff))
        outf = (fname[:-1]) + (fname[-1][:-extlen] + extra + '.npy')
    outf = File(outf)
    
    if not force and outf.exists() and outf.stat().mtime >= fname.stat().mtime:
        # don't need to remake it
        return numpy.load(str(outf))
    #~ else:
        #~ print("Failed:", force, outf.exists(),
                #~ outf.stat().mtime >= fname.stat().mtime if outf.exists() else '',
                #~ outf.stat().mtime if outf.exists() else '',
                #~ fname.stat().mtime if outf.exists() else '')
        
    
    # load file
    if template == None:
        template = 'startlocs/cg/RW.yml'
    
    with open(str(template),'r') as f:
        resparams = yaml.load(f)
    rvecs = []
    for rdict in resparams['structure']:
        name = rdict['type']
        atoms = rdict['atoms']
        anames, masses = zip(*sorted([(k, v['mass']) for k,v in atoms.items()]))
        
        res = Res(anames, masses)
        for atom in res:
            adict = atoms[atom.name]
            x,y,z = adict['xyz']
            atom.x = Vec(x,y,z)
        rvecs.append(res)
    
    atoms = [a for r in rvecs for a in r]
    atompairs = [(rvecs[i]['CA'], rvecs[j]['CA']) for i,j in ijs]
    
    dists = []
    with open(str(fname), 'r') as f:
        xyz = xyzfile.XYZreader(f)
        for n, f in enumerate(xyz):
            if n < cutoff: continue
            if printevery and (n % printevery) == 0 and n > 0:
                print(n, end = ' ')
                sys.stdout.flush()
            f.into(atoms)
            dists.append([(a1.x - a2.x).mag() for a1,a2 in atompairs])
    
    if printevery:
        print('Done.')
        sys.stdout.flush()
    
    dists = numpy.array(dists, dtype=float)
    
    numpy.save(str(outf), dists)
    return dists

if __name__ == "__main__":
    cutoff=100
     
    import sys, FRETvals
    if 2 <= len(sys.argv) <= 3:
        fname = sys.argv[1]
        template = sys.argv[2] if len(sys.argv) == 3 else None
        
        getRijs(fname, cutoff=cutoff, ijs=FRETvals.ijs, template=template, printevery=100)
    else:
        print("Wrong number of arguments (%d)" % len(sys.argv))
