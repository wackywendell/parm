from simw import *
from fpath import File, Path, Dir
import xyzfile
import yaml
import numpy
import sys

def getDIHdists(fname, bins, cutoff, template=None, outf=None, force=False, printevery=None):
    fname = File(fname)
    if not fname.exists():
        raise NotImplementedError
        
    if outf is None:
        extlen = len(fname.extension) + 1
        mn, mx = ('%.2f' % min(bins)), ('%.2f' % max(bins))
        #print('min-max:', mn,mx)
        extra = '-' + 'dihedrals' + ('-%s-%s_%d_%d' % (mn,mx, len(bins), cutoff))
        #print(extra)
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
    
    dists = []
    
    atomquads = [[r['CA'] for r in rs] 
                for rs in zip(rvecs, rvecs[1:], rvecs[2:], rvecs[3:])]
    
    with open(str(fname), 'r') as f:
        xyz = xyzfile.XYZreader(f)
        for n, f in enumerate(xyz):
            if n < cutoff: continue
            if printevery and (n % printevery) == 0 and n > 0:
                print(n, end = ' ')
                sys.stdout.flush()
            f.into(atoms)
            dists.extend([dihedral.getang(
                                a2.x - a1.x, a3.x - a2.x, a4.x - a3.x)
                            for a1, a2, a3, a4 in atomquads])
    if printevery:
        print('Done.')
        sys.stdout.flush()
    
    hst, _ = numpy.histogram(dists, bins)
    
    numpy.save(str(outf), hst)
    return hst

if __name__ == "__main__":
    pass
    #~ xs = numpy.linspace(0,30,101)
    #~ xsfine = numpy.linspace(0, 30, 201)
    #~ cutoff=10000
    #~ 
    #~ import sys
    #~ if 2 <= len(sys.argv) <= 3:
        #~ fname = sys.argv[1]
        #~ template = sys.argv[2] if len(sys.argv) == 3 else None
        #~ 
        #~ getCOMdists(fname, xs, cutoff=cutoff, template=template)
        #~ getCOMdists(fname, xsfine, cutoff=cutoff, template=template)
    #~ else:
        #~ print("Wrong number of arguments (%d)" % len(sys.argv))
