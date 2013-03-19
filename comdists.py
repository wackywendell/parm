from simw import *
from fpath import File, Path, Dir
import xyzfile
import yaml
import numpy

def getCOMdists(fname, bins, cutoff, template=None, outf=None, force=False):
    fname = File(fname)
        
    if outf is None:
        extlen = len(fname.extension) + 1
        mn, mx = int(min(bins)+0.5), int(max(bins)+0.5)
        #print('min-max:', mn,mx)
        extra = '-' + 'COMHist' + ('-%d-%d_%d_%d' % (mn,mx, len(bins), cutoff))
        #print(extra)
        outf = (fname[:-1]) + (fname[-1][:-extlen] + extra + '.npy')
    outf = File(outf)
    
    if not force and outf.exists() and outf.stat().mtime >= fname.stat().mtime:
        # don't need to remake it
        return numpy.load(str(outf))
    else:
        print("Failed:", force, outf.exists,
                outf.stat().mtime >= fname.stat().mtime,
                outf.stat().mtime, fname.stat().mtime)
        
    
    # load file
    if template == None:
        template = 'startlocs/%s-eq.yaml' % (fname[-1][:3])
    
    with open(str(template),'r') as f:
        resparams = yaml.load(f)
    resdict, = resparams['atoms']
    anames, masses_indices = zip(*sorted(resdict['sizes'].items()))
    masses, radii, indices = zip(*masses_indices)
    res1 = Res(anames, masses)
    res2 = Res(anames, masses)
    
    atoms = list(res1) + list(res2)
    
    dists = []
    with open(str(fname), 'r') as f:
        xyz = xyzfile.XYZreader(f)
        for n, f in enumerate(xyz):
            if n < cutoff: continue
            f.into(atoms)
            dists.append((res1.atomvec.com() - res2.atomvec.com()).mag())
    
    hst, _ = numpy.histogram(dists, bins)
    
    numpy.save(str(outf), hst)
    return hst

if __name__ == "__main__":
    xs = numpy.linspace(0,30,101)
    xsfine = numpy.linspace(0, 30, 201)
    cutoff=10000
    
    import sys
    if 2 <= len(sys.argv) <= 3:
        fname = sys.argv[1]
        template = sys.argv[2] if len(sys.argv) == 3 else None
        
        getCOMdists(fname, xs, cutoff=cutoff, template=template)
        getCOMdists(fname, xsfine, cutoff=cutoff, template=template)
    else:
        print("Wrong number of arguments (%d)" % len(sys.argv))
