import os, os.path
from fpath import File, Path, Dir
import numpy

def loadTSV(fname):
    fname = File(fname)
    with fname.open('r') as f:
        l = f.readline()
        while l[0] == '#':
            l = f.readline()
        colnames = l.strip('\n').split('\t')
        matr = numpy.loadtxt(f)
    assert len(matr.T) == len(colnames)
    return dict(zip(colnames, matr.T))
    
def makeHist(fname, colname, bins, cutoff=0, outf=None, force=False):
    fname = File(fname)
    if outf is None:
        extlen = len(fname.extension) + 1
        extra = '-' + colname + 'Hist' + ('-%d-%d_%d_%d' % (min(bins), max(bins), len(bins), cutoff))
        outf = (fname[:-1]) + (fname[-1][:-extlen] + extra + '.npy')
    outf = File(outf)
    
    if not force and outf.exists() and outf.stat().mtime > fname.stat().mtime:
        # don't need to remake it
        return numpy.load(str(outf))
    
    # need to make histogram, so read in the big file
    fullf = loadTSV(fname)
    col = fullf[colname]
    hst, _ = numpy.histogram(col[cutoff:], bins)
    
    numpy.save(str(outf), hst)
    return hst
        
    
    
    
