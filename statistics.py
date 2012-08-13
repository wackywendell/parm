import gzip
import simw as sim

class Stat(object):
    AllStats = {}
    def __init__(self, calcfunc, headersfunc=None, name=None):
        self.name = calcfunc.func_name if name is None else name
        self.calculate = calcfunc
        if headersfunc is not None: self.headers = headersfunc
        if name in self.AllStats and self.AllStats[self.name] != self:
            raise ValueError("Stat cannot register " + repr(calcfunc) + 
                    ", key " + repr(name) + " already taken by " +
                    repr(self.AllStats[name]))
        self.AllStats[self.name] = self

    def calculate(self, collec, reslists):
        raise NotImplementedError
    
    def headers(self, reslists):
        return [self.name]
    
    def headerfunc(self, newfunc):
        self.headers = newfunc

def RStat(func):
    name = func.func_name
    def headers(reslists):
        return ([name] if len(reslists) <= 1 else
            ['%s_%d' % (name, n) for n in range(len(reslists))])
    
    s = Stat(func, headers)
    #@s.headerfunc
    #s.headers = headers
    return s

@Stat
def T(collec, reslists):
    return [collec.temp()]

@Stat
def E(collec, reslists):
    return [collec.energy()]

@Stat
def L(collec, reslists):
    return [collec.angmomentum().mag()]

@Stat
def Kinetic(collec, reslists):
    return [collec.kinetic()]

@Stat
def COMV(collec, reslists):
    return [collec.comv().mag()]

@RStat
def Rg(collec, reslists):
    return [sim.calc_Rg(reslist) for reslist in reslists]

def Rijs(ijs):
    def Rijs(collec, reslists):
        return [(reslist[i]['CA'].x - reslist[j]['CA'].x).mag()
                    for reslist in reslists for i,j in ijs]
    def headers(reslists):
        ijstrs = ['%d-%d' % ij for ij in ijs]
        if len(reslists) == 1: return ijstrs
        else: return [ijstr + '_' + str(n) 
            for n,reslist in enumerate(reslists) for ijstr in ijstrs]
    return Stat(Rijs, headers, 'Rijs')

#Rijs([(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
#            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)])


class TableWriter(object):
    def __init__(self, fname, headers=None, gzip=True, sep='\t'):
        self.fname = fname
        self.gzip = gzip
        self.sep = sep
        self.headers = headers
        self.headers_written = False
        assert all([isinstance(s, str) for s in headers]), headers
    
    def write_row(self, cols):
        txt = '\n' + '\t'.join(cols)
        mode = 'a'
        if not self.headers_written:
            txt = '\t'.join(self.headers) + txt
            self.headers_written = True
            mode = 'w'
        
        with gzip.open(self.fname, mode) as f:
            f.write(txt)

class StatWriter(TableWriter):
    def __init__(self, fname, collec, reslists, stats=None, gzip=True,
                    contin=False):
        if stats == None:
            names, stats = zip(*sorted(Stat.AllStats.items()))
        headers = ['time'] + [h for s in stats for h in s.headers(reslists)]
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
        self.stats = stats
        if contin: self.headers_written=True
    
    def update(self, t):
        def tostr(x):
            if isinstance(x, float): return '%.6g' % x
            return str(x)
        statstrs = ['%d' % t] + [tostr(x) 
                        for s in self.stats
                        for x in s.calculate(self.collec, self.reslists)]
        assert len(statstrs) == len(self.headers)
        self.write_row(statstrs)

def loadTable(fname):
    import numpy
    arr = numpy.genfromtxt(fname, names=True)
    names = arr.dtype.names
    return [(nm, arr[nm]) for nm in names]
