import gzip

class Stat(Object):
    AllStats = {}
    def __init__(self, calcfunc, headersfunc=None, name=None):
        if name is None: name = calcfunc.func_name
        self.calculate = calcfunc
        if headersfunc is not None: self.headers = headersfunc
        if name in self.AllStats and self.AllStats[name] != self:
            raise ValueError("Stat cannot register " + repr(stat) + 
                    ", key " + repr(name) + " already taken by " +
                    repr(cls.AllStats[name]))
        self.AllStats[name] = self

    def calculate(self, collec, reslists):
        raise NotImplementedError
    
    def headers(self, reslists):
        return [self.name]
    
    def headerfunc(self, newfunc):
        self.headers = newfunc

def RStat(func):
    s = Stat(func)
    @s.headerfunc
    def headers(self, reslists):
        return ([self.name] if len(reslists) <= 1 else
            ['%s_%d' % (self.name, n) for n in range(reslists)])
    return s

@Stat
def T(self, collec, reslists):
    return [collec.temp()]

@Stat
def E(self, collec, reslists):
    return [collec.energy()]

@Stat
def L(self, collec, reslists):
    return [collec.angmomentum().mag()]

@Stat
def Kinetic(self, collec, reslists):
    return [collec.kinetic()]

@Stat
def COMV(self, collec, reslists):
    return [collec.comv().mag()]

@RStat
def Rg(self, collec, reslists):
    return [sim.calc_Rg(reslist) for reslist in reslists]

def Rijs(ijs):
    def Rijs(self, collec, reslists):
        return [(reslist[i]['CA'] - reslist[j]['CA']).mag()
                    for reslist in reslists for i,j in ijs]
    def headers(self, reslists):
        ijstrs = ['%d-%d' % ij for ij in ijs]
        if len(reslists) == 1: return ijs
        else: return [ijstr + '_' + str(n) 
            for n,reslist in enumerate(reslists) for ijstr in ijstrs]
    return Stat(Rijs, headers, 'Rijs')

class TableWriter(object):
    def __init__(self, fname, headers=None, gzip=True, sep='\t'):
        self.fname = fname
        self.gzip = gzip
        self.sep = sep
        self.headers = headers
        self._headers_written = False
    
    def write_row(self, cols):
        txt = '\n' + '\t'.join(cols)
        mode = 'a'
        if not self._headers_written:
            txt = '\t'.join(self.headers) + txt
            self._headers_written = True
            mode = 'w'
        
        with gzip.open(self.fname, mode) as f:
            f.write(txt)

class StatWriter(TableWriter):
    def __init__(self, fname, collec, reslists, stats=None, gzip=True):
        if stats == None:
            names, stats = zip(*sorted(Stat.AllStats))
        headers = [h for h in s.headers(reslists) for s in stats]
        TableWriter.__init__(self, fname, headers, gzip)
        self.collec = collec
        self.reslists = reslists
    
    def update(self):
        def tostr(x):
            if isinstance(x, float): return '%.6g' % x
            return str(x)
        statstrs = [tostr(x) 
                        for x in s.calculate(self.collec, self.reslists)
                        for s in stats]
        assert len(statstrs) == len(self.headers)
        self.write_row(statstrs)

def loadTable(fname):
    import numpy
    arr = numpy.genfromtxt(fname, names=True)
    names = arr.dtype.names
    return [(nm, arr[nm]) for nm in names]
