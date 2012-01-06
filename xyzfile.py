from __future__ import print_function
import re, math, numpy as np
from itertools import izip
import hashlib, logging
#~ import pyparsing as parse

class XYZwriter:
    def __init__(self, f):
        self.file = f
    
    def writeframe(self, reslist, comment='', com=None):
        Natoms = sum(len(r) for r in reslist)
        print(Natoms, file=self.file)
        print(comment, file=self.file)
        for r in reslist:
            for atom in r:
                elem = atom.element
                if com is None:
                    x,y,z = tuple(atom.x)
                else:
                    x,y,z = tuple(atom.x - com)
                print("%s %.3f %.3f %.3f" % (elem,x,y,z), file=self.file)
        self.file.flush()
    
    def size(self):
        location = self.file.tell()
        self.file.seek(0,2)
        size = self.file.tell()
        self.file.seek(location)
        return size
    
    def close(self):
        self.file.close()

class XYZError(Exception):
    pass

def md5_file(f, block_size = None):
    loc = f.tell()
    f.seek(0)
    md5 = hashlib.md5()
    if block_size is None:
        block_size = md5.block_size
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    f.seek(loc)
    return md5.hexdigest()

#~ parse.ParserElement.setDefaultWhitespaceChars(" \t")
#~ real = parse.Word(parse.nums) + parse.Optional("." + parse.Word(parse.nums))
#~ real = real.setParseAction(lambda toks: float("".join(toks)))
#~ intparse = parse.Word(parse.nums).setName('int').setParseAction(lambda toks: int(toks[0]))
#~ space = parse.Suppress(' ')
#~ endln = parse.Suppress(parse.LineEnd())
#~ line = parse.Group((parse.Word(parse.alphas) + real + real + real + endln))
#~ manylines = parse.Group(parse.OneOrMore(line))
#~ 
#~ firstline = intparse + endln
#~ secondline = parse.Suppress('time ') + parse.restOfLine + endln
#~ fullframe = firstline + secondline + manylines
#~ fileparse = parse.OneOrMore(fullframe)

class XYZreader:
    firstline = re.compile("([0-9]*)")
    commentline = re.compile("time ([0-9.]*)")
    regline = re.compile("(\w{1,2})\s+([\-0-9.]+)\s+([\-0-9.]+)\s+([\-0-9.]+)")
    
    def __init__(self, f):
        self.file = f
        self._size = None
        self._md5 = None
    
    def __match__(self, regex, *types):
        line = self.file.readline()
        if line == '':
            raise StopIteration, "Reached end of file."
        line = line.strip()
        match = regex.match(line)
        if match is None and types:
            raise XYZError("Failed to parse at line " + repr(line))
        groups = match.groups()
        if len(groups) != len(types):
            raise XYZError("Wrong number of types given.")
        return tuple(t(val) for t,val in zip(types,groups))
    
    def readframe(self):
        num, = self.__match__(self.firstline, int)
        t, = self.__match__(self.commentline, float)
        #~ try:
            #~ lines = [self.regline.match(self.file.readline()).groups()
                #~ for i in range(num)
            #~ ]
        #~ except AttributeError:
            #~ raise StopIteration("Frame broken")

        lines = [self.file.readline().strip().split(' ') for i in range(num)]
        lines = [(e,(float(x),float(y),float(z))) for e,x,y,z in lines]
        return Frame(lines, t)
        
        lines = "".join([self.file.readline() for i in range(num)])
        return Frame()
        
    def __iter__(self):
        location = self.file.tell()
        while True:
            yield self.readframe()
        self.file.seek(location)
    
    def close(self):
        self.file.close()
    
    def all(self):
        frames = Frames(iter(self))
        for n,f in enumerate(frames):
            f.indx = n
        return frames
    
    def md5(self):
        if self._md5 is None:
            logging.debug('calculating md5...')
            self._md5 = md5_file(self.file, 2**16)
            logging.debug('md5: %s' % self._md5)
        return self._md5
    
    def size(self):
        if self._size is None:
            location = self.file.tell()
            self.file.seek(0,2)
            self._size = self.file.tell()
            self.file.seek(location)
        return self._size

class Velreader(xyzreader):
    def readframe(self):
        f = xyzreader.readframe(self)
        return VelFrame(f)
        
class Frame:
    def __init__(self, locs, t=None, indx=None):
        self.t = t
        self.indx = None
        self.elems, self.locs = zip(*locs)
        self.__locarray = None
    
    def __iter__(self):
        return iter(self.locs)
    
    def into(self, atoms, check=True):
        if not check:
            for a, loc in zip(atoms, self.locs):
                #~ print(type(a),a)
                a.x.set(*loc)
            return
        
        for a, elem, loc in izip(atoms, self.elems, self.locs):
            if a.element != elem:
                raise TypeError, ("Element mismatch for atom %d (%s) and line %s" 
                                    % (num, atom.element, repr(line)))
            a.x.set(*loc)
    
    @property
    def time(self):
        return self.t
    
    @property
    def locarray(self):
        if self.__locarray is None:
            #~ self.__constructed += 1
            #~ print("constructing", self.time, self.__constructed)
            self.__locarray = np.array(list(self.locs))
        #~ else:
            #~ print("FOUND")
        return self.__locarray

def VFrame(Frame):
    def into(self, atoms, check=True):
        if not check:
            for a, loc in zip(atoms, self.locs):
                a.v.set(*loc)
            return
        
        for a, elem, loc in izip(atoms, self.elems, self.locs):
            if a.element != elem:
                raise TypeError, ("Element mismatch for atom %d (%s) and line %s" 
                                    % (num, atom.element, repr(line)))
            a.v.set(*loc)

class Frames(list):
    def into(self, atoms, check=True):
        for f in self:
            f.into(atoms, check)
            yield f.time

class FramePairs(tuple):
    """Should be Framepairs(FrameList, VFrameList)"""
    def into(self, atoms, check=True):
        for f, vf in zip(*self):
            assert vf.time = f.time
            f.into(atoms, check)
            vf.into(atoms, check)
            yield f.time

import re
refirst = '(\d*)\n'
resecond = 'time (\d*\.?\d+)(.*)\n'
reline = '(\w+) (\d*\.?\d+) (\d*\.?\d+) (\d*\.?\d+)\n'
remany = '(?:(\w+) (\d*\.?\d+) (\d*\.?\d+) (\d*\.?\d+)\n)+'


if __name__ == '__main__':
    s = open('test/T1-1K.xyz','r').read()
    print('got s')
    #~ v1 = re.match(refirst, s)
    #~ v12 = re.match(refirst + resecond, s)
    #~ v13 = re.match(refirst + resecond + reline + reline, s, re.MULTILINE)
    #~ 
    #~ print('found', len(v13.groups()))
    
    firstv = firstline.parseString(s)
    firstsecv = (firstline + secondline).parseString(s)
    v13 = (firstline + secondline + line).parseString(s)
    vs = fullframe.parseString(s)
    
    if False:
        import cProfile as profile
        from simw import atomvec
        f=open('/home/wendell/idp/data/T3-200K.xyz','r')
        x=XYZreader(f)
        #~ print(x.readframe())
        #~ profile.runctx('for l in x: pass', globals(), locals())
        atoms = atomvec([1]*1013)
        #~ profile.runctx('for t,l in x: Frame(l,t).into(atoms,False)', globals(), locals())
        frames = [Frame(l,t) for t,l in x]
        profile.runctx('for f in frames: f.into(atoms,False)', globals(), locals())
