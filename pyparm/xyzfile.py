import re, math, numpy as np
import hashlib, logging
import pyparm.d3 as sim

class XYZwriter:
    def __init__(self, f, usevels = True):
        self.file = f
        self.usevels = usevels
    
    def writeframe(self, atoms, com=None, box=None, **kwargs):
        if 'time' not in kwargs:
            raise ValueError("Time must be in keyword arguments, or .xyz will not be readable")
        Natoms = len(atoms)
        print(Natoms, file=self.file)
        comment = " ".join(["=".join((k,str(v))) for k,v in list(kwargs.items())])
        print(comment, file=self.file)
        
        lines = []
        for a in atoms:
            elem = a.element
            if com is None:
                x,y,z = tuple(a.x) if box is None else tuple(box.diff(a.x, sim.Vec.Zero()))
            else:
                x,y,z = tuple(a.x - com)  if box is None else tuple(box.diff(a.x, com))
            
            line = [elem] + ['%.4f' % coord for coord in (x,y,z)]
            if self.usevels:
                vx,vy,vz = tuple(a.v)
                line.extend(['%.4f' % coord for coord in (vx,vy,vz)])
            if hasattr(a, 'sigma'):
                line.append('%.6f' % a.sigma)
            
            print(' '.join(line), file=self.file)
                
        
        # do all the writing at once
        #print('\n'.join(lines), file=self.file)
        self.file.flush()
    
    def writefull(self, t, atoms, collec, com=None):
        cdict = {
            'time':t,
            'E':collec.energy(),
            'T':collec.temp(),
            'K':collec.kinetic(),
            'L':collec.angmomentum().mag(),
            'v':collec.comv().mag(),
            }
        
        if hasattr(collec, 'interactions'):
            for name,interac in list(collec.interactions.items()):
                cdict[name + 'E'] = interac.energy(collec.getbox())
            for i in list(collec.interactions.values()):
                if isinstance(i, sim.BondPairs):
                    cdict[name + 'mean'] = i.mean_dists()
                    cdict[name + 'std'] = i.std_dists()
                if isinstance(i, sim.AngleTriples):
                    cdict[name + 'mean'] = i.mean_dists()
                    cdict[name + 'std'] = i.std_dists()
                
        if com is None: com = collec.com()
        self.writeframe(atoms, com, **cdict)
        return cdict
    
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
    firstline = re.compile("^([0-9]*)$")
    commentline = re.compile("time[= ][0-9.]*|[\w-]+=.?( [\w-]+=.?)*|[0-9.]*")
    regline = re.compile("(\w{1,2})\s+([\-0-9.]+)\s+([\-0-9.]+)\s+([\-0-9.]+)")
    
    def __init__(self, f):
        self.file = f
        self._size = None
        self._md5 = None
        self.lineno = 1
    
    def __match__(self, regex, *types):
        line = self.file.readline()
        if line == '':
            raise StopIteration("Reached end of file.")
        line = line.strip()
        match = regex.match(line)
        if match is None:
            rname = "regex '{}'".format(repr(regex.pattern))
            if regex == self.firstline:
                rname = "first frame line"
            elif regex == self.commentline:
                rname = "comment line" 
            raise XYZError("Failed to parse {} at line {}: {}".format(rname, self.lineno, repr(line)))
        self.lineno += 1
        groups = match.groups()
        if len(types) != 0 and len(groups) != len(types):
            raise XYZError("Wrong number of types given.")
        if len(types) == 0:
            return line
        return tuple(t(val) for t,val in zip(types,groups))
    
    def readframe(self):
        num, = self.__match__(self.firstline, int)
        comment = self.__match__(self.commentline)
        if '=' in comment:
            commentdict = dict([pair.split('=') for pair in comment.split(' ')])
        elif ' ' in comment:
            key,time = comment.split(' ')
            commentdict = {key:time}
        elif comment.strip() != '':
            try:
                commentdict = {'time':float(comment.strip())}
            except ValueError:
                commentdict = {'comment':comment}
        else:
            commentdict = dict()

        lines = [(i,self.file.readline().strip().replace('\t',' ').split(' ')) for i in range(num)]
        
        firstn, firstline = lines[0]
        len0 = len(firstline)
        badlines = [(n,l) for (n,l) in lines if len(l) not in (4,5,7,8) or len(l) != len0]
        if len(badlines) > 0:
            badn, badline = badlines[0]
            if badline == ['']:
                e = "Coordinate line expected at end of file, line {} (Frame started at {})" .format(
                    self.lineno + badn, self.lineno-2)
            else:
                e = "bad coordinate line {} (Frame started at {}): {}" .format(
                        self.lineno + badn, self.lineno-2, str(badline))
            raise ValueError(e)
            
        self.lineno += num
        
        if len0 == 4:
            lines = [(e,(float(x),float(y),float(z))) for n,(e,x,y,z) in lines]
        elif len0 == 5:
            lines = [(e,(float(x),float(y),float(z)), float(s)) for n,(e,x,y,z,s) in lines]
        elif len0 == 7:
            lines = [(e,(float(x),float(y),float(z)),
                            (float(vx),float(vy),float(vz))) 
                    for n,(e,x,y,z,vx,vy,vz) in lines]
        elif len0 == 8:
            lines = [(e,(float(x),float(y),float(z)),
                            (float(vx),float(vy),float(vz)),
                            float(s)) 
                    for n,(e,x,y,z,vx,vy,vz, s) in lines]
        else:
            raise ValueError
            
        if 'time' in commentdict:
            f = Frame(lines, float(commentdict['time']))
            del commentdict['time']
        else:
            f = Frame(lines, t=None)
        f.comment = comment
        f.values = commentdict
        return f
        
    def __iter__(self):
        location = self.file.tell()
        while True:
            yield self.readframe()
        self.file.seek(location)
    
    def into(self, atoms, check=True):
        for f in self:
            f.into(atoms, check)
            vals = f.values
            if f.sigmas is not None:
                vals = dict(f.values)
                if 'sigmas' in f.values:
                    raise ValueError("sigmas already in f.values, don't know what to do")
                vals['sigmas'] = f.sigmas
            yield f.time, vals
    
    def _convert_time_line(self, line):
        if '=' not in line:
            key,time = line.split(' ')
            return {key:time}
        try:
            return dict([pair.split('=') for pair in line.split(' ')])
        except:
            print('Line:', line, [pair.split('=') for pair in line.split(' ')])
            raise
    
    def vdicts(self):
        #~ print('Getting vdicts...')
        location = self.file.tell()
        self.file.seek(0)
        vlines = [l for l in self.file if 'time' in l]
        vdicts = [self._convert_time_line(l) for l in vlines]
        self.file.seek(location)
        #~ print('Got vdicts.')
        return vdicts
    
    def close(self):
        self.file.close()
    
    def all(self):
        #~ print('Getting frames...')
        frames = Frames(iter(self))
        for n,f in enumerate(frames):
            f.indx = n
        #~ print('Got frames.')
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
        
class Frame:
    dtype = np.float64
    
    def __init__(self, locs, t=None): #, indx=None):
        self.t = t
        #~ self.indx = indx
        resized = list(zip(*locs))
        if len(resized) == 2:
            self.elems, self.locs = resized
            self.vels = None
            self.sigmas = None
        elif len(resized) == 3:
            self.elems, self.locs, vels = resized
            try:
                iter(vels[0])
                self.vels = vels
                self.sigmas = None
            except TypeError:
                self.vels = None
                self.sigmas = vels
        elif len(resized) == 4:
            self.elems, self.locs, self.vels, self.sigmas = resized
                
        else:
            raise ValueError("Resized wrong length (%d)" % len(resized))
        
        if self.vels is not None:
            self.vels = np.array(self.vels, dtype=self.dtype)
        if self.sigmas is not None:
            self.sigmas = np.array(self.sigmas, dtype=self.dtype)
        self.locs = np.array(self.locs, dtype=self.dtype)
    
    def __iter__(self):
        return iter(self.locs)
    
    def __len__(self):
        return len(self.locs)
    
    def _setx(self, atm, loc):
        x,y,z = loc
        atm.x = (float(x), float(y),float(z))
    def _setv(self, atm, vel):
        x,y,z = vel
        atm.v = (float(x), float(y),float(z))
    
    def into(self, atoms, check=True):
        if self.vels is not None:
            for a, elem, loc, vel in zip(atoms, self.elems, self.locs, self.vels):
                if check and hasattr(a, 'element') and a.element != elem:
                    raise TypeError("Element mismatch for Atom %s at %s: %s is not %s" 
                                        % (a.name, loc, a.element, elem))
                self._setx(a, loc)
                self._setv(a, vel)
        else:
            for a, elem, loc in zip(atoms, self.elems, self.locs):
                if check and hasattr(a, 'element') and a.element != elem:
                    raise TypeError("Element mismatch for Atom %s at %s: %s is not %s" 
                                        % (a.name, loc, a.element, elem))
                self._setx(a, loc)
    
    @property
    def time(self):
        return self.t

class Frames(list):
    def into(self, atoms, check=True):
        for f in self:
            f.into(atoms, check)
            yield f.time, f.values

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
        from simw import AtomVec
        f=open('/home/wendell/idp/data/T3-200K.xyz','r')
        x=XYZreader(f)
        #~ print(x.readframe())
        #~ profile.runctx('for l in x: pass', globals(), locals())
        atoms = AtomVec([1]*1013)
        #~ profile.runctx('for t,l in x: Frame(l,t).into(atoms,False)', globals(), locals())
        frames = [Frame(l,t) for t,l in x]
        profile.runctx('for f in frames: f.into(atoms,False)', globals(), locals())
