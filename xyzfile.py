from __future__ import print_function

class XYZ:
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
    
    def close(self):
        self.file.close()
