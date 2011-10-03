from __future__ import print_function
from simw import *

natoms = 20

def makeline(av, spacing=1, svec=Vec(0,0,0)):
    tot = av.N()
    startx = -float(spacing)*(tot-1)/2
    yval = .03
    zval = .03
    for i in range(tot):
        av.get(i).x = svec + Vec(startx + i*spacing,yval,zval)
        yval *= -1
        if(yval < 0): zval *= -1

masses = fvector([1]*natoms)
#~ masses=fvector([1,10000])
av = atomvec((natoms+2)*1.2,masses)
makeline(av)

av2 = atomvec((natoms+2)*1.2,masses)
makeline(av2,1,Vec(0,2,0))

atomgroups = avector([av])

bondinter = spring(100, .9)
bonds = intraMolNNPair(atomgroups,bondinter)
bondanginter = bondangle(100, 2.2)
bondang = intraMolNNTriple(atomgroups,bondanginter)
dihedralinter = dihedral(fvector([1,2.18,.6]))
torsion = intraMolNNQuad(atomgroups,dihedralinter)

ljinter = LJcutoff(30, 3, 7.5)
ljforce = interMolPair(atomgroups,ljinter)

intervec = ivector([bonds,bondang,torsion,ljforce])

collec = collectionSol(0,0,atomgroups, intervec)
#~ collec = collection(atomgroups, intervec)
collec.setForces()
collec.seed()

#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
t=0
tot = 20
steps = 20000
xlocsteps = 50
dt = float(tot) / steps
tlist = []
xlist = []
ylist = []
vlist = []
E=[]
print("running...")
printperc=10
toprint = printperc
for i in range(steps/xlocsteps):
    for i in range(xlocsteps):
        collec.timestep(dt)
        t+=dt
    tlist.append(t)
    xlist.append([a.x.getx() for a in itern(av,av.N())])
    #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
    yval = xlist[-1][0] - xlist[-1][1];
    ylist.append(yval);
    vs = [a.v.mag() for a in itern(av,av.N())]
    vlist.append(vs)
    E.append([collec.Energy(),ljforce.energy()])#,collec.kinetic()])
    if(100.0*t/tot) > toprint:
        print(toprint,'%')
        toprint += printperc
#~ print(xlist[0],xlist[-1])

print("importing...")
import matplotlib.pyplot as plt
import numpy as np
print("plotting...")
print(vlist[-1])
Emean = np.mean(E,0)
Estd = np.std(E)
print("E dev:", Emean, Estd,Estd/Emean)
#~ plt.plot(tlist, xlist,'.-')
#~ plt.title('x')
#~ plt.show()
#~ plt.plot(tlist, ylist,'.-')
#~ plt.title('separation')
#~ plt.show()
plt.plot(tlist, vlist,'.-')
plt.title('velocity')
plt.show()
plt.plot(tlist, E,'.-')
plt.title('Energy')
plt.show()
print("done!")


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
