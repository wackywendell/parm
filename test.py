# encoding: UTF-8
from __future__ import print_function
from simw import *
from sim3d import randcol, Window, key
from math import sqrt

tot = 200
dt = .0001 #dt = float(tot) / steps
steps = int(tot / dt)
showsteps = 200 #steps/200
natoms = 10
damping = 0
Temp = 0
fps = 40

def rotate(x, y):
    while True:
        yield (x,0)
        yield (0,y)
        yield (-x,0)
        yield (0,-y)

def makeline(av, spacing=1, svec=Vec(0,0,0)):
    tot = av.N()
    startx = -float(spacing)*(tot-1)/2
    startvec = svec + Vec(startx,0,0)
    yval = .3
    zval = .3
    xspace = sqrt(spacing**2  - (yval)**2 - (zval)**2)
    for i, vs in zip(range(0,tot), rotate(yval, zval)):
        y,z = vs
        av.get(i).x = startvec + Vec(i*xspace,y,z)

spacing = 5
masses = fvector([1,10]*(natoms/2))
#~ masses=fvector([1,10000])
av = atomvec(masses)
makeline(av,spacing,Vec(0,0,0))
for i in range(1,av.N()):
    r = av.get(i).x - av.get(i-1).x
    print(r, r.mag())

av2 = atomvec(masses)
makeline(av2,1,Vec(0,2,0))

atomgroups = avector([av])

bondinter = spring(300, spacing)
bonds = intraMolNNPair(atomgroups,bondinter)
bondanginter = bondangle(300, 1.5708)
bondang = intraMolNNTriple(atomgroups,bondanginter)
dihedralinter = dihedral(fvector([30,-60]))
torsion = intraMolNNQuad(atomgroups,dihedralinter)

ljinter = LJcutoff(30, 3, 7.5)
ljforce = interMolPair(atomgroups,ljinter)
ljintraaction = LJcutoff(100, 1.5*spacing, 7.5)
#~ ljintraaction = LJforce(1000, .890899*5)
ljintra = intraMolPairs(atomgroups,ljintraaction,3)
#~ ljintra = intraMolNNPair(atomgroups,ljintraaction)

intervec = ivector([bonds, bondang,torsion,ljintra])
#~ intervec = ivector([bonds,bondang,torsion,ljforce,ljintra])


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
t=0
sigma = sqrt(4*Temp*damping/dt)
collec = collectionSol(damping,sigma,atomgroups, intervec)
#~ collec = collection(atomgroups, intervec)
collec.setForces()
collec.seed()

E = Stat("E",collec.Energy, "E")
Elj = Stat("L-J E",ljintra.energy, "E_{\mathrm{L-J}}")
Ebond = Stat("Bond E", bonds.energy, "E_{\mathrm{bonds}}")
Eang = Stat("Angle E", bondang.energy, "E_{\mathrm{angle}}")
Etor = Stat("Torsion E", torsion.energy, "E_{\mathrm{torsion}}")
K = Stat("K", collec.kinetic, "K")
#~ Emod = Stat("Mod E",lambda: collec.Energy() + ljintra.energy(), "E_{\mathrm{mod}}")
#~ Emod = Stat("Mod E",lambda: (ljintra.energy() + bonds.energy() + 
            #~ bondang.energy() + torsion.energy() + collec.kinetic())
                        #~ , "E_{\mathrm{mod}}")
#~ stats = [(E,Elj,Emod)]
#~ stats = [(E,Ebond,Eang,Etor,Elj,K)]
stats = [(E,Elj)]
L=natoms*spacing
zoom = 1.1


import numpy as np
window = Window(600)
def update(actualtimestep):
    global window,t,L,stats
    for i in range(steps/showsteps):
        collec.timestep(dt)
        t+=1
    if window.keys[key.BRACKETLEFT]:
        #~ print('left')
        L /= zoom*actualtimestep
    if window.keys[key.BRACKETRIGHT]:
        #~ print('right')
        L *= zoom*actualtimestep
    c = collec.com()
    locs = [(a.x - c)/L + Vec(.5,.5,.5) for a in av]
    window.setcol((1,1,1),1)
    window.drawline(locs, thickness=6)
    window.setcol((0,.6,0),.8)
    for loc in locs:
        window.drawsphere(tuple(loc), spacing/float(L))
    for tup in stats:
        for s in tup:
            s.update(t*dt)
    #~ print(t, ljintra.energy())
    #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
    window.caption("{:.1f}% :: ".format(100.0*t/steps) + 
                "; ".join(str(s) for s in tup for tup in stats))
    if t >= steps:
        raise StopIteration

window.run(update, fps)
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
import numpy as np
print("plotting...")

for tup in stats:
    Stat.plots(*tup)

print("done!")


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
