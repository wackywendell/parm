# encoding: UTF-8
from __future__ import print_function
from simw import *
from sim3d import randcol, Window, key
from math import sqrt

tot = 500
dt = .0003 #dt = float(tot) / steps
steps = int(tot / dt)
showsteps = 200
natoms = 10
damping = 0.5
Temp = 10
#~ damping = .8
#~ Temp = 10

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
av = atomvec(masses)
makeline(av)
av2 = atomvec(masses)
makeline(av2,1,Vec(0,2,0))

atomgroups = avector([av])

bondinter = spring(100, .9)
bonds = intraMolNNPair(atomgroups,bondinter)
bondanginter = bondangle(100, 2.2)
bondang = intraMolNNTriple(atomgroups,bondanginter)
dihedralinter = dihedral(fvector([1,1.18,.6]))
torsion = intraMolNNQuad(atomgroups,dihedralinter)

ljinter = LJcutoff(30, 3, 7.5)
ljforce = interMolPair(atomgroups,ljinter)

intervec = ivector([bonds,bondang])
#~ intervec = ivector([bonds,bondang,torsion,ljforce])


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
t=0
dt = float(tot) / steps
sigma = sqrt(4*Temp*damping/dt)
collec = collectionSol(damping,sigma,atomgroups, intervec)
#~ collec = collection(atomgroups, intervec)
collec.setForces()
collec.seed()


tlist = []
E=[]
K=[]
T=[]
Rg=[]
L=natoms
zoom = 1.1


import numpy as np
window = Window(600)
def update(actualtimestep):
    global window,t,L
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
    for a in av:
        loc = (a.x - c)/L + Vec(.5,.5,.5)
        window.drawsphere(tuple(loc), 1/float(L))
    tlist.append(t*dt)
    #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
    curE = collec.Energy()
    E.append(curE)#,collec.kinetic()])
    curT = collec.Temp()
    T.append(curT)
    K.append(collec.kinetic())
    Emean = float(np.mean(E))
    Tmean = float(np.mean(T))
    Rg.append(collec.gyradius())
    stats = (curE, float(100*np.std(E))/Emean, curT, Tmean, float(100*np.std(T))/Tmean, 
        int(100.0*t/steps))
    window.caption("E={0:.3f},{1:.3f}%; T={2:.3}/{3:.3},{4:.3}%; {5}%".format(
        *stats))
    if t >= steps:
        raise StopIteration

begin = int(t * .1)
tlist = tlist[begin:]
T = T[begin:]
E = E[begin:]
Rg = Rg[begin:]

window.run(update, 20)
#~ print(xlist[0],xlist[-1])
#~ exit();
print("importing...")
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
rc('text', usetex=True)
import numpy as np
print("plotting...")

def plotline(ts, num, *args, **kwargs):
    plt.plot((ts[0], ts[-1]), (num,num), *args, **kwargs)

def plot(lst, name):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    sovm = std / mean
    print("{name}={mean:.2f}, σ={std:.3f} ({sovm:.2g})".format(**locals()))
    plt.plot(tlist, lst, 'b.-')
    label = u'$\overline {} = {:.3f}$'.format(name,mean)
    stdlabel = '$\sigma={0:.3f}\;({1:.4f})$'.format(std, sovm)
    plotline(tlist, mean, 'r-', label=label, linewidth=2)
    plotline(tlist, mean-std, 'r-', label = stdlabel)
    plotline(tlist, mean+std, 'r-')
    plt.title(name)
    plt.legend()
    plt.show()

Pairs = [(E, "E"),(T, "T")]#,(Rg, "Rg"), (K, "K")]

for xs, name in Pairs:
    plot(xs, name)

print("done!")


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
