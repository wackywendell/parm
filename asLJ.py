# encoding: UTF-8
from __future__ import print_function
print("Importing...")

from simw import *
from math import sqrt
from simpdb import Resvec, make_structure
from xyzfile import XYZ
from Bio.PDB import PDBParser
import numpy as np

damping = 0.5
Temp = 3

tot = 50000
dt = .02 #dt = float(tot) / steps
#~ tot = steps * dt
steps = int(tot / dt)
showsteps = 10000
bondspring = 1000
anglespring = 1000
LJepsilon = 1
numres = None
#~ damping = 0.5
#~ Temp = 2
pdbfile = '/home/wendell/idp/pdb/aS.pdb'
loadfile='../blengths/stats.pkl'
moviefile = '../T%d-%dK.xyz' % (Temp, int(tot / 1000))
#~ moviefile = '../T10-50K.mp4'
size = 600
fps = 40

####################
print("parsing...")
avecs = Resvec.from_pdb('aS', pdbfile, loadfile, numchains=1)


if numres > 0:
    oldavecs = avecs
    avecs = avecs[:numres]
    print("%d Residues found, %d used, with %d atoms. Making structure..."
                % (len(oldavecs), len(avecs), sum(len(r) for r in avecs)))

print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))


interactions, trackers = make_structure(avecs,
            bondspring, anglespring, LJepsilon)

LJ = interactions[-1]

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)
intervec = ivector(interactions)
trackvec = tvector(trackers)

t=0
collec = collectionSol(dt, damping, Temp, atomgroups, intervec, trackvec)
collec.seed()
collec.setForces()

tlist = []
E=[]
LJE=[]
K=[]
T=[]
Rg=[]
c = collec.com()
L = max((a.x - c).mag() * 2 for av in atomgroups for a in av)
#~ L=numres * 5
zoom = 1.03

# see https://secure.wikimedia.org/wikipedia/en/wiki/CPK_coloring
def f(r,g,b):
    return r/255.,g/255.,b/255.
CPKcolors = {'C':f(32,32,32),'N':f(32,96,255),'O':f(238,32,16),
             'S':f(255,255,48)}
#rasmol
colors = {'C':f(211,211,211),'N':f(135,206,230),'O':f(255,0,0),
             'S':f(255,255,0)}
#from http://life.nthu.edu.tw/~fmhsu/rasframe/COLORS.HTM#cpkcolors
colors = {'C':f(200,200,200),'N':f(143,143,255),'O':f(240,0,0),
             'S':f(255,200,50)}
del f

print("Running... output to " + str(moviefile))
xyz = XYZ(open(moviefile, 'w'))
    
t = 0
from datetime import datetime
starttime = datetime.now()
try:
    for i in range(showsteps):
        for j in range(steps/showsteps):
            collec.timestep()
            t+=1
        
        c = collec.com()
        xyz.writeframe(avecs, 'time %.2f' % (t * dt), c)
        tlist.append(t*dt)
        #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
        curE = collec.Energy()
        E.append(curE)#,collec.kinetic()])
        LJE.append(LJ.energy())
        curT = collec.Temp()
        T.append(curT)
        K.append(collec.kinetic())
        Emean = float(np.mean(E))
        Tmean = float(np.mean(T))
        Rg.append(collec.gyradius())
        #~ stats = (curE, float(100*np.std(E))/Emean, curT, Tmean, float(100*np.std(T))/Tmean, 
            #~ 100.0*t/steps)
except KeyboardInterrupt:
    pass
    
begin = int(t * .1)
tlist = tlist[begin:]
T = T[begin:]
E = E[begin:]
Rg = Rg[begin:]

tottime = datetime.now() - starttime

print("%s, %.3f fps" % (tottime, t / tottime.total_seconds()))

#~ print(xlist[0],xlist[-1])
exit();
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

Pairs = [(E, "E"),(T, "T"),(LJE, "LJ E")]#,(Rg, "Rg"), (K, "K")]

for xs, name in Pairs:
    plot(xs, name)

print("done!")


#~ for a in itern(av, av.N()):
    #~ print(tuple(itern(a.x,3)))
