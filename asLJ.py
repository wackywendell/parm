# encoding: UTF-8
from __future__ import print_function
print("Importing...")

from simw import *
from sim3d import randcol, Window, key, Vidwriter
from math import sqrt
from simpdb import Resvec, make_structure
from Bio.PDB import PDBParser
import numpy as np

tot = 10000
dt = .02 #dt = float(tot) / steps
steps = int(tot / dt)
showsteps = 2000
#~ damping = 0
#~ Temp = 0
bondspring = 100
anglespring = 100
LJsigma = 1
LJepsilon = 1
numres = None
damping = 0.5
Temp = 1
pdbfile = '/home/wendell/idp/pdb/aS.pdb'
loadfile='../blengths/stats.pkl'
#~ moviefile = None
moviefile = '../T1-10k-d.5.mp4'
size = 600
fps = 40

####################
print("parsing...")
pdbp = PDBParser()
aS = pdbp.get_structure('AS', pdbfile)
chain = list(aS.get_chains())[0]

print("Making structure...")
if numres <= 0:
    numres = None
avecs, interactions = make_structure(chain.child_list[:numres], loadfile,
            bondspring, anglespring, LJsigma, LJepsilon)

print(sum([len(av) for av in avecs]), "atoms", len(avecs), "residues")

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)
intervec = ivector(interactions)

t=0
dt = float(tot) / steps
#collectionSol(const flt dt, const flt damping, const flt desiredT, 
                #~ vector<atomgroup*> groups=vector<atomgroup*>(),
                #~ vector<interaction*> interactions=vector<interaction*>());
collec = collectionSol(dt, damping, Temp, atomgroups, intervec)
collec.setForces()
collec.seed()


tlist = []
E=[]
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

print("Running...")
if moviefile:
    vidwriter = Vidwriter(moviefile, size, fps)
    window = Window(size, outfile=vidwriter.inputfile, fps=fps)
else:
    window = Window(size, fps)

def update(actualtimestep):
    global window,t,L, actualtime
    for i in range(steps/showsteps):
        collec.timestep()
        t+=1
    if window.keys[key.BRACKETLEFT]:
        #~ print('left')
        L /= zoom*actualtimestep
    if window.keys[key.BRACKETRIGHT]:
        #~ print('right')
        L *= zoom*actualtimestep
    c = collec.com()
    for av in avecs:
        for a in av:
            n = a.name[:1]
            window.setcol(colors.get(a.name[:1], (1,1,1)))
            loc = (a.x - c)/L + Vec(.5,.5,.5)
            window.drawsphere(tuple(loc), 1.6/float(L))
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
        100.0*t/steps)
    window.caption("E={0:.3f},{1:.3f}%; T={2:.3}/{3:.3},{4:.3}%; {5:.1f}%".format(
        *stats))
    if t >= steps:
        raise StopIteration

begin = int(t * .1)
tlist = tlist[begin:]
T = T[begin:]
E = E[begin:]
Rg = Rg[begin:]

from datetime import datetime
starttime = datetime.now()
window.run(update)
tottime = datetime.now() - starttime

print("%s, %.3f fps" % (tottime, t / tottime.total_seconds()))

if moviefile:
    vidwriter.close()
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
