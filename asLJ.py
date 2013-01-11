# encoding: UTF-8

print("Importing...")

from simw import *
from math import sqrt
import simpdb
from xyzfile import XYZwriter, XYZreader
from Bio.PDB import PDBParser
import numpy as np
from os.path import expanduser
from optparse import OptionParser

mydir = expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1)
parser.add_option('-D', '--damping', type=float, dest='damping', default=.5)
parser.add_option('-t', '--time', type=float, dest='time', default=1)
parser.add_option('-d', '--dt', type=float, dest='dt', default=.02)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=1000)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test-T{T}-{t}K.xyz'))
opts,args = parser.parse_args()

steps = int(opts.time * 1000 / opts.dt + .5)
bondspring = 1000
anglespring = 1000
LJepsilon = 1
pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'), 
                            t=format(opts.time, '.4g'),
                            N=opts.numres)
startxyz = XYZreader(open(opts.startfile,'r')) if opts.startfile else None
showsteps = int(float(steps) / opts.showsteps+.5)

print('steps:', steps, opts.showsteps, showsteps, showsteps*opts.dt)

####################
print("Importing PDB...")
avecs = simpdb.Resvec.from_pdb('aS', pdbfile, loadfile, numchains=1)

if startxyz is not None:
    print("Importing xyz file", opts.startfile)
    frames = startxyz.all()
    frames[-1].into([a for res in avecs for a in res])

if opts.numres > 0:
    oldavecs = avecs
    avecs = avecs[:opts.numres]
    print("%d Residues found, %d used, with %d atoms. Making structure..."
                % (len(oldavecs), len(avecs), sum(len(r) for r in avecs)))

print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))


bonds = simpdb.make_bonds(avecs, bondspring)
angles = simpdb.make_angles(avecs, anglespring)
LJ,neighbors = simpdb.make_LJ(avecs, 2.5, LJepsilon)
interactions = ivector([bonds, angles, LJ])
trackers = tvector([neighbors])

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

t=0
collec = collectionSol(opts.dt, opts.damping, opts.temp, atomgroups, interactions, trackers)
collec.seed()
collec.setForces()

tlist = []
E=[]
LJE=[]
K=[]
T=[]
Rg=[]

print("Running... output to " + str(moviefile))
xyz = XYZwriter(open(moviefile, 'w'))
    
#~ print('Energy:', collec.Energy())

t = 0
from datetime import datetime
starttime = datetime.now()
try:
    for i in range(opts.showsteps):
        curlim = (i+1) * opts.time * 1000 / opts.showsteps - opts.dt/2
        while t * opts.dt < curlim:
            collec.timestep()
            t+=1
            #~ if t % (100) == 0:
                #~ print('t:', t*opts.dt)
        
        c = collec.com()
        xyz.writeframe(avecs, 'time %d' % int(t * opts.dt+.5), c)
        tlist.append(t*opts.dt)
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
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

print("%s, %.3f fps" % (tottime, t / totsecs))

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
    label = '$\overline {} = {:.3f}$'.format(name,mean)
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
