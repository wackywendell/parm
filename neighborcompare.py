# encoding: UTF-8
from __future__ import print_function
print("Importing...")

from simw import *
from math import sqrt
import simpdb
from xyzfile import XYZwriter, XYZreader
from Bio.PDB import PDBParser
import numpy as np
from os.path import expanduser
from optparse import OptionParser

def average(lst):
    return sum(lst) / len(lst)

def average_squared(lst):
    return sqrt(sum([s*s for s in lst]) / len(lst))

mydir = expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1)
parser.add_option('-D', '--damping', type=float, dest='damping', default=.5)
parser.add_option('--dsteps', type=int, dest='dsteps', default=1000)
parser.add_option('-t', '--time', type=float, dest='time', default=1)
parser.add_option('-d', '--dt', type=float, dest='dt', default=.01)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=1000)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test/T{T}-{t}K-{L}.xyz'))
parser.add_option('-L', '--label', dest = 'label', default=None)
parser.add_option('-X', dest='simplified', action='store_true')
opts,args = parser.parse_args()

steps = int(opts.time * 1000 / opts.dt + .5)
bondspring = 5000
anglespring = 30
LJepsilon = 1
chargescreen = 9 # 9 angstrom Debye length
chargek = 6.954

pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'
label = opts.label if opts.label else (
    'simplified' if opts.simplified else 'neighbors')
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'),
    t=format(opts.time, '.4g'), N=opts.numres, L=label, d=opts.dt)
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

print(len(list(simpdb.Resvec.all_bonds(avecs))), 'bonds, ',
    len(list(simpdb.Resvec.all_angles(avecs))), 'angles')

bonds = simpdb.make_bonds(avecs, bondspring)
angles = simpdb.make_angles(avecs, anglespring)
charges = simpdb.make_charges(avecs, chargescreen, k=chargek)
dihedrals = simpdb.make_dihedrals(avecs, anglespring*10)

if not opts.simplified:
    LJ,neighbors = simpdb.make_LJ(avecs, LJepsilon, 2)
    trackers = tvector([neighbors])
    print('ignoring', LJ.list().ignore_size())
else:
    LJ = simpdb.make_LJ_simple(avecs, 2.5, LJepsilon)
    trackers = tvector([])
    print('ignoring', LJ.ignore_size())
    
interactions = ivector([bonds, angles, LJ, dihedrals])
#~ interactions = ivector([bonds, angles, LJ])

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

t=0
collec = collectionSol(opts.dt, opts.damping, opts.temp, atomgroups, interactions, trackers)
collec.seed(0)
collec.setForces()

tlist = []
E=[]
LJE=[]
K=[]
T=[]
Rg=[]
bstds=[]
astds=[]



print("Running... output to " + str(moviefile))
xyz = XYZwriter(open(moviefile, 'w'))

print('Initialized: LJ_E=%.3f; E=%.3f; cos=%.3f' % (LJ.energy(), collec.energy(), dihedrals.mean_dists()))
if opts.damping == 0:
    print("Damping for", opts.dsteps, "steps")
    collec.changeT(.5, opts.temp)
    for i in range(opts.dsteps):
        collec.timestep()
    collec.changeT(opts.damping, opts.temp)
    

xyz.writeframe(avecs, 'time 0', collec.com())
print('At 0: LJ_E=%.3f; E=%.3f' % (LJ.energy(), collec.energy()))
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
        
        curE, curLJE, curT = collec.energy(), LJ.energy(), collec.temp()
        c = collec.com()
        E.append(curE)
        LJE.append(curLJE)
        T.append(curT)
        comment = 'time %d' % int(t * opts.dt+.5)
        curE = collec.energy()
        Emean = float(np.mean(E))
        Tmean = float(np.mean(T))
        output = 'At %.2f: LJ_E=%.3f; E=%.3f; T=%.2f; cos=%.4f' % (t * opts.dt, curLJE, curE, curT, dihedrals.mean_dists())
        if not opts.simplified:
            updates, npairs = neighbors.which(), neighbors.numpairs()
            output += '; updates: %d (%.2f, %d)' % (updates, float(t)/updates, npairs)
        print(output)
        if opts.damping == 0:
            #~ print("bond std: %.5f; angle std: %.5f" % 
                    #~ (bonds.std_dists(), angles.std_dists()))
            bstds.append(bonds.std_dists())
            astds.append(angles.std_dists())
            print("DeltaE: %.5f; bond std: %.5f; angle std: %.5f; Tavg=%.2f" % 
                    (np.std(E)/np.mean(E), average_squared(bstds), average_squared(astds), Tmean))
        xyz.writeframe(avecs, comment, c)
        tlist.append(t*opts.dt)
        #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
        K.append(collec.kinetic())
        Rg.append(collec.gyradius())
        #~ stats = (curE, float(100*np.std(E))/Emean, curT, Tmean, float(100*np.std(T))/Tmean, 
            #~ 100.0*t/steps)
except KeyboardInterrupt:
    pass
    
if opts.damping == 0:
    print("std(E): %.5g | std/mean: %.4g" % (np.std(E), np.std(E)/np.mean(E)))
    print("bond stds:", average_squared(bstds), "angle stds:", average_squared(astds))

begin = int(t * .1)
tlist = tlist[begin:]
T = T[begin:]
E = E[begin:]
Rg = Rg[begin:]

tottime = datetime.now() - starttime
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

print("%s, %.3f fps" % (tottime, t / totsecs))
print("done!")
