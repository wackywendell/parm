#!/usr/bin/env python
# encoding: UTF-8
from __future__ import print_function

from simw import *
from math import sqrt
import simpdb, statistics
from xyzfile import XYZwriter, XYZreader
from Bio.PDB import PDBParser
import sys, os, os.path
import numpy as np
from os.path import expanduser
from optparse import OptionParser
from collections import defaultdict

mydir = expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1.0)
parser.add_option('-D', '--damping', type=float, default=.001)
parser.add_option('-t', '--time', type=float, default=1.0)
parser.add_option('-d', '--dt', type=float, default=.01)
parser.add_option('--rampsteps', type=int, default=0)
parser.add_option('--rampdamp', type=float, default=.5)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=0)
parser.add_option('-Z', '--showsize', type=float, dest='showsize', default=.25)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test/T{T}.xyz'))
parser.add_option('-g', '--statfile', dest='statfile', 
            default=None)
parser.add_option('-C', '--continue', dest='cont', action='store_true')
parser.add_option('-K', '--chargek', dest='chargek', type=float, 
            default=1.0)
parser.add_option('-k', '--keepcharged', dest='subtract', 
                                action='store_false', default=True)
parser.add_option('-P', '--hydrophobk', dest='hphobk', type=float,default=2.5)
parser.add_option('--seed', dest='seed', action='store_false',
            default=True)
parser.add_option('--startdt', type=float, default=1e-8)
parser.add_option('--startsteps', type=int, default=0)
parser.add_option('--startdamp', type=float, default=20000.0)
parser.add_option('--startfactor', type=float, default=2.0)
parser.add_option('-A', '--atomsize', type=float, default=1.0)
parser.add_option('--nodihedral', dest='dihedral', action='store_false', default=True)
parser.add_option('-H','--hydrogen', dest='hydrogen', action='store_true', default=False)
parser.add_option('--hphobsigma', type=float, default=5.2)
parser.add_option('--hmix', type=float, default=0)
parser.add_option('-3', '--ph3', dest='ph3', action='store_true', default=False)

parser.add_option('--LJES', dest='LJESratio', type=float, default=None)

opts,args = parser.parse_args()

steps = int(opts.time * 1000 / opts.dt + .5)
showtime = opts.time * 1000 / opts.showsteps if opts.showsteps else 1000*opts.showsize
showsteps = int(showtime / opts.dt + .5)
bondspring = 5000
anglespring = 20000
chargescreen = 9

# comes from using 293K, angstrom, and e_charge as standard units
# (e_charge^2)/(4pi 80 electric_constant angstrom boltzmann (293K))
#                    ^ 80 comes from relative permittivity of water
chargek = 7.1288713 * opts.chargek
#this also gives 64.07 fs as the base time unit
# sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom

LJe = .516 * opts.hphobk if opts.hphobk > 0 else 0
# .516 comes from Ashbaugh: 0.300kcal/mol / (avogadro⋅boltzmann⋅293K)

if opts.LJESratio is not None:
    if opts.LJESratio == float('inf'):
        LJe = .516  # for attractive 
        chargek = 0
    elif opts.LJESratio == 0:
        LJe = 0
        chargek = 7.1288713
    else:
        sqratio = sqrt(opts.LJESratio)
        LJe = .516 * sqratio
        chargek = 7.1288713 / sqratio
LJepsilon = 1 # for repulsive
#Hphobs = [[0,LJe],[LJe,2.5*LJe]]

Hphobsigma = opts.hphobsigma # monomer distance from random walk (44 A) to CGMD (8.5) comparison
hydroindex = simpdb.hydroindex7
if opts.ph3:
    print("Using pH 3")
    hydroindex = simpdb.hydroindex3
if opts.hmix:
    import random
    def newval(v):
        newv = v + ((random.random()-0.5) * 2 * opts.hmix)
        return max([0,min([1,newv])])
    hydroindex = dict([(k, newval(v)) for k,v in hydroindex.items()])
print('New Hydrophobicity indices:')
print(*['%s: %.3f' % (r,v) for r,v in hydroindex.items()], sep=', ')

# experiments run at 293K
# sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom -> time units
t0 = 6.4069e-14

pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats-H.pkl'
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'), 
                            t=format(opts.time, '.4g'),
                            N=opts.numres)

if opts.cont:
    moviefile = opts.startfile
print('outputting every %d steps (%d time units) to %s.' % (showsteps, showtime, moviefile))

if opts.statfile is None:
    statfile = str(moviefile).rpartition('.')[0] + '.tsv.gz'
else:
    statfile = opts.statfile.format(T=format(opts.temp, '.3g'), 
                            t=format(opts.time, '.4g'),
                            N=opts.numres)

    
#~ showsteps = int(float(steps) / opts.showsteps+.5)

#~ print('steps:', steps, opts.showsteps, showsteps, showsteps*opts.dt)

####################
print("Importing PDB...")
numres = opts.numres if opts.numres > 0 else None
avecs = simpdb.Resvec.from_pdb('aS', pdbfile, loadfile, H=opts.hydrogen, 
                        numchains=1, numres=numres)

# Import startfile
if not opts.startfile:
    startxyz = None
else:
    lastframe = None
    try:
        startxyz = XYZreader(open(opts.startfile,'r'))
        while True: lastframe = startxyz.readframe()
    except StopIteration:
        pass
    except IOError as e:
        if e.errno == 2 and opts.cont:
            print("File does not exist, continuing anyway.")
            startxyz = None
            opts.cont = False
        else: raise
    if lastframe is None:
        print("Startfile was empty; using PDB locations.")
        opts.cont = False
        startxyz = None

#### Read startfile into atoms
if startxyz is not None:
    print("Importing xyz file", opts.startfile)
    lastframe.into([a for res in avecs for a in res])
    startedat = lastframe.time
    print("Last time is", startedat)
    startxyz.close()
    del startxyz

print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))

#~ if opts.numres > 0:
    #~ oldavecs = avecs
    #~ avecs = avecs[:opts.numres]
    #~ print("%d Residues found, %d used, with %d atoms. Making structure..."
                #~ % (len(oldavecs), len(avecs), sum(len(r) for r in avecs)))
#~ else:
    #~ print("%d Residues found with %d atoms. Making structure..."
                #~ % (len(avecs), sum(len(r) for r in avecs)))


angles = simpdb.make_angles(avecs, anglespring)
bonds = simpdb.make_bonds(avecs, bondspring)
#~ angles = None
if opts.hydrogen:
    LJ,neighbors = simpdb.make_LJ(avecs, LJepsilon, opts.atomsize,
                                            2.0, simpdb.AliceSizes)
else:
    LJ,neighbors = simpdb.make_LJ(avecs, LJepsilon, opts.atomsize,
                                            2.0, simpdb.RichardsSizes)
dihedrals = (simpdb.make_dihedrals(avecs, anglespring) 
            if (opts.dihedral and not opts.hydrogen) else None)
charges = simpdb.make_charges(avecs, chargescreen, k=chargek, 
                    subtract=opts.subtract, ph=(3 if opts.ph3 else 7.4)) if chargek > 0 else None
hphob, HPneighbors = (
    simpdb.make_hydrophobicity(avecs, LJe, #Hphobs, #instead of LJe
                        sigma=Hphobsigma, hydroindex=hydroindex) 
    if LJe > 0
    else (None,None))
interactions = dict(bond=bonds, angle=angles, LJ=LJ, electric=charges, 
                dihedral=dihedrals, HPhob=hphob)
interactions = dict([(k,v) for k,v in interactions.items() if v is not None])

intervec = ivector(list(interactions.values()))
#print(opts.dihedral, opts.hydrogen, (opts.dihedral and not opts.hydrogen))
print("interactions:", ", ".join(interactions.keys()))
#~ if opts.chargek: print('Charged', opts.chargek)
#~ else: print('No charges.')

trackers = tvector([neighbors, HPneighbors]) if HPneighbors is not None else tvector([neighbors]) 

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

t=0
#print(opts.dt, opts.damping, opts.temp) #INFO
collec = collectionSol(opts.dt, opts.damping, opts.temp, atomgroups, intervec, trackers)
collec.interactions = interactions
if opts.seed: collec.seed()
else: collec.seed(1)
collec.setForces()

# Stage 1: Run really, really slowly, increasing dt every now and then.
# Only needed occasionally when parameters have changed and forces might 
# be very strong.
if opts.startsteps > 0:
    print("Running startsteps...", 'E:', collec.energy())
    damp, dt = opts.startdamp, opts.startdt
    while dt < opts.dt:
        print('factor: %.2f; E: ' % (opts.dt / dt), end='')
        collec.changeT(dt, damp, opts.temp)
        for i in range(opts.startsteps):
            collec.timestep()
            if i % (opts.startsteps/5) == 0:
                print('  %.6g' % collec.energy(), end='')
                sys.stdout.flush()
        print()
        dt *= opts.startfactor
    
    collec.timestep()
    #~ print('E:', collec.energy()) #INFO
    print("Finished startsteps")#, opts.dt, opts.damping, opts.temp) #INFO
    collec.changeT(opts.dt, opts.damping, opts.temp)

# Stage 2: Run at normal speed for a while to 'ramp up' to the right temperature
# only needed for low damping.
if opts.rampsteps > 0:
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('comv: %.3g -> %.3g   L: %.3g -> %.3g' % (oldcomv, collec.comv().mag(),
            oldL, collec.angmomentum().mag()))
    print("Damping for", opts.rampsteps, "steps")
    collec.changeT(opts.dt, opts.rampdamp, opts.temp)
    for i in range(opts.rampsteps):
        collec.timestep()
        if i % (opts.rampsteps / 5) == 0:
            print("E:", collec.energy())
    collec.changeT(opts.dt, opts.damping, opts.temp)
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('comv: %.3g -> %.3g   L: %.3g -> %.3g' % (oldcomv, collec.comv().mag(),
            oldL, collec.angmomentum().mag()))
    
mode = 'a' if opts.cont else 'w'
xyz = XYZwriter(open(moviefile, mode))
if not opts.cont: xyz.writefull(0, avecs, collec)

def printlist(lst, name):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    sovm = 100*std / mean
    last = float(lst[-1])
    lstd = 100*(last - mean) / std
    print("{name:10s}={mean:8.2f}, σ={std:8.3f} ({sovm:7.3g}%) [{last:8.2f} ({lstd:7.3g}%)]".format(**locals()))

print("Running... output to " + str(moviefile))

print("Data output to " + str(statfile))
statistics.Rijs([(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)])
table = statistics.StatWriter(statfile, collec, [avecs], contin=opts.cont)
valtracker = defaultdict(list)

# t is current step num
t = 0 if not opts.cont else int(startedat / opts.dt +.5)
curlim = t + showsteps
print('Starting.', t, curlim, 'E:', collec.energy())
from datetime import datetime
starttime = datetime.now()
try:
    while t < steps:
        while t < curlim:
            collec.timestep()
            t+=1
            #~ if t % (100) == 0:
            #~ print('t:', t*opts.dt, collec.energy())
        curlim += showsteps
        c = collec.com()
        values = xyz.writefull(int(t * opts.dt+.5), avecs, collec)
        table.update(int(t * opts.dt+.5))
        #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
        print('------', int(t*opts.dt+.5))
        for k,v in sorted(values.items()):
            #~ print('k:', k)
            valtracker[k].append(v)
            if k is 'time': continue
            printlist(valtracker[k], k)
        #~ for res in avecs:
            #~ if res.resname != 'GLU': continue
            #~ a1,a2 = res['HE2'], res['OE2']
            #~ print(a1.name, a2.name, (a1.x - a2.x).mag())
except KeyboardInterrupt:
    pass

tottime = datetime.now() - starttime
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

print("%s, %.3f fps" % (tottime, t / totsecs))
