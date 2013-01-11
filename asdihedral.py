#!/usr/bin/env python
# encoding: UTF-8


from simw import *
from math import sqrt
import simpdb
from xyzfile import XYZwriter, XYZreader
from Bio.PDB import PDBParser
import numpy as np
from os.path import expanduser
from optparse import OptionParser
from collections import defaultdict

mydir = expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1)
parser.add_option('-D', '--damping', type=float, default=.001)
parser.add_option('-t', '--time', type=float, default=1)
parser.add_option('-d', '--dt', type=float, default=.01)
parser.add_option('--rampsteps', type=int, default=0)
parser.add_option('--rampdamp', type=float, default=.5)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=1000)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test/T{T}-{t}K.xyz'))
parser.add_option('-C', '--continue', dest='cont', action='store_true')
parser.add_option('-K', '--chargek', dest='chargek', type=float, 
            default=0)
parser.add_option('--startdt', type=float, default=.000000001)
parser.add_option('--startsteps', type=int, default=0)
parser.add_option('--startdamp', type=float, default=20000)
parser.add_option('-A', '--atomsize', type=float, default=1.0)
parser.add_option('--nodihedral', dest='dihedral', action='store_false')
parser.add_option('-H','--hydrogen', dest='hydrogen', action='store_true')

opts,args = parser.parse_args()

steps = int(opts.time * 1000 / opts.dt + .5)
showtime = opts.time * 1000 / opts.showsteps
showsteps = int(showtime / opts.dt + .5)
bondspring = 5000
anglespring = 200
LJepsilon = 1
chargescreen = 9

# experiments run at 293K
# sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom -> time units
t0 = 6.4069e-14
# then we can get a damping coefficient
# from Allen and Tildesley, \xi = k_B T / m D
# those are mass, D of water, right? 18.01528 g/mol, 2.2 e-9 m^2/s 
# D is from http://hyperphysics.phy-astr.gsu.edu/hbase/tables/liqprop.html
# then we get \xi = 5.727e13 hz
# or 3.94 / t0
# boltzmann⋅293K/((18.01528 g/mol) 2.2e−9 m^2/s)⋅avogadro sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom
damp0 = 3.94
if opts.damping < 0: opts.damping = damp0

# comes from using 300K, angstrom, and e_charge as standard units
# (e_charge^2)/(4pi 80 electric_constant angstrom boltzmann (310.65K))
#                    ^ 80 comes from relative permittivity of water
chargek = 6.954 * opts.chargek
#this also gives 62.22 fs as the base time unit
# sqrt((u/avogadro)/(boltzmann⋅310.65K))⋅angstrom

pdbfile = mydir + 'pdb/aS.pdb'
loadfile= mydir + 'blengths/stats.pkl'
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'), 
                            t=format(opts.time, '.4g'),
                            N=opts.numres)
if opts.cont:
    moviefile = opts.startfile
    
startxyz = XYZreader(open(opts.startfile,'r')) if opts.startfile else None
#~ showsteps = int(float(steps) / opts.showsteps+.5)

#~ print('steps:', steps, opts.showsteps, showsteps, showsteps*opts.dt)

####################
print("Importing PDB...")
avecs = simpdb.Resvec.from_pdb('aS', pdbfile, loadfile, H=opts.hydrogen, numchains=1)

if startxyz is not None:
    print("Importing xyz file", opts.startfile)
    frames = startxyz.all()
    frames[-1].into([a for res in avecs for a in res])
    startedat = frames[-1].time
    startxyz.close()
    del startxyz

if opts.numres > 0:
    oldavecs = avecs
    avecs = avecs[:opts.numres]
    print("%d Residues found, %d used, with %d atoms. Making structure..."
                % (len(oldavecs), len(avecs), sum(len(r) for r in avecs)))
else:
    print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))


bonds = simpdb.make_bonds(avecs, bondspring)
angles = simpdb.make_angles(avecs, anglespring)
LJ,neighbors = simpdb.make_LJ(avecs, LJepsilon, opts.atomsize, 2.0)
dihedrals = simpdb.make_dihedrals(avecs, anglespring) if opts.dihedral else None
charges = simpdb.make_charges(avecs, chargescreen, k=chargek) if opts.chargek else None
interactions = (bonds, angles, LJ, charges, dihedrals)

interactions = ivector([i for i in interactions if i is not None])
if opts.chargek: print('Charged', opts.chargek)
else: print('No charges.')

trackers = tvector([neighbors])

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

t=0
collec = collectionSol(opts.dt, opts.damping, opts.temp, atomgroups, interactions, trackers)
collec.seed()
collec.setForces()

# Stage 1: Run really, really slowly, increasing dt every now and then.
# Only needed occasionally when parameters have changed and forces might 
# be very strong.
factor=3
if opts.startsteps > 0:
    print("Running startsteps...")
    damp, dt = opts.startdamp, opts.startdt
    while dt < opts.dt:
        print('factor', opts.dt / dt)
        collec.changeT(damp, dt)
        for i in range(opts.startsteps):
            collec.timestep()
            if i % (opts.startsteps/5) == 0:
                print('E:', collec.energy())
        dt *= factor

    collec.changeT(opts.damping, opts.dt)
    print("Finished startsteps.")
    
# Stage 2: Run at normal speed for a while to 'ramp up' to the right temperature
# only needed for low damping.
if opts.rampsteps > 0:
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('comv:', oldcomv, '->', collec.comv().mag(),
            'L:', oldL, '->', collec.angmomentum().mag())
    print("Damping for", opts.rampsteps, "steps")
    collec.changeT(opts.rampdamp, opts.temp)
    for i in range(opts.rampsteps):
        collec.timestep()
    collec.changeT(opts.damping, opts.temp)
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('comv:', oldcomv, '->', collec.comv().mag(),
            'L:', oldL, '->', collec.angmomentum().mag())

mode = 'a' if opts.cont else 'w'
xyz = XYZwriter(open(moviefile, mode))
if not opts.cont: xyz.writefull(0, avecs, collec)

valtracker = defaultdict(list)

def printlist(lst, name):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    sovm = 100*std / mean
    print("{name}={mean:.2f}, σ={std:.3f} ({sovm:.2g}%)".format(**locals()))

print("Running... output to " + str(moviefile))

# t is current step num
t = 0 if not opts.cont else int(startedat / opts.dt +.5)
curlim = t + showsteps
print('Starting.', t, curlim)
from datetime import datetime
starttime = datetime.now()
try:
    while t < steps:
        while t < curlim:
            collec.timestep()
            t+=1
            #~ if t % (100) == 0:
                #~ print('t:', t*opts.dt)
        curlim += showsteps
        c = collec.com()
        values = xyz.writefull(int(t * opts.dt+.5), avecs, collec)
        #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
        print('------', int(t*opts.dt+.5))
        for k,v in list(values.items()):
            valtracker[k].append(v)
            if k is 'time': continue
            printlist(valtracker[k], k)
except KeyboardInterrupt:
    pass

tottime = datetime.now() - starttime
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

print("%s, %.3f fps" % (tottime, t / totsecs))
