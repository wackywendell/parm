#!/usr/bin/env python
# encoding: UTF-8


from simw import *
from math import sqrt
import simpdb, statistics, util
import FRETs
from xyzfile import XYZwriter, XYZreader
import sys, os, os.path
from datetime import datetime
import numpy as np
import os.path
from optparse import OptionParser
from collections import defaultdict
from sys import stdout, stderr

mydir = os.path.expanduser('~/idp/')
parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='temp', default=1.0)
parser.add_option('-D', '--damping', type=float, default=.001)
parser.add_option('-t', '--time', type=float, default=1.0)
parser.add_option('-d', '--dt', type=float, default=.01)
parser.add_option('--rampsteps', type=int, default=0)
parser.add_option('--rampdamp', type=float, default=.5)
parser.add_option('--rampdt', type=float, default=None)
parser.add_option('-S', '--showsteps', type=int, dest='showsteps', default=0)
parser.add_option('-Z', '--showsize', type=float, dest='showsize', default=.25)
parser.add_option('-N', '--numres', type=int, dest='numres', default=0)
parser.add_option('-s', '--startfile', dest='startfile', default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=(mydir + 'test/{i}-a{a}-T{T}.xyz'))
parser.add_option('-g', '--statfile', dest='statfile', 
            default=None)
parser.add_option('-C', '--continue', dest='cont', action='store_true')
parser.add_option('-K', '--chargek', dest='chargek', type=float, 
            default=1.0)
parser.add_option('-k', '--keepcharged', dest='subtract', 
                                action='store_false', default=True)
parser.add_option('-a', '--alpha', type=float,default=1.2)
parser.add_option('--seed', dest='seed', action='store_false',
            default=True)
parser.add_option('--relaxdt', type=int, default=12)
parser.add_option('--relaxsteps', type=int, default=0)
parser.add_option('--relaxdamp', type=float, default=200.0)
parser.add_option('-A', '--atomsize', type=float, default=None)
parser.add_option('--nodihedral', dest='dihedral', action='store_false', default=True)
parser.add_option('-H','--hydrogen', dest='hydrogen', action='store_true', default=False)
parser.add_option('--dihk', type=float, default=215.0)
parser.add_option('--hphobsigma', type=float, default=4.8) # monomer distance from random walk (44 A) to CGMD (8.5) comparison is 5.2
parser.add_option('--hmix', type=float, default=0)
parser.add_option('--norij', dest='rij', action='store_false', default=True)
parser.add_option('--hset', type='choice', default='KD',
        choices = ['monera', 'sharpcorrected','KD','aWW','hydroavg'])
parser.add_option('--hscale', type='choice', default='zeroone',
        choices=['None', 'max', 'minmax', 'zeroone'])
parser.add_option('--hjoin', type='choice', default='arithmetic',
        choices=['arithmetic', 'arithmeticzero','geometriczero','max','maxzero'])
parser.add_option( '-i', dest='protein', type=str)#, default='aS', type='choice',
                    #choices=['aS','bS','gS','tau'])

print(" ".join(sys.argv))
opts,args = parser.parse_args()
if opts.atomsize is None:
    opts.atomsize = 1 if opts.hydrogen else 1.0 / (2**(1.0/6.0))

steps = int(opts.time * 1000 / opts.dt + .5)
showtime = opts.time * 1000 / opts.showsteps if opts.showsteps else 1000*opts.showsize
showsteps = int(showtime / opts.dt + .5)
chargescreen = 9

# comes from using 293K, angstrom, and e_charge as standard units
# (e_charge^2)/(4pi 80 electric_constant angstrom boltzmann (293K))
#                    ^ 80 comes from relative permittivity of water
chargek = 7.1288713 * opts.chargek
#this also gives 64.07 fs as the base time unit
# sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom

LJe = chargek / opts.hphobsigma * opts.alpha if opts.alpha > 0 else 0
# LJe = .516 comes from Ashbaugh: 0.300kcal/mol / (avogadro⋅boltzmann⋅293K)

LJepsilon = 1 # for repulsive

if not hasattr(simpdb, opts.hset):
    print("Hydrophobicity scale %s not found" % opts.hset, file=sys.stderr)
    exit(1)


# experiments run at 293K
# sqrt((u/avogadro)/(boltzmann⋅293K))⋅angstrom -> time units
t0 = 6.4069e-14


# Make hydrophobicity coeffs
hset = dict(getattr(simpdb, opts.hset))

#-----------------------------------------------------------------------
# Load PDB
if hasattr(FRETs, opts.protein + 'pdb'):
    pdbfile = getattr(FRETs, opts.protein + 'pdb')
    ijs, _, _ = getattr(FRETs, opts.protein)
else:
    pdbfile = opts.protein
    ijs = []
loadfile= mydir + 'blengths/stats-H.pkl'
protname, _ = os.path.splitext(os.path.basename(opts.protein))
moviefile = opts.xyzfile.format(T=format(opts.temp, '.3g'), 
                            t=format(opts.time, '.4g'),
                            i=protname,
                            a=format(opts.alpha, '.4g'),
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
avecs = simpdb.Resvec.from_pdb('protein', pdbfile, loadfile, H=opts.hydrogen, 
                        numchains=1, numres=numres)

# Import startfile
if not opts.startfile:
    startxyz = None
else:
    print("Reading xyz file", opts.startfile)
    lastframe = None
    try:
        startxyz = XYZreader(open(opts.startfile,'r'))
        while True:
            frame = startxyz.readframe()
            t = float(frame.time)
            if not math.isnan(t) and not math.isinf(t):
                lastframe = frame
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
    if float(startedat) >= opts.time: exit
    startxyz.close()
    del startxyz

print("%d Residues found with %d atoms. Making structure..."
                % (len(avecs), sum(len(r) for r in avecs)))


########################################################################
# Set up interactions
box = InfiniteBox()

angles = simpdb.make_angles(avecs, opts.temp, usestd=True)
bonds = simpdb.make_bonds(avecs, opts.temp, usestd=True)
#~ angles = None
if opts.hydrogen:
    LJ,neighbors = simpdb.make_LJ(box, avecs, LJepsilon, opts.atomsize,
                                            2.0, simpdb.AliceSizes)
    print("made LJ with Alice's atom sizes x %.2f" % opts.atomsize)
else:
    LJ,neighbors = simpdb.make_LJ(box, avecs, LJepsilon, opts.atomsize,
                                            2.0, simpdb.RichardsSizes)
neighbors.update_list()
print("Made LJ, energy", LJ.energy(box))
dihedrals = (simpdb.make_dihedrals(avecs, opts.dihk) 
            if (opts.dihk > 0 and not opts.hydrogen) else None)
charges = simpdb.make_charges(avecs, chargescreen, k=chargek, 
                    subtract=opts.subtract, ph=7.4) if chargek > 0 else None
hphob, HPneighbors = (
    simpdb.make_hydrophobicity(box, avecs, LJe,
                        hydroindex=hset, scale=opts.hscale, method=opts.hjoin,
                        sigma=opts.hphobsigma)
    if LJe > 0
    else (None,None))

interactions = dict(bond=bonds, angle=angles, LJ=LJ, electric=charges, 
                dihedral=dihedrals, HPhob=hphob)
interactions = dict([(k,v) for k,v in list(interactions.items()) if v is not None])

intervec = ivector(list(interactions.values()))
#print(opts.dihedral, opts.hydrogen, (opts.dihedral and not opts.hydrogen))
print("interactions:", ", ".join(list(interactions.keys())))
#~ if opts.chargek: print('Charged', opts.chargek)
#~ else: print('No charges.')

trackers = tvector([neighbors, HPneighbors]) if HPneighbors is not None else tvector([neighbors]) 

if len(avecs) < 20:
    print(", ".join([r.resname for r in avecs]))

atomgroups = avector(avecs)

t=0
#print(opts.dt, opts.damping, opts.temp) #INFO
collec = collectionSol(box, opts.dt, opts.damping, opts.temp, atomgroups, intervec, trackers)
collec.interactions = interactions
if opts.seed: collec.seed()
else: collec.seed(1)
collec.setForces()

mode = 'a' if opts.cont else 'w'
xyz = XYZwriter(open(moviefile, mode), rijpairs = ijs)

########################################################################
# Stage 1: Run really, really slowly, increasing dt every now and then.
# Only needed occasionally when parameters have changed and forces might 
# be very strong.
if opts.relaxsteps > 0 and not opts.cont:
    print("Running startsteps...", 'E:', collec.energy())
    for dtfac in range(opts.relaxdt, -1, -1):
        dt = opts.dt / (2**dtfac)
        print('factor: %d; E: ' % (opts.dt / dt), end='')
        collec.changeT(dt, opts.relaxdamp, opts.temp)
        for i in range(opts.relaxsteps):
            collec.timestep()
            if i % (opts.relaxsteps/5) == 0:
                print('  %.6g' % collec.energy(), end='')
                sys.stdout.flush()
                xyz.writeframe(avecs, collec.com(),
                                stage='relax',
                                dt="%.2g" % dt)
        print()
    
    collec.timestep()
    #~ print('E:', collec.energy()) #INFO
    print("Finished startsteps")#, opts.dt, opts.damping, opts.temp) #INFO
    sys.stdout.flush()
    collec.changeT(opts.dt, opts.damping, opts.temp)

########################################################################
# Stage 2: Run at normal speed for a while to 'ramp up' to the right temperature
# only needed for low damping.
K, E, T = collec.kinetic(), collec.energy(), collec.temp()
lastE = float('inf')
frac = K / E
i = 0
if opts.rampsteps > 0 and not opts.cont:
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('comv: %.3g -> %.3g   L: %.3g -> %.3g' % (oldcomv, collec.comv().mag(),
            oldL, collec.angmomentum().mag()))
    print("Damping for", opts.rampsteps, "step increments")
    collec.changeT((opts.rampdt if opts.rampdt else opts.dt), opts.rampdamp, opts.temp)
    while (i < 5 
            or not .95 < (float(T) / float(opts.temp)) < 1.05
            or lastE > E):
        i+=1
        for j in range(opts.rampsteps):
            collec.timestep()
        lastE = E
        K, E, T = collec.kinetic(), collec.energy(), collec.temp()
        frac = K / E
        print(  "E:", '%10d' % collec.energy(), 
                "K:", '%10d' % collec.kinetic(), 
                'T:', '%10.2f' % collec.temp(),
                'K/E:', '%5.2f' % frac)
        sys.stdout.flush()
        xyz.writeframe(avecs, collec.com(),
                            stage='ramp',
                            damp="%.2g" % opts.rampdamp)
    collec.changeT(opts.dt, opts.damping, opts.temp)
    oldcomv = collec.comv().mag()
    collec.resetcomv()
    oldL = collec.angmomentum().mag()
    collec.resetL()
    print('Ramp finished!', 
            'comv: %.3g -> %.3g   L: %.3g -> %.3g' % (oldcomv,
            collec.comv().mag(), oldL, collec.angmomentum().mag()))
    print("E:", '%.2f' % collec.energy(), 'T:', '%.2f' % collec.temp())
    sys.stdout.flush()
    
if not opts.cont: xyz.writefull(0, avecs, collec)

def printlist(lst, name):
    mean = float(np.mean(lst,0))
    std = np.std(lst)
    sovm = 100*std / mean if mean > 0 else float('nan')
    last = float(lst[-1])
    lstd = 100*(last - mean) / std if std > 0 else float('nan')
    print("{name:10s}={mean:8.2f}, σ={std:8.3f} ({sovm:7.3g}%) [{last:8.2f} ({lstd:7.3g}%)]".format(**locals()))

########################################################################
##
## Starting Actual Run
##
########################################################################

print("Running... output to " + str(moviefile))

print("Data output to " + str(statfile))
if (not ijs or not opts.rij):
    print('Rijs are none', not ijs, not opts.rij)
    rijs = None
else:
    print('Rijs are not none', ijs, not ijs, not opts.rij)
    rijs = (statistics.Rijs, ijs)
    pairstr = ', '.join([str(i) + '-' + str(j) for i,j in ijs])
    print("Writing rijs", pairstr)
table = statistics.StatManager(statfile, collec, [avecs],
            stats=(statistics.COMV, statistics.E, statistics.Kinetic,
                    statistics.L, statistics.Rg, statistics.T),
            groups=(statistics.Dihedrals, statistics.CALocs, rijs), 
            contin=opts.cont)
valtracker = defaultdict(list)

def writestats(t):
    values = xyz.writefull(int(t * opts.dt+.5), avecs, collec)
    table.update(int(t * opts.dt+.5), write=True)
    #~ ylist.append([a.x.gety() for a in itern(av,av.N())])
    
    endtime, interval = util.get_eta(t-startt, steps-startt, starttime)
    totinterval = util.to_dhms(endtime - starttime)
    
    print('------ ', int(t*opts.dt+.5), 
        ' Ends in (', util.interval_str(*interval), ' / ', 
            util.interval_str(*totinterval), ') on ', 
                            endtime.strftime('(%a %b %d, %H:%M:%S)'),
        sep = '')
        
    for k,v in sorted(values.items()):
        valtracker[k].append(v)
        if k is 'time': continue
        printlist(valtracker[k], k)
    sys.stdout.flush()

# t is current step num
startt = t = 0 if not opts.cont else round(startedat / opts.dt)
if not opts.cont: table.update(t, write=True)
curlim = t + showsteps
print('Starting.', t, curlim, 'E: %.2f' % collec.energy())
starttime = datetime.now()
writestats(t)
try:
    while t < steps:
        while t < curlim:
            collec.timestep()
            t+=1
            #~ if t % (100) == 0:
            #~ print('t:', t*opts.dt, collec.energy())
        curlim += showsteps
        writestats(t)
        
except KeyboardInterrupt:
    pass

tottime = datetime.now() - starttime
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

print("%s, %.3f fps" % (tottime, (t-startt) / totsecs))
