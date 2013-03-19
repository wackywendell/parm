#!/usr/bin/env python
# encoding: UTF-8

from optparse import OptionParser
from itertools import count
import sys, datetime, collections
import yaml

from simw import *
import simpdb, xyzfile

def errprint(*args, **kw):
    kwargs = {'file':sys.stderr}
    kwargs.update(kw)
    print(*args, **kwargs)
    kwargs['file'].flush()

parser = OptionParser("Runs a simple simulation.")
parser.add_option('-T', '--temperature', type=float, dest='T', default=1.0)
parser.add_option('-D', '--damping', type=float, default=.001)
parser.add_option('-t', '--time', type=float, default=1.0)
parser.add_option('-d', '--dt', type=float, default=.01)

parser.add_option('--rampsteps', type=int, default=0)
parser.add_option('--rampdamp', type=float, default=.5)
parser.add_option('--rampdt', type=float, default=0.0)
                
parser.add_option('-S', '--showcount', type=int, default=0, 
                    help='Number of output lines (overrides showsteps)')
parser.add_option('-Z', '--showsteps', type=float, default=.25,
                    help='Time between output lines (in thousands)')
parser.add_option('-g', '--statfile', dest='statfile', 
            default=None)
parser.add_option('-x', '--xyzfile', dest='xyzfile', 
            default=None)
parser.add_option('-a', '--alpha', dest='alpha', type=float, default=1.0)
parser.add_option('-R', '--Rijs', type='choice', choices=['none','aS','tau'])

errprint(" ".join(sys.argv))
opts, args = parser.parse_args()

if opts.Rijs == 'aS':
    from FRETvals import ijs as FRETijs
elif opts.Rijs == 'tau':
    from tauFRET import ijs as FRETijs
else:
    FRETijs = []
if FRETijs:
    ijstr = ', '.join(['%d-%d' % ij for ij in FRETijs])
    print("Tracking ijs", ijstr)

steps = int(opts.time * 1000 / opts.dt + .5)
showtime = opts.time * 1000 / opts.showcount if opts.showcount else 1000*opts.showsteps
showsteps = int(showtime / opts.dt + .5)

assert len(args) == 1
arg, = args

with open(arg,'r') as f:
    parameters = yaml.load(f)

trackers = []
interactions = {}
constraints = []

########################################################################
## Prep long-range interactions and make atoms
rvecs = []

#--------------------
# CA-CA LJ
sigmas = [atomdict['LJAttractFixedRepulse']['sigma'] * atomdict['LJAttractFixedRepulse']['cut']
            for resdict in parameters['structure']
            for atomdict in resdict['atoms'].values() if 'LJAttractFixedRepulse' in atomdict
        ]
if sigmas:
    innerradius = maxsigmacut = max(sigmas)
    outerradius = innerradius * 2.5
    neighbors = neighborlist(innerradius, outerradius)
    trackers.append(neighbors)
    interactions['LJ'] = LJ = LJAttractFixedRepulse(neighbors)

#--------------------
# Charges
if 'basics' in parameters and 'Charges' in parameters['basics']:
    esparams = parameters['basics']['Charges']
    interactions['Charges'] = charges = Charges(esparams['screeninglength'], esparams['epsilon'])

#-----------------------------------------------------------------------
# Make Atoms
for n,resdict in enumerate(parameters['structure']):
    name = resdict['type']
    atoms = resdict['atoms']
    anames, masses = zip(*sorted([(k, v['mass']) for k,v in atoms.items()]))
    
    res = Res(anames, masses)
    
    for atom in res:
        adict = atoms[atom.name]
        x,y,z = adict['xyz']
        atom.x = Vec(x,y,z)
        #~ print(name, atom.name, adict.keys())
        if 'LJAttractFixedRepulse' in adict.keys():
            params = adict['LJAttractFixedRepulse']
            neweps = [e*opts.alpha for e in params['epsilons']]
            args = [neweps] + [params[k] for k in ['repeps','sigma','indx', 'cut']]
            #~ print('LJish args:',*args)
            LJat = LJAttractFixedRepulseAtom(atom, *args)
            #~ print('LJishAtom sigma:', LJat.sigma)
            LJ.add(LJat)
        if 'Charge' in adict.keys():
            charges.add(atom, adict['Charge'])
    
    rvecs.append(res)

#--------------------
## print out

errprint('Added:', ', '.join(["%d to %s" % (i.size(),n) for n,i in interactions.items()]))

#-----------------------------------------------------------------------
## Make bonds and angles
interactions['bond'] = bonds = bondpairs()
totignored = 0
for n1,aname1,n2,aname2,x0,k in parameters['bonds']:
    try:
        a1,a2 = rvecs[n1][aname1], rvecs[n2][aname2]
    except KeyError:
        print(n1, aname1, n2, aname2)
        raise
    #k = max((k, 200.0))
    bonds.add(k,x0,a1,a2)
    neighbors.ignore(a1,a2)
    totignored += 1
    
errprint("Made", len(parameters['bonds']), "bonds,", end=' ')

interactions['angles'] = angles = angletriples()
for n1,aname1,n2,aname2,n3,aname3,q0,k in parameters['angles']:
    a1,a2,a3 = rvecs[n1][aname1], rvecs[n2][aname2], rvecs[n3][aname3]
    angles.add(k,q0,a1,a2,a3)
    neighbors.ignore(a1,a3)
    totignored += 1
    
errprint(len(parameters['angles']), "angles,", end=' ')

interactions['dihedrals'] = dihedral = dihedrals()
for n1,aname1,n2,aname2,n3,aname3,n4,aname4,coscoeff,sincoeff,usepow in parameters['dihedrals']:
    a1,a2 = rvecs[n1][aname1], rvecs[n2][aname2]
    a3,a4 = rvecs[n3][aname3], rvecs[n4][aname4]
    coscoeff = [opts.T*n for n in coscoeff]
    sincoeff = [opts.T*n for n in sincoeff]
    dihedral.add(coscoeff,sincoeff,a1,a2,a3,a4,usepow)
    neighbors.ignore(a1,a4) # CAᵢ - CAᵢ₊₃
    totignored += 1

errprint(len(parameters['dihedrals']), "dihedrals,", end=' ')


#-----------------------------------------------------------------------
## Create simulator
errprint("Ignoring", neighbors.ignore_size(), "CA-CA pairs.")
atomgroups = [r.atomvec for r in rvecs]
collec = collectionSol(opts.dt, opts.damping, opts.T,
                            atomgroups, list(interactions.values()),
                            trackers, constraints)
collec.seed()
collec.setForces()

########################################################################
## Run simulation

def statistics(t, c=collec):
    d = dict(
        time=t,
        E=collec.energy(),
        KE=collec.kinetic(),
        T=collec.temp(),
        #Distance=(rvecs[0]['CA'].x - rvecs[1]['CA'].x).mag() - x0
        RG=collec.gyradius()
        )
    for name, i in interactions.items():
        d['V' + name] = i.energy()
    return d

#-----------------------------------------------------------------------
# Ramping
if opts.rampsteps > 0:
    collec.resetcomv()
    collec.resetL()
    errprint("Damping for", opts.rampsteps, "steps:\n    E: ", end='')
    collec.changeT(opts.rampdt or opts.dt, opts.rampdamp, opts.T)
    for i in range(opts.rampsteps):
        collec.timestep()
        if i % (opts.rampsteps / 5) == 0:
            errprint("%.3f.." % collec.energy(), end='')
    collec.changeT(opts.dt, opts.damping, opts.T)
    collec.resetcomv()
    collec.resetL()
    errprint("Done.")

startt = t = 0
curlim = t + showsteps
#errprint('Starting', t, curlim, 'E: %8.3f' % collec.energy())
starttime = datetime.datetime.now()
valtracker = collections.defaultdict(list)

if opts.xyzfile:
    global xyz
    #~ print(' '.join(rvecs[0].names), '::', ' ' .join(rvecs[1].names))
    xyz = xyzfile.XYZwriter(open(opts.xyzfile, 'w'), usevels=False, printnames=True)

def statprintline(stats):
    global opts, xyz, collec
    if opts.statfile:
        allstats = Rijs()
        allstats.update(stats)
        ks, vs = zip(*sorted(allstats.items()))
        with open(opts.statfile, 'a') as f:
            errprint(*['%.6g' % v for v in vs], file=f, sep='\t')
    if opts.xyzfile:
        xyz.writeframe(rvecs, com=collec.com(), **stats)

def Rijs():
    return dict([('r' + str(i) + '-' + str(j), 
                        (rvecs[i]['CA'].x - rvecs[j]['CA'].x).mag())
                for i,j in FRETijs])
        

if opts.statfile:
    stats = statistics(0)
    allstats = Rijs()
    allstats.update(stats)
    with open(opts.statfile, 'w') as f:
        ks, vs = zip(*sorted(allstats.items()))
        errprint(*ks, file=f, sep='\t')
    statprintline(stats)

#-----------------------------------------------------------------------
## Actual run
try:
    while t < steps:
        while t < curlim:
            collec.timestep()
            t+=1
        curlim += showsteps
        c = collec.com()
        #~ values = xyz.writefull(int(t * opts.dt+.5), avecs, collec)
        
        stats = statistics(t*opts.dt)
        statprintline(stats)
        #~ table.update(int(t * opts.dt+.5), write=True)
        
        endtime, interval = simpdb.get_eta(t, steps, starttime)
        totinterval = simpdb.to_dhms(endtime - starttime)
        
        #days, hr, mn, sec = interval
        
        errprint('------ ', int(t*opts.dt+.5), 
            ' Ends in (', simpdb.interval_str(*interval), ' / ', 
                simpdb.interval_str(*totinterval), ') on ', 
                                endtime.strftime('(%a %b %d, %H:%M:%S)'),
            sep = '')
        for k,v in sorted(stats.items()):
            valtracker[k].append(v)
            if k is 'time': continue
            simpdb.printlist(valtracker[k], k, file=sys.stderr)
        sys.stderr.flush()
except KeyboardInterrupt:
    errprint('Keyboard Interrupt, stopped.')

#-----------------------------------------------------------------------
## Finish
tottime = datetime.datetime.now() - starttime
totsecs = 24*60*60 * tottime.days + tottime.seconds + (tottime.microseconds / 1000000.)

errprint("%s, %.3f fps" % (tottime, (t-startt) / totsecs))

