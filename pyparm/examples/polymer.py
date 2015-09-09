# -*- coding: utf-8 -*-

# A basic example of a polymer simulation.

# Statistics such as energy, (instantaneous) temperature, and (instantaneous) pressure are
# output to "polymer-example.npz", and coordinates are output to "polymer-example.xyz".
# These xyz files can be read by VMD.

# This simulation demonstrates usage of simcli.py, as well as the XYZreader.

# This was written for Python 3.4+; it may take some small effort to use a lower version of Python.

# ------------------------------------------------------------------------------
# std library
import sys
import argparse

# Third-party imports
import numpy as np

# pyparm imports
import pyparm.d3 as sim
import pyparm.util as util
# ------------------------------------------------------------------------------
# Parameters
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

scigroup = parser.add_argument_group('Scientific Parameters')
scigroup.add_argument('-T', '--temp', type=float, default=0.01, help='Simulation temperature')
scigroup.add_argument('-N', '--npoly', type=int, default=20, help='Number of polymers')
scigroup.add_argument('-n', '--perpoly', type=int, default=10,
    help='Number of beads per polymer')
scigroup.add_argument('-d', '--density', type=float, default=0.05,
    help='Number of polymers per unit volume')
scigroup.add_argument('-a', '--angle', type=float, default=120,
    help='Bond angle, in degrees, with 180 being a straight chain')
scigroup.add_argument('-t', '--time', type=float, default=1000.0,
    help='Number of time units to run the simulation')

simgroup = parser.add_argument_group('Simulation Parameters')
simgroup.add_argument('--damping', type=float, default=1.0,
    help='The damping coefficient of the integrator, related to viscosity')
simgroup.add_argument('--dt', type=float, default=1e-2, help='Timestep')
simgroup.add_argument('--mass', type=float, default=1.0, help='Mass per bead')
simgroup.add_argument('--width', type=float, default=1.0, help='Width of beads')
simgroup.add_argument('--anglek', type=float, default=100.0, help='Strength of bond angle')
simgroup.add_argument('--bondk', type=float, default=100.0, help='Strength of bond potential')

outgroup = parser.add_argument_group('Output Parameters')
outgroup.add_argument('--xyz', default='polymer-example.xyz')
outgroup.add_argument('--stats', default='polymer-example.npz')
parser.add_argument('--statt', type=float, default=10.0,
    help="how often to save statistics (time units)")
parser.add_argument('--xyzt', type=float, default=1.0,
    help="how often to output a frame to the XYZ file (time units)")
outgroup.add_argument('--printn', type=int, default=200)

opts = parser.parse_args()

# ------------------------------------------------------------------------------
# Setting up the simulation
L = (opts.npoly / opts.density)**(1.0/3.0)
N = opts.npoly * opts.perpoly
sigmas = [1.]*N

box = sim.OriginBox(L)
atoms = sim.AtomVec([1.]*N)
# the NeighborList, for keeping track of what atoms are near what other atoms
neighbors = sim.NeighborList(box, atoms, 0.4)
repulse = sim.Repulsion(atoms, neighbors)
bonds = sim.BondPairs()
angles = sim.AngleTriples()

# ------------------------------------------------------------------------------
# Initial Conditions
# Now we have created our Atom, but we need to add our atoms to it. We do that
# in a way that prevents overlap

print('Placing atoms...', end=' ')
E0 = 0
angle = (opts.angle % 180) * np.pi/180
lasta = None
preva = None
for n, a, s in zip(range(N), atoms, sigmas):
    E = E0 + 10
    a.v = sim.rand_vec()  # from a gaussian distribution
    repulse.add(sim.EpsSigExpAtom(a, 1.0, s, 2.0))
    if n % opts.perpoly != 0:
        bonds.add(opts.bondk, s, lasta, a)
        neighbors.ignore(lasta, a)
    if n % opts.perpoly > 1:
        angles.add(opts.anglek, angle, preva, lasta, a)
        neighbors.ignore(preva, a)
    while E > E0 + 0.1:
        if n % opts.perpoly > 1:
            dx = sim.rand_vec()
            lastdx = lasta.x - preva.x
            lastdx /= np.linalg.norm(lastdx)
            perp = np.cross(dx, lastdx)
            perp /= np.linalg.norm(perp)
            dx = (lastdx * np.cos(np.pi - angle) +
                  perp * np.sin(np.pi - angle))
            a.x = lasta.x + dx
        elif n % opts.perpoly > 0:
            dx = sim.rand_vec()
            dx /= np.linalg.norm(dx)
            a.x = lasta.x + dx
        else:
            a.x = box.rand_loc()
        neighbors.update_list()
        E = repulse.energy(box) + bonds.energy(box) + angles.energy(box)
    E0 = E
    lastx = a.x
    preva = lasta
    lasta = a
    if n % 10 == 0:
        print(N - n, end=', ')
        sys.stdout.flush()
print('Done.')

# the integrator
# We use a simple velocity-verlet integrator, which is time-reversible and
# NVE ensemble
# i.e., it preserves number of atoms, volume of box, and energy
collec = sim.CollectionSol(box, atoms, opts.dt, opts.damping, opts.temp,
    [repulse, bonds, angles], [neighbors])

# subtract center-of-mass velocity from all particles
collec.reset_com_velocity()
# scale all velocities to get an instantaneous temperature T = T0, at least at the beginning
collec.scale_velocities_to_temp(opts.temp)

################################################################################
# Data Analysis

data_functions = {
    'E': collec.energy,
    'T': collec.temp,
    'U': lambda: repulse.energy(box),
    'P': collec.pressure
}

data_arrays = {k: [] for k in data_functions}
data_arrays['t'] = []


def take_data(time):
    """Take each measurement in data_functions at time 'time', and store it in
    data_arrays"""
    data_arrays['t'].append(time)
    for k, f in data_functions.items():
        data_arrays[k].append(f())


def write_data():
    np.savez_compressed(opts.stats, **data_arrays)

# ------------------------------------------------------------------------------
# XYZ file


def write_xyz(time):
    # print the current line to file
    with open(opts.xyz, 'a') as f:
        print(N, file=f)
        # xyz format for VMD requires a line here, and ignores it; I put the time here.
        print(time, file=f)
        for a in atoms:
            x = box.diff(a.x, sim.vec())
            print('C', *x, file=f)

# empty out the file
with open(opts.xyz, 'w') as f:
    f.truncate()

# ------------------------------------------------------------------------------
# TCL file
# this is just to tell VMD to show the box and show the atoms as the right size.
# Not necessary if you're not using VMD for visualization.

tcl_str = """\
mol modstyle 0 0 VDW 1 32
mol rep VDW 1 20

set cell [pbc set {{{L} {L} {L}}} -all];
pbc box -toggle -center origin -color red;
set natoms [atomselect 0 "name C";];
$natoms set radius {r};""".format(L=L, r=sigmas[0]/2.0)

with open(opts.xyz[:-4] + '.tcl', 'w') as f:
    print(tcl_str, file=f)

xyz_m = 0
data_m = 0
print_m = 1

# for tracking and printing our progress
progress = util.Progress(opts.time)
total_steps = int(opts.time / opts.dt + 0.5)

for step in range(total_steps+1):
    if step > 0:
        collec.timestep()
    time = step*opts.dt
    
    if time > xyz_m * opts.xyzt:
        write_xyz(time)
        xyz_m += 1
        
    if time > data_m * opts.statt:
        take_data(time)
        write_data()
        data_m += 1
    
    if time > print_m * opts.time / opts.printn:
        print(progress.eta_str(time), 'T: %6.2f' % collec.temp())
        print_m += 1
    
print("Done.")
