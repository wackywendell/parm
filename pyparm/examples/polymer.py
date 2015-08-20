# A basic example of a polymer simulation.

# Statistics such as energy, (instantaneous) temperature, and (instantaneous) pressure are
# output to "polymer-example.npz", and coordinates are output to "polymer-example.xyz".
# These xyz files can be read by VMD.

# This simulation demonstrates usage of simcli.py, as well as the XYZreader.

# This was written for Python 3.4+; it may take some small effort to use a lower version of Python.

#---------------------------------------------------------------------------------------------------
#std library
import argparse

# Third-party imports
import numpy as np
from numpy import pi, sqrt, array

# pyparm imports
import pyparm.d3 as parm
import pyparm.util as util
from pyparm.xyzfile import XYZwriter, XYZreader
from pyparm.simcli import Simulation,StatSet
#---------------------------------------------------------------------------------------------------
# Parameters
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

scigroup = parser.add_argument_group('Scientific Parameters')
scigroup.add_argument('-T', '--temp', type=float, default=1.0, help='Simulation temperature')
scigroup.add_argument('-N', '--npoly', type=int, default=20, help='Number of polymers')
scigroup.add_argument('-n', '--perpoly', type=int, default=10, 
    help='Number of beads per polymer')
scigroup.add_argument('-d', '--density', type=float, default=0.1, 
    help='Number of polymers per unit volume')
scigroup.add_argument('-a', '--angle', type=float, default=160, 
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
parser.add_argument('--xyzt', type=float, default=10.0,
    help=   "how often to output a frame to the XYZ file (time units)")
outgroup.add_argument('--printn', type=int, default=200)

opts = parser.parse_args()

#---------------------------------------------------------------------------------------------------
# Setting up the simulation
L = (opts.npoly / opts.density)**(1.0/3.0)

box = sim.OriginBox(L)
atoms = sim.AtomVec(opts.npoly * opts.perpoly)
neighbors = sim.NeighborList(box, atoms, 0.4) # the NeighborList, for keeping track of what atoms are near what other atoms
LJ = sim.LJgroup(atoms, neighbors)
collec = sim.CollectionVerlet(box, atoms, dt, [LJ], [neighbors]) # the integrator
# We use a simple velocity-verlet integrator, which is time-reversible and NVE ensemble
# i.e., it preserves number of atoms, volume of box, and energy

#---------------------------------------------------------------------------------------------------
# Initial Conditions
# Now we have created our Atom, but we need to add our atoms to it. We do that in a way that prevents overlap

E0 = 0
for a,s in zip(atoms, sigmas):
    E = E0 + 10
    a.v = sim.randVec() # from a gaussian distribution
    LJ.add(sim.LJatom(1, s, a))
    while E > E0 + 0.1:
        a.x = box.randLoc()
        neighbors.update_list()
        E = LJ.energy(box)
    E0 = E

collec.resetcomv() # subtract center-of-mass velocity from all particles
collec.scaleVelocitiesT(T0) # scale all velocities to get an instantaneous temperature T = T0, at least at the beginning

####################################################################################################
# Data Analysis

data_functions = {
    'E' : collec.energy,
    'T' : collec.temp,
    'U' : lambda: LJ.energy(box),
    'P' : collec.pressure
    }

data_arrays = {k : np.zeros((data_n,)) for k in data_functions}
data_arrays['t'] = np.zeros((data_n,))

def take_data(idx, time):
    """Take each measurement in data_functions at time 'time', and store it in data_arrays at index idx"""
    data_arrays['t'][idx] = time
    for k, f in data_functions.items():
        data_arrays[k][idx] = f()

def write_data(idx):
    delim = '\t'
    keys = sorted(data_arrays.keys())
    header = delim.join(keys)
    data = {k:data_arrays[k][:idx+1] for k in keys}
    np.savez_compressed(data_file, **data)

#---------------------------------------------------------------------------------------------------
# XYZ file

element_names = ['C']*N1 + ['O']*N2
def write_xyz(time):
    # print the current line to file
    with open(xyz_file, 'a') as f:
        print(N, file=f)
        print(time, file=f) # xyz format for VMD requires a line here, and ignores it; I put the time here.
        for e,a in zip(element_names, atoms):
            x = box.diff(a.x, sim.Vec())
            print(e, *x, file=f)

# empty out the file
with open(xyz_file, 'w') as f:
    f.truncate()

#---------------------------------------------------------------------------------------------------
# TCL file
# this is just to tell VMD to show the box and show the atoms as the right size.
# Not necessary if you're not using VMD for visualization.

tcl_str = """\
mol modstyle 0 0 VDW 1 32
mol rep VDW 1 20

set cell [pbc set {{{L} {L} {L}}} -all];
pbc box -toggle -center origin -color red;
set natoms [atomselect 0 "name C";];
$natoms set radius {r1};
set natoms [atomselect 0 "name O";];
$natoms set radius {r2};""".format(L=L, r1=sigmas[0]/2.0, r2=sigmas[-1]/2.0)

with open(tcl_file, 'w') as f:
    print(tcl_str, file=f)

xyz_m = 0
data_m = 0
print_m = 1

progress = util.Progress(total_time) # for tracking and printing our progress

for step in range(total_steps+1):
    if step > 0: collec.timestep()
    time = step*dt
    
    if time > xyz_m * xyz_dt:
        write_xyz(time)
        xyz_m += 1        
        
    if time > data_m * data_dt:
        take_data(data_m, time)
        write_data(data_m)
        data_m += 1
    
    if time > print_m * print_dt:
        print(progress.eta_str(time))
        print_m += 1
    
print("Done.")
