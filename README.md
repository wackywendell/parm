ParM
====

#### The PARticle Molecular dynamics library

Purpose
----

The purpose of this library is to provide methods to produce *any* MD
simulation, with any type of integrator, any type of particles, and any
type of force in 2D or 3D.

It is currently single-core only, and provides objects for many
different integrators, and many different forces, which can be composed
together for simulations from a few particles bouncing around in a box
to a full protein simulation.

Examples
----

### C++

There are several examples in C++ in the `src/bin` folder, all with good comments:

* `LJatoms.cpp` (\ref LJatoms.cpp): a simulation of a  Lennard-Jones simulation (2D, or 3D)
* `packer.cpp` (\ref packer.cpp): Generates packings (2D or 3D)
* `hardspheres.cpp` (\ref hardspheres.cpp): for a more unconventional example of collision-driven dynamics

### Python

See `pyparm/examples/LJ.py` (\ref LJ.py) for an example of a simple Lennard-Jones
simulation, with data analysis included.

See `pyparm/packmin.py` for an example of how to make a packing.

Basic Concepts
----

 * `Vec`: this is a "vector" in the physics sense, having either 2 or 3
 dimensions.

 *  `Atom`: an `Atom` is the basic unit of the simulation; it
 represents a particle with mass,
    position, velocity, etc.

 *  `AtomGroup`: an `AtomGroup` is a set of atoms, grouped together for
 the sake of utility.
    * Note that `AtomVec` is a concrete type, and `AtomGroup` is an
    abstract class

 *  `Interaction`: an `Interaction` is a definition of a force and
 energy, such as a Lennard-Jones potential (use
 `NListed<IEpsSigCutAtom, LJAttractRepulsePair>`), springs for
 bonds (`BondPairs`), etc.

    * Note that the neighbor list has been "abstracted" to work with
    many potentials; to use it, you create a `NeighborList`, then use
    `NListed<FooAtom, FooPair>` as the Interaction

 *  `Box`: a box is either infinite (`InfiniteBox`) or periodic
 (`OriginBox`), and takes care of the boundary conditions

 *  `Collection`: a grouping together of a `Box`, `AtomGroup`, and
 `Interaction`s, with an integrator (such as velocity Verlet,
 `CollectionVerlet`, or brownian motion, `CollectionSol`).

Standard Steps
----

1.  Make a `Box`. `OriginBox` is a good standard choice for periodic boundary conditions.
2.  Make an `AtomVec`.
3.  Set masses, positions and velocities
    1. Positions might be set using `box.rand_loc()`, for a random point inside the box
    2. Velocities may be set using `rand_vec()`, which generates a Gaussian distribution, like the expected Boltzmann distribution
3.  Make `Interaction`s. For neighborlisted interactions, make
`NeighborList` first, then make interactions.
    1. `NListed<EpsSigAtom, LJRepulsivePair>` (`LJRepulsive` in Python) is a repulsive Lennard-Jones Interaction
    2. `NListed<EpsSigExpAtom, RepulsionPair>` (`Repulsion` in Python) is a Repulsion or Harmonic Interaction (exponent can be chosen)
    3.   Add atoms / pairs to Interaction
4.  Make a `Collection`. Note that the NeighborList has to be added to
`trackers`
    1. `CollectionVerlet` is a good NVE Collection
5.  Run `Collection.timestep()` many, many times
    1.   Use methods such as `Collection.kinetic_energy()` or
    `Collection.temp()` to get statistics
    2. Or use `tracker`s like `RsqTracker` to track running statistics
6.  Write output to files

Dependencies
----

 - STL, the C++ Standard Template Library
 - Boost: primarily for random numbers, also a few odds and ends
 - (optional) SWIG: for generating Python bindings
 - (optional) Python: for generating Python bindings
    - Known to compile for python 3.2-3.4, and probably with 2.6-2.7

Python
----

This library includes a `sim.i` file for use with SWIG for generating
Python bindings.

#### To Generate Python module

Run `make 2d` or `make 3d` to generate bindings for a 2D or 3D library;
run `make pyparm` to generate bindings for both.

#### Using the Python module

Use `import pyparm.d2 as sim` or `import pyparm.d3 as sim` to import
the module. Then use it freely.

Other Notes
----

### Lennard-Jones

This module uses the equation
\f$V\left(r\right)=\varepsilon\left(1-\frac{\sigma^{6}}{r^{6}}\right)^{2}\f$

The other standard form is
\f$V\left(r\right)=4\varepsilon\left(\frac{\sigma^{\prime12}}{r^{12}}-\frac{\sigma^{12}}{r^{6}}\right)\f$

To convert, use \f$\sigma=2^{\frac{1}{6}}\sigma^{\prime}\f$.
