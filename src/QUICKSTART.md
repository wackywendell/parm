# Quickstart

## C++

Based on \ref LJatoms.cpp

### Setup

Include the necessary header files:

#### Make an `AtomVec`

Create 20 atoms, with mass 1.0:

~~~{.cpp}
uint Natoms = 20;
boost::shared_ptr<AtomVec> atomptr(new AtomVec(Natoms, 1.0));
~~~

We use `boost::shared_ptr` for the major objects, to prevent memory leaks / segmentation faults and
integrate more smoothly with Python. If you are not familiar with `shared_ptr`, its a wrapper around
a pointer such that multiple copies can be made, and the object will be deleted when the last
`shared_ptr` to it is deleted. See [a
reference](http://www.cplusplus.com/reference/memory/shared_ptr/) for more details.

#### Make a `Box`

~~~{.cpp}
const flt L = 20;
boost::shared_ptr<OriginBox> obox(new OriginBox(L));
~~~

Here, we'll use a periodic box with sides of length `L`.

#### Make at least one `Interaction`

We'll use the Lennard-Jones Interaction, \f$V\left(r\right)=\varepsilon\left(1-\frac{\sigma^{6}}{r^{6}}\right)^{2}\f$, where \f$\sigma\f$ is
the size of each particle.

For computational purposes, we "truncate and shift" this potential at \f$2.5\sigma\f$, a standard
practice. We will also use a Neighbor list, which keeps track of which atoms are in the vicinity
of each other.

We create that `Interaction`, and give its `NeighborList` a variable name:

~~~{.cpp}
const flt sigma = 1.;
const flt sigcut = 2.5;
boost::shared_ptr<NListed<EpsSigCutAtom, LennardJonesCutPair> > 
	LJ(new NListed<EpsSigCutAtom, LennardJonesCutPair>(obox, atomptr, 0.4*(sigcut*sigma)));
boost::shared_ptr<NeighborList> neighbor_list = LJ->neighbor_list();
~~~{.cpp}

Note the use of templating: `NListed<A, P>` is a template for an `Interaction` that *can* be 
used with a NeighborList, and can be used with a number of different potentials, or custom ones.

#### Initialize the Atoms

~~~{.cpp}
for (uint i=0; i < atomptr->size(); i++){
	atoms[i].x = obox->rand_loc(); // random location in the box
	atoms[i].v = rand_vec(); // from a Gaussian
	atoms[i].f = Vec();
	atoms[i].a = Vec();
	
	// Add it to the Lennard-Jones potential
	LJ->add(EpsSigCutAtom(epsilon, sigma, atoms.get_id(i), sigcut));
}

// Update the neighbor list
neighbor_list->update_list(true);
~~~

Now we have given every `Atom` a position and a random velocity.

#### Make a `Collection`

A `Collection` is a class that represents an integrator. There are various kinds of integrators:
constant energy, constant temperature, energy minimization, etc. Here, we choose a Velocity Verlet
integrator, which is constant energy:

~~~{.cpp}
const flt dt = 0.0001;
CollectionVerlet collec = CollectionVerlet(boost::static_pointer_cast<Box>(obox), atomptr, dt);
~~~

#### Link the `Interaction`, `NeighborList`, and `Collection`

We made our Collection, but it needs to know what potential to use with the atoms, and that it needs
to update the `NeighborList`:

~~~{.cpp}
collec.add_tracker(neighbor_list);
collec.add_interaction(LJ);
~~~

### Running the simulation

Let's run it for 500,000 timesteps, outputting the current energy, kinetic_energy energy, and potential
energy as we go:

~~~{.cpp}
for(uint i=0; i<500; i++){
	for(uint j=0; j<1000; j++){
		collec.timestep();
	}
	writefile(outfile, atoms, *obox);
	cout << (500-i) << " E: " << collec.energy() << " K: " << collec.kinetic_energy() 
		<< " U: " << LJ->energy(*obox) << "\n";
}
~~~

That's all! To see this example fully fleshed out, see \ref LJatoms.cpp