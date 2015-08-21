#include <iostream>
#include <fstream>

#include "vecrand.hpp"
#include "interaction.hpp"
#include "collection.hpp"

// Some constants
// Note that "flt" is a "floating point number", defaults to "double"
// but can be compiled as "long double" if needed
const flt sigma = 1.0;
const flt sigmal = 1.4;
const flt epsilon = 1.0;
const uint Ns = 40;
const uint Nl = 40;
const flt phi0 = 0.001;

// Some algorithmic constants
const flt dt = 0.1;
// the final pressure
const flt P0 = 1e-8;
// The initial pressure
const flt startP = 1e-4;

// The algorithm works by minimizing H = U(x_i...) + PV w.r.t. x_i and κ=log V.
// The final overlap is thus approximately related to P, and the higher the P, the
// faster it converges.
// In this simulation, we start with a high P, and then relax to lower Ps.

// These are the "quitting" parameters, see below for use
const flt P_frac = 1e-4;
#ifndef LONGFLOAT
const flt force_max = 1e-14;
#else
const flt force_max = 1e-18;
#endif

void writefile(AtomVec& atoms, OriginBox& obox);

int main(){
    cout << "Float size: " << sizeof(flt) << " epsilon: " << std::numeric_limits<flt>::epsilon() << "\n";
    assert(std::numeric_limits<flt>::epsilon() < force_max*10);

    // new random seed each time, for velocities and placement
    seed();
    
    // Volume of the spheres
    // Each sphere is σ^d * π/(2d), i.e. π σ^2/4 for 2D, π σ^3/6 for 3D
    // Total volume of the spheres takes a constant N out front
    const flt Vs = (Ns * pow(sigma, NDIM) + Nl * pow(sigmal, NDIM)) * M_PI_2 / NDIM;
    
    // Initial length of the box is [(volume of spheres) / phi0]^(1/d)
    const flt L = pow(Vs / phi0, OVERNDIM);
    cout << "Using L = " << L << "\n";
    
    // Create the bounding box (sides of length L), and a "vector" of Natoms atoms
    boost::shared_ptr<OriginBox> obox(new OriginBox(L));

    // Create a vector of the masses of the atoms
    // We just use all mass 1, because this is a packing
    boost::shared_ptr<AtomVec> atomptr(new AtomVec(Nl + Ns, 1.0));
    AtomVec & atoms = *atomptr;
    
    // Harmonic Interaction
    // Its called "Repulsion" for historical reasons
    // It takes a pointer to the box, a pointer to the atoms, and a "skin depth" for the NeighborList
    boost::shared_ptr<NListed<EpsSigExpAtom, RepulsionPair> > 
    boost::shared_ptr<NListed<EpsSigExpAtom, RepulsionPair> >
    boost::shared_ptr<NListed<EpsSigExpAtom, RepulsionPair> >
        hertzian(new NListed<EpsSigExpAtom, RepulsionPair>(obox, atomptr, 0.1*sigma));
    boost::shared_ptr<NeighborList> nl = hertzian->neighbor_list();
    // ^ this is the Interaction
    
    // Note that NListed is a class template; its an Interaction that
    // takes various structs as template parameters, and turns them into
    // a neighbor Interaction
    // See interaction.hpp for a whole bunch of them
    // Also note that the NeighborList is not the same as the
    // "neighborlisted" Interaction; multiple interactions can use the
    // same NeighborList
    
    // Now we run through all the atoms, set their positions / velocities,
    // and add them to the Repulsion Interaction (i.e., to the neighbor list)
    for (uint i=0; i < atoms.size(); i++){
        atoms[i].x = obox->rand_loc(); // random location in the box
        atoms[i].v = Vec::Zero(); // A zero-vector
        atoms[i].f = Vec::Zero();
        atoms[i].a = Vec::Zero();

        flt sig = i < Ns ? sigma : sigmal;
        
        // Add it to the Repulsion potential
        hertzian->add(EpsSigExpAtom(atoms.get_id(i), epsilon, sig, 2));
        //                                                          ^ exponent for the harmonic Interaction
    }

    // force an update the NeighborList, so we can get an accurate energy
    nl->update_list(true);

    cout << "Starting. Neighborlist contains " << nl->numpairs() << " / " <<
        (atoms.size()*(atoms.size()-1))/2 << " pairs\n";
    
    //Now we make our "Collection"
    CollectionNLCG collec = CollectionNLCG(obox, atomptr, dt, P0);
    
    // This is very important! Otherwise the NeighborList won't update!
    collec.add_tracker(nl);
    // And add the Interaction
    collec.add_interaction(hertzian);
    
    writefile(atoms, *obox);
    
    //Print out total energy, kinetic energy, and potential energy
    cout << "H: " << collec.hamiltonian() << " K: " << collec.kinetic()
        << " U: " << hertzian->energy(*obox) << " phi: " << (Vs/obox->V()) << "\n";
    
    // Run the simulation! And every _ steps, write a frame to the .xyz
    // file and print out the energies again
    uint i = 0;
    for(flt curP=startP; curP>P0; curP/=10){
        cout << "P: " << curP << "\n";
        collec.set_pressure(curP);
        while (true) {
            for(uint j=0; j<1000; j++){
                collec.timestep();
            }
            i++;
            writefile(atoms, *obox);

            flt pdiff = collec.pressure() / curP - 1.0;
            flt force_err = 0;

            for(uint k=0; k<atoms.size(); k++){
                flt fmag = atoms[k].f.mag();
                if(fmag > force_err) force_err = fmag;
            }

            cout.precision(sizeof(flt));
            cout << i << " H: " << collec.hamiltonian() << " K: " << collec.kinetic() 
                << " U: " << hertzian->energy(*obox) << " phi: " << (Vs/obox->V()) << "\n";
            cout.precision(6);
            cout << "        Pdiff: " << pdiff << " force_err: " << force_err << "\n";

            if(abs(pdiff) < P_frac and force_err < force_max){
                cout << "Done!\n";
                break;
            }
        }
    }
};

void writefile(AtomVec& atoms, OriginBox& obox){
    // The .xyz format is simple:
    // Line 1: [Number of atoms]
    // Line 2: [Comment line, whatever you want, left blank here]
    // Line 3: [element type, here C for carbon]\t[x]\t[y]\t[z]
    // Line 4: [element type, here C for carbon]\t[x]\t[y]\t[z]
    // ...
    // Line N+1: [element type, here C for carbon]\t[x]\t[y]\t[z]
    // Line N+2: [element type, here C for carbon]\t[x]\t[y]\t[z]
    // Line N+3: [Number of atoms]
    // Line N+4: [Comment line, whatever you want, left blank here]
    // Line N+5: [element type, here C for carbon]\t[x]\t[y]\t[z]
    // ...
    // And so on. each set of N atoms/coordinates corresponds to a "frame",
    // which then make a movie. There must be the same N in each frame for VMD.
    // Note that if this is compiled as a 2D simulation, it will leave out
    // the z coordinate, and VMD can't handle that.
    
    ofstream outf;
    outf.open("packing.xyz", ios::out);
    outf.precision(24);
    
    outf << atoms.size() << endl;
    outf << "L=" << obox.L() << endl; // blank line for comment
    for(uint i=0; i<atoms.size(); i++){
        if(i < Ns){
            outf << "C";
        } else {
            outf << "O";
        }
        Vec normloc = obox.diff(Vec::Zero(), atoms[i].x);
        for(uint j=0; j<NDIM; j++){
            outf << "\t" << normloc[j];
        }
        outf << endl;
    };

    // Unnecessary extra:
    // Write a "tcl" file with the box boundaries
    // the "tcl" file is made specifically for VMD
    ofstream pbcfile;
    pbcfile.open("packing.tcl", ios::out);
    pbcfile << "set cell [pbc set {";
    for(uint j=0; j<NDIM; j++){
        pbcfile << obox.box_shape()[j] << " ";
    }
    pbcfile << "} -all];\n";
    pbcfile << "pbc box -toggle -center origin -color red;\n";

    pbcfile << "set natoms [atomselect 0 \"name C\";];\n"
            << "$natoms set radius " << (sigma/2.0) << ";\n"
            << "set natoms [atomselect 0 \"name O\";];\n"
            << "$natoms set radius " << (sigmal/2.0) << ";\n";


    // Now you should be able to run "vmd -e LJatoms-pbc.tcl LJatoms.xyz"
    // and it will show you the movie and also the bounding box
    // if you have .vmdrc in that same directory, you should also be able
    // to toggle the box with the "b" button
};
