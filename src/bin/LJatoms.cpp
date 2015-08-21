#include <iostream>
#include <fstream>

#include "vecrand.hpp"
#include "interaction.hpp"
#include "collection.hpp"

// Some constants
// Note that "flt" is a "floating point number", defaults to "double"
// but can be compiled as "long double" if needed
const flt sigma = 1.0;
const flt epsilon = 1.0;
const flt sigcut = 2.5;
const flt phi = 0.3;
const flt dt = 0.0001;
// uint is just short for unsigned int
const uint Natoms = 400;

// Function for writing a frame to an open file
void writefile(ofstream& outf, AtomVec& atoms, Box& bx);

int main(){
    // new random seed each time, for velocities and placement
    seed();
    
    const flt Vs = Natoms * pow(sigma, NDIM) * M_PI_2 / NDIM;
    const flt L = pow(Vs / phi, OVERNDIM);
    cout << "Using L = " << L << "\n";
    
    // Create the bounding box (sides of length L), and a "vector" of Natoms atoms
    boost::shared_ptr<OriginBox> obox(new OriginBox(L));
    boost::shared_ptr<AtomVec> atomptr(new AtomVec(Natoms, 1.0));
    AtomVec & atoms = *atomptr;
    
    // LJ Interaction
    // We'll cut it off at 2.5σ, and have a NeighborList out to 1.4 times that
    boost::shared_ptr<NListed<LJatomcut, LJCutPair> > 
        LJ(new NListed<LJatomcut, LJCutPair>(obox, atomptr, 0.4*(sigcut*sigma)));
    boost::shared_ptr<NeighborList> nl = LJ->nlist();
    // ^ this is the Interaction
    
    // Note that NListed is a class template; its an Interaction that
    // takes various structs as template parameters, and turns them into 
    // a neighbor Interaction
    // See Interaction.hpp for a whole bunch of them
    // Also note that the NeighborList is not the same as the
    // "neighborlisted" Interaction; multiple interactions can use the
    // same NeighborList
    
    // Now we run through all the atoms, set their positions / velocities,
    // and add them to the LJ Interaction (i.e., to the neighbor list)
    for (uint i=0; i < atoms.size(); i++){
        if((i+1) % 50 == 0) cout << (Natoms - (i+1)) / 50 << ", "; cout.flush();
        
        // we track energy to see if things are overlapping
        flt E0 = LJ->energy(*obox);
        atoms[i].x = obox->rand_loc(); // random location in the box
        atoms[i].v = randVec(); // from a Gaussian
        atoms[i].f = Vec::Zero();
        atoms[i].a = Vec::Zero();
        
        // Add it to the Lennard-Jones potential
        LJ->add(LJatomcut(epsilon, sigma, atoms.get_id(i), sigcut));
        // force an update the NeighborList, so we can get an accurate energy
        nl->update_list(true);
        // If we're overlapping too much, try a new location
        while(LJ->energy(*obox) > E0 + epsilon/2.0){
            atoms[i].x = obox->rand_loc(); // random location in the box
            nl->update_list(true);
        }
    }

    cout << "Starting. Neighborlist contains " << nl->numpairs() << " / " <<
        (atoms.size()*(atoms.size()-1))/2 << " pairs.\n";
    
    //Now we make our "Collection"
    CollectionVerlet collec = CollectionVerlet(boost::static_pointer_cast<Box>(obox), atomptr, dt);
    
    // This is very important! Otherwise the NeighborList won't update!
    collec.add_tracker(nl);
    // And add the Interaction
    collec.add_interaction(LJ);
    
    // subtract off center of mass velocity, and set a total energy
    collec.reset_com_velocity();
    collec.scale_velocities_to_energy(Natoms / 4.0);
    
    // We'll put the energies to stdout and the coordinates to a new file.
    // VMD is good for 3D visualization purposes, and it can read .xyz files
    // see 'writefile' function
    
    ofstream outfile;
    outfile.open("LJatoms.xyz", ios::out);
    writefile(outfile, atoms, *obox);
    
    //Print out total energy, kinetic energy, and potential energy
    cout << "E: " << collec.energy() << " K: " << collec.kinetic_energy() 
        << " U: " << LJ->energy(*obox) << "\n";
    // Run the simulation! And every _ steps, write a frame to the .xyz
    // file and print out the energies again
    for(uint i=0; i<500; i++){
        for(uint j=0; j<1000; j++){
            collec.timestep();
        }
        writefile(outfile, atoms, *obox);
        cout << (500-i) << " E: " << collec.energy() << " K: " << collec.kinetic_energy() 
            << " U: " << LJ->energy(*obox) << "\n";
    }
    
    // Unnecessary extra:
    // Write a "tcl" file with the box boundaries
    // the "tcl" file is made specifically for VMD
    ofstream pbcfile;
    pbcfile.open("LJatoms-pbc.tcl", ios::out);
    pbcfile << "set cell [pbc set {";
    for(uint j=0; j<NDIM; j++){
        pbcfile << obox->box_shape()[j] << " ";
    }
    pbcfile << "} -all];\n";
    pbcfile << "pbc box -toggle -center origin -color red;\n";
    // Now you should be able to run "vmd -e LJatoms-pbc.tcl LJatoms.xyz"
    // and it will show you the movie and also the bounding box
    // if you have .vmdrc in that same directory, you should also be able
    // to toggle the box with the "b" button
    
    cout << "Done. Neighborlist contains " << nl->numpairs() << " / " <<
        (atoms.size()*(atoms.size()-1))/2 << " pairs.\n";
};

void writefile(ofstream& outf, AtomVec& atoms, Box& bx){
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
    
    outf << atoms.size() << endl;
    outf << endl; // blank line for comment
    for(uint i=0; i<atoms.size(); i++){
        outf << "C";
        Vec normloc = bx.diff(Vec::Zero(), atoms[i].x);
        for(uint j=0; j<NDIM; j++){
            outf << "\t" << normloc[j];
        }
        outf << endl;
    };
};
