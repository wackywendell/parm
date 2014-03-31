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
const flt L = 4.0;
const flt dt = 0.0001;
// uint is just short for unsigned int
const uint Natoms = 10;

// Function for writing a frame to an open file
void writefile(ofstream& outf, atomvec& atoms, Box& bx);

int main(){
    // new random seed each time, for velocities and placement
    seed();
    
    // Create the bounding box (sides of length L), and a "vector" of Natoms atoms
    boost::shared_ptr<OriginBox> obox(new OriginBox(L));
    boost::shared_ptr<atomvec> atomptr(new atomvec(Natoms, 1));
    atomvec & atoms = *atomptr;
    
    // LJ interaction
    // We'll cut it off at 2.5σ, and have a neighborlist out to 1.4 times that
    boost::shared_ptr<neighborlist> nl(new neighborlist(obox, sigcut*sigma, 1.4*(sigcut*sigma)));
    boost::shared_ptr<NListed<LJatomcut, LJCutPair> > LJ(new NListed<LJatomcut, LJCutPair>(nl));
    // ^ this is the interaction
    
    // Note that NListed is a class template; its an interaction that
    // takes various structs as template parameters, and turns them into 
    // a neighbor interaction
    // See interaction.hpp for a whole bunch of them
    // Also note that the neighborlist is not the same as the
    // "neighborlisted" interaction; multiple interactions can use the
    // same neighborlist
    
    // Now we run through all the atoms, set their positions / velocities,
    // and add them to the LJ interaction (i.e., to the neighbor list)
    for (uint i=0; i < atoms.size(); i++){
        // we track energy to see if things are overlapping
        flt E0 = LJ->energy(*obox);
        atoms[i].x = obox->randLoc(); // random location in the box
        atoms[i].v = randVec(); // from a Gaussian
        atoms[i].f = Vec();
        atoms[i].a = Vec();
        
        // Add it to the Lennard-Jones potential
        LJ->add(LJatomcut(epsilon, sigma, &atoms[i], sigcut));
        // force an update the neighborlist, so we can get an accurate energy
        nl->update_list(true);
        // If we're overlapping too much, try a new location
        while(LJ->energy(*obox) > E0 + epsilon/2.0){
            atoms[i].x = obox->randLoc(); // random location in the box
            nl->update_list(true);
        }
    }
    
    //Now we make our "collection"
    collectionVerlet collec = collectionVerlet(boost::static_pointer_cast<Box>(obox), atomptr, dt);
    
    // This is very important! Otherwise the neighborlist won't update!
    collec.addTracker(nl);
    // And add the interaction
    collec.addInteraction(LJ);
    
    // subtract off center of mass velocity, and set a total energy
    collec.resetcomv();
    collec.scaleVelocitiesE(Natoms / 4.0);
    
    // We'll put the energies to stdout and the coordinates to a new file.
    // VMD is good for 3D visualization purposes, and it can read .xyz files
    // see 'writefile' function
    
    ofstream outfile;
    outfile.open("LJatoms.xyz", ios::out);
    writefile(outfile, atoms, *obox);
    
    //Print out total energy, kinetic energy, and potential energy
    cout << "E: " << collec.energy() << " K: " << collec.kinetic() 
        << " U: " << LJ->energy(*obox) << "\n";
    // Run the simulation! And every _ steps, write a frame to the .xyz
    // file and print out the energies again
    for(uint i=0; i<500; i++){
        for(uint j=0; j<1000; j++){
            collec.timestep();
        }
        writefile(outfile, atoms, *obox);
        cout << "E: " << collec.energy() << " K: " << collec.kinetic() 
            << " U: " << LJ->energy(*obox) << "\n";
    }
    
    // Unnecessary extra:
    // Write a "tcl" file with the box boundaries
    // the "tcl" file is made specifically for VMD
    ofstream pbcfile;
    pbcfile.open("LJatoms-pbc.tcl", ios::out);
    pbcfile << "set cell [pbc set {";
    for(uint j=0; j<NDIM; j++){
        pbcfile << obox->boxshape()[j] << " ";
    }
    pbcfile << "} -all];\n";
    pbcfile << "pbc box -toggle -center origin -color red;\n";
    // Now you should be able to run "vmd -e LJatoms-pbc.tcl LJatoms.xyz"
    // and it will show you the movie and also the bounding box
    // if you have .vmdrc in that same directory, you should also be able
    // to toggle the box with the "b" button
};

void writefile(ofstream& outf, atomvec& atoms, Box& bx){
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
        Vec normloc = bx.diff(Vec(), atoms[i].x);
        for(uint j=0; j<NDIM; j++){
            outf << "\t" << normloc[j];
        }
        outf << endl;
    };
};
