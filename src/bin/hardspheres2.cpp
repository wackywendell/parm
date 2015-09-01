#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include <string>
#include <vector>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>

#include "vecrand.hpp"
#include "interaction.hpp"
#include "collection.hpp"
using namespace std;

/** Hard-Sphere Simulation

Runs a hard-sphere simulation, keeping track of MSDs.

Several steps:
1. Place particles randomly
2. "Relax" them with a soft-sphere potential, so they are not on top of 
        each other
3. Run hard-sphere simulation for 'tottime' steps, without taking data, 
        to equilibrate
4. Run hard-sphere simulation for 'tottime' steps, taking data

Several outputs:
1. A .xyz file, which has locations for the particles, for use with VMD
2. a .tcl file, with Atom sizes, for use with VMD
3. a .msd file, with MSDs for the particles

**/

// Some constants
// Note that "flt" is a "floating point number", defaults to "double"
// but can be compiled as "long double" if needed
//const flt phi = 0.27;
const flt sigma = 1.0;
//const flt dt0 = 0.0001;
//const flt dt1 = 1.0;
// uint is just short for unsigned int
//const uint Natoms = 24;

// for the "relaxing" stage
const flt epsilon = 1.0;
const flt sigcut = 1.0;

// for the hard-sphere run
const flt T = 1.0;
//const flt dt = 0.02; // time between velocity reshuffling
//const flt tottime = 200000; // total number of steps to run
//const flt sizeratio = 2.0;

// Data parameters
const flt nMSDs = 400; // this is the goal; duplicates removed, so it will be less

// Other parameters
const uint printn = 100;

// Function for writing a frame to an open file (see below)
void writefile(ofstream& outf, AtomVec& atoms, Box& bx);
// Function for writing box size / atomsize to file (see below)
void writetcl(const char *fname, OriginBox& box, vector<flt> &atomsizes);

int main(int argc, char **argv){
    int c;
    
    flt phi=.40;
    flt sizeratio=2.0;
    flt dt=0.02;
    int Natoms=40;
    int tottime=200000;
    string outname = "hardspheres.msd";

    opterr = 0;
    while ((c = getopt (argc, argv, "n:f:s:d:t:o:")) != -1) switch (c){
        case 'n': //Number of particles
            Natoms = atol(optarg);
            break;
        case 'f': //Phi -Packing fraction
            phi = atof(optarg);
            break;
        case 's': //Size ratio of large particle to small particle
            sizeratio = atof(optarg);
            break;
        case 'd': //dt of integration of simulation
            dt = atol(optarg);
            break;
        case 't': // total time of simulation
            tottime = atof(optarg);
            break;
        case 'o':
            outname = optarg;
            break;
        default:
            abort ();
    }
  	
    
    //~ srand ( time(NULL) );
    //~ srand(10000);
    ofstream myfile;
    myfile.open(outname.c_str());
    cout << "Running..." << endl;

    string xyzname=outname.substr(0,outname.size()-4)+".xyz";
    string tclname=outname.substr(0,outname.size()-4)+".tcl";

    // new random seed each time, for velocities and placement
    seed();
    
    vector<flt> atomsizes(Natoms, sigma);
    atomsizes[0] = sigma * sizeratio;
    vector<flt> atommasses(Natoms, 1);
    atommasses[0] = pow(sizeratio, NDIM); // NDIM is "number of dimensions" (i.e., 3)
    
    flt Vs = 0;
    for(uint i=0; i<atomsizes.size(); ++i){
        Vs += pow(atomsizes[i], NDIM) * M_PI / 2 / NDIM;
    };
    
    flt L = pow(Vs/phi, OVERNDIM); // OVERNDIM is "1/number of dimensions" (i.e., 1.0/3.0)
    assert(L > (atomsizes[0] + atomsizes[1]));
    // if the box is not wider than two atoms, this algorithm will have trouble
    
    cout << "phi " << (Vs/pow(L,NDIM)) << "\n";
    
    // Create the bounding box (sides of length L), and a "vector" of Natoms atoms
    boost::shared_ptr<OriginBox> obox(new OriginBox(L));
    boost::shared_ptr<Box> boxptr = boost::static_pointer_cast<Box>(obox);
    boost::shared_ptr<AtomVec> atomptr(new AtomVec(atommasses));
    AtomVec & atoms = *atomptr;
    
    
    ////////////////////////////////////////////////////////////////////
    // The relaxing stage
    
    // Repulsion Interaction, for the relaxing stage
    boost::shared_ptr<NeighborList> nl(new NeighborList(obox, sigcut*sigma, 1.4*(sigcut*sigma)));
    boost::shared_ptr<NListed<EpsSigExpAtom, RepulsionPair> > hertz(
                        new NListed<EpsSigExpAtom, RepulsionPair>(nl));
    // ^ this is the Interaction
    
    // A Repulsion Interaction has energy 
    //    Uij = {1/2.5 * |sigma_ij - r_ij|^2.5      r_ij < sigma_ij
    //           0                                  r_ij >= sigma_ij
    
    // Note that NListed is a class template; its an Interaction that
    // takes various structs as template parameters, and turns them into 
    // a neighbor Interaction
    // See interaction.hpp for a whole bunch of them
    // Also note that the NeighborList is not the same as the
    // "neighborlisted" Interaction; multiple interactions can use the
    // same NeighborList
    
    // Now we run through all the atoms, set their positions / velocities,
    // and add them to the hertzian Interaction (i.e., to the neighbor list)
    for (uint i=0; i < atoms.size(); i++){
        // we track energy to see if things are overlapping
        atoms[i].x = obox->rand_loc(); // random location in the box
        atoms[i].v = Vec::Zero();
        atoms[i].f = Vec::Zero();
        atoms[i].a = Vec::Zero();
        
        flt cursigma = sigma;
        if(i==0) (cursigma = sigma*sizeratio);
        // Add it to the potential
        hertz->add(EpsSigExpAtom(&atoms[i], epsilon, cursigma, sigcut));
    }
    // force an update the NeighborList, just to make sure
    nl->update_list(true);
    
    //Now we make our "Collection"
    CollectionOverdamped collec0 = CollectionOverdamped(boxptr, atomptr, 1e-3);
    // CollectionOverdamped uses "overdamped" dynamics, meaning that 
    // dx/dt = γ f
    // here we use γ = 1, because it doesn't matter
    // Note that this Collection always goes down the potential energy gradient
    
    // This is very important! Otherwise the NeighborList won't update!
    collec0.add_tracker(nl);
    // And add the Interaction
    collec0.add_interaction(hertz);
    
    // subtract off center of mass velocity, and set a total energy
    collec0.reset_com_velocity();
    
    // Potential energy per Atom
    int swtch = 0; //Switch to help terminate simulations of unattainable packing fractions
    flt U = collec0.potentialenergy()/Natoms;
    flt U1 = U; //U at prior timestep
    flt U0 = U; //Initial potential energy
    cout    << "Energy0:  " << U << "\n";
    // note that U is > 0 if any atoms overlap, and U = 0 when no atoms overlap
    // so we run
    
    while(swtch == 0){
	if (U==0){
		swtch=1;
		}	
	else if (U<=U0){
		for(uint i=0; i<1000; ++i) collec0.timestep();
		U1=U;
		U = collec0.potentialenergy()/Natoms;
		cout    << "Energy:  " << U << "\n";
		}
	else if (U>.99*U1){
		swtch=1;
		}
    }
    cout    << "Energyf:  " << U << "\n";
    
    if (U==0){
    ////////////////////////////////////////////////////////////////////
    // The equilibration stage
    
    // This is the Brownian motion, hard-sphere collider
    CollectionCDBD collec = CollectionCDBD(obox, atomptr, dt, T, atomsizes);
    
    cout << "Equilibrating... \n";
    for(uint i=0; i<printn; ++i){
        for(uint j=0; j<(tottime/printn); ++j) collec.timestep();
        cout << (printn-i) << ", " << std::flush;
    }
    cout << "Done.\n";
    
    ////////////////////////////////////////////////////////////////////
    // The data stage
    
    // Now we need to make the RsqTracker (which tracks r squared vs. time)
    // we need to tell it which Δts to track, so we do that here.
    // Note that a set is a Collection that is sorted and removes duplicates,
    // so that's why we start with a set.
    // So 
    set<uint> nset;
    
    // Log-spaced Δts
    flt maxlog = log(tottime);
    for(uint n=0; n<nMSDs; n++){
        flt newlog = (maxlog * n) / (nMSDs-1);
        uint newn = ceil(exp(newlog));
        if(newn <= 0) newn = 1;
        if(newn > (uint)tottime) newn = tottime;
        nset.insert(newn);
    }
    
    vector<uint> MSDns(nset.begin(), nset.end());
    cout << "Using " << MSDns.size() << " MSDns, [" << MSDns.front() << "-" << MSDns.back() << "]\n";
    
    boost::shared_ptr<RsqTracker> rsqtracker(new RsqTracker(atomptr, MSDns));
    boost::shared_ptr<StateTracker> rsqptr = boost::static_pointer_cast<StateTracker>(rsqtracker);
    collec.add(rsqptr);
    
    // VMD is good for 3D visualization purposes, and it can read .xyz files
    // see 'writefile' function
    ofstream xyzfile;
    xyzfile.open(xyzname.c_str(), ios::out);
    writefile(xyzfile, atoms, *obox);
    writetcl(tclname.c_str(), *obox, atomsizes); // see below
                                                       // unnecessary extra
    
    // Now we really start running
    cout << "Running... \n";
    for(uint i=0; i<printn; ++i){
        for(uint j=0; j<(tottime/printn); ++j) collec.timestep();
        
        cout << (printn-i) << ", " << std::flush;
        writefile(xyzfile, atoms, *obox);
    }
    cout << "Done, saving MSDs\n";
    
    // Save the MSDs to a file
    ofstream msdfile;
    msdfile.open(outname.c_str(), ios::out);
    // Retrieve the time-averaged r^2 values for each Atom for each Δt
    vector<vector<flt> > MSDmeans = rsqtracker->means();
    
    // This will be a tab-separated file, with the first column being 
    // Δt in time units (not timesteps),
    // second column <r(t) - r(t-Δt)>² for the large particle 
    // in units of the small particle diameters,
    // third column <r(t) - r(t-Δt)>² for the first small particle,
    // 4th column <r(t) - r(t-Δt)>² for the second small particle, etc.
    
    
    for(uint i=0; i<MSDns.size(); i++){
        msdfile << (MSDns[i] * dt);
        for(vector<flt>::iterator it=MSDmeans[i].begin(); it<MSDmeans[i].end(); it++){
            msdfile << '\t' << *it;
        }
        msdfile << "\n";
    } 

// Save the MFDs to a file so that we can compute alpha
    ofstream mfdfile;
    mfdfile.open("hardspheres.mfd", ios::out);
    // Retrieve the time-averaged r^4 values for each Atom for each Δt
    vector<vector<flt> > MFDmeans = rfrtracker->means();  

    // This will be a tab-separated file, with the first column being 
    // Δt in time units (not timesteps),
    // second column <r(t) - r(t-Δt)>^4 for the large particle 
    // in units of the small particle diameters,
    // third column <r(t) - r(t-Δt)>^4 for the first small particle,
    // 4th column <r(t) - r(t-Δt)>^4 for the second small particle, etc.

    for(uint i=0; i<MSDns.size(); i++){
        mfdfile << (MSDns[i] * dt);
        for(vector<flt>::iterator it=MFDmeans[i].begin(); it<MFDmeans[i].end(); it++){
            mfdfile << '\t' << *it;
        }
        mfdfile << "\n";
    } 
    
    // Now you should be able to run "vmd -e hardspheres-pbc.tcl hardspheres.xyz"
    // and it will show you the movie and also the bounding box
    // if you have .vmdrc in that same directory, you should also be able
    // to toggle the box with the "b" button

    cout << "Done." << endl;
    return 0;
    };
};

void writetcl(const char *fname, OriginBox& box, vector<flt> &atomsizes){
    // Unnecessary extra:
    // Write a "tcl" file with the box boundaries and the Atom sizes
    // the "tcl" file is made specifically for VMD
    ofstream pbcfile;
    pbcfile.open(fname, ios::out);
    // large Atom size
    pbcfile << "set natoms [atomselect 0 \"name O\";];\n";
    pbcfile << "$natoms set radius " << (atomsizes[0]/2.0) << ";\n";
    // small Atom size
    pbcfile << "set natoms [atomselect 0 \"name C\";];\n";
    pbcfile << "$natoms set radius " << (atomsizes[1]/2.0) << ";\n\n";
    
    pbcfile << "set cell [pbc set {";
    for(uint j=0; j<NDIM; j++){
        pbcfile << box.box_shape()[j] << " ";
    }
    pbcfile << "} -all];\n";
    pbcfile << "pbc box -toggle -center origin -color red;\n";
    // Now you should be able to run "vmd -e hardspheres-pbc.tcl hardspheres.xyz"
    // and it will show you the movie and also the bounding box
    // if you have .vmdrc in that same directory, you should also be able
    // to toggle the box with the "b" button
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
        if(i == 0){outf << "O";}
        else{outf << "C";};
        Vec normloc = bx.diff(Vec::Zero(), atoms[i].x);
        for(uint j=0; j<NDIM; j++){
            outf << "\t" << normloc[j];
        }
        outf << endl;
    };
};
