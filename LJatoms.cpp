#include <iostream>

#include "vecrand.hpp"
#include "interaction.hpp"
#include "collection.hpp"

const flt sigma = 1.0;
const flt epsilon = 1.0;
const flt sigcut = 2.5;
const flt L = 8.0;
const flt dt = 0.001;

int main(){
    seed(); // new random seed each time, for velocities and placement
    
    OriginBox obox(L);
    atomvec atoms(10, 1);
    
    // LJ interaction
    // We'll cut it off at 2.5σ, and have a neighborlist out to 1.4 times that
    neighborlist nl = neighborlist(&obox, sigcut*sigma, 1.4*(sigcut*sigma));
    NListed<LJatomcut, LJAttractPair> LJ = NListed<LJatomcut, LJAttractPair>(&nl);
    
    for (uint i=0; i < atoms.size(); i++){
        flt E0 = LJ.energy(&obox);
        atoms[i].x = obox.randLoc(); // random location in the box
        atoms[i].v = randVec(); // from a Gaussian
        atoms[i].f = Vec();
        atoms[i].a = Vec();
        
        // Add it to the Lennard-Jones potential
        LJ.add(LJatomcut(epsilon, sigma, &atoms[i], sigcut));
        // force an update the neighborlist, so we can get an accurate energy
        nl.update_list(true);
        // If we're overlapping too much, try a new location
        while(LJ.energy(&obox) > E0 + epsilon/2.0){
            atoms[i].x = obox.randLoc(); // random location in the box
            nl.update_list(true);
        }
    }
    
    vector<atomgroup*> allatoms;
    allatoms.push_back(&atoms);
    collectionVerlet collec = collectionVerlet(&obox, dt, allatoms);
    // This is very important! Otherwise the neighborlist won't update!
    collec.addTracker(&nl);
    collec.addInteraction(&LJ);
    
    // remove center of mass velocity, and set a temperature of 2
    collec.resetcomv();
    collec.scaleVelocitiesT(2.0);
    
    cout << "E: " << collec.energy() << " K: " << collec.kinetic() 
        << " U: " << LJ.energy(&obox) << "\n";
    for(uint i=0; i<20; i++){
        for(uint j=0; j<1000; j++){
            collec.timestep();
        }
        cout << "E: " << collec.energy() << " K: " << collec.kinetic()
            << " U: " << LJ.energy(&obox) << "\n";
    }
}
