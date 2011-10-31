#ifndef COLLECTION_H
#define COLLECTION_H

#include "interaction.hpp"
#include "vecrand.hpp"

class collection {
    /* A group of atomgroups and interactions, meant to encapsulate an 
     * entire simulation.
     * 
     * Adds general simulation time-stepping as well as statistical tracking.
     */
    protected:
        vector<atomgroup*> groups;
        vector<interaction*> interactions;
        
    public:
        collection(vector<atomgroup*> groups=vector<atomgroup*>(),
                            vector<interaction*> interactions=vector<interaction*>());
        
        //Timestepping
        virtual void setForces();
        virtual void timestep()=0;
        
        //Stats
        flt Energy();
        flt Temp();
        flt kinetic();
        Vec com();
        Vec comv();
        flt gyradius(); // Radius of gyration
        ~collection(){};
};

class collectionSol : public collection {
    /** for use in solution, with damped forces and random forces
      * 
      * Treats atom.f as the "configurational force", and atom.a as
      * the acceleration due to atom.f + damping + random force.
      * 
      * 1) find positions: (x0,v0,a0,f0) -> (x, v0, a0, f0)
      *         sets atom.x from previous force, velocity, and position
      *         x = x0 + dt v0 + 1/2 dt^2 a0
      * 2) intermediate v: (v0,f0) -> (v1, f0)
      *         v1 = damped(v0) +  dt/2 f0/m
      *                 where damped(vo) = v0*(1-h*damping/2m + (h*damping/m)^2/4)
      * 3) setForces(): (x, f0) -> (x,f)
      *         reset and set the forces
      *         f = grad(V(x))
      * 4) setAccel(): (f, a0) -> (f, a1)
      *         a1 = f/m + gaussian
      *         set acceleration given the forces, adding in random pieces, no damping
      * 5) Finish v and a: (v1, f, a1) -> (v, f, a)
      *         v = v1 + dt/2 a1
      *         a = a1 - damping * v
      *     
      
     **/
    protected:
        bivariateGauss gauss;
        flt dt;
        flt damping;
        flt desT; // desired temperature
        flt sigmar, sigmav, corr; // note that this is sigmar/sqrt(T/m), same for sigmav
                                  // corr is unitless, and correct
        flt c0, c1, c2; // from Allen and Tildesley
        void setCs();
    
    public:
        collectionSol(const flt dt, const flt damping, const flt desiredT, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>());
        void timestep();
        void seed(uint n){gauss.seed(n);};
        void seed(){gauss.seed();};
};
#endif
