#ifndef COLLECTION_H
#define COLLECTION_H

#include "vecrand.hpp"
#include "interaction.hpp"
#include "constraints.hpp"
#include <vector>

//~ #include <complex>
//~ typedef complex<double> complx;

class collection {
    /* A group of atomgroups and interactions, meant to encapsulate an 
     * entire simulation.
     * 
     * Adds general simulation time-stepping as well as statistical tracking.
     */
    protected:
        vector<atomgroup*> groups;
        vector<interaction*> interactions;
        vector<statetracker*> trackers;
        vector<constraint*> constraints;
        metagroup atoms;
        void update_trackers();
        void update_constraints();
        
    public:
        collection(vector<atomgroup*> groups=vector<atomgroup*>(),
            vector<interaction*> interactions=vector<interaction*>(),
            vector<statetracker*> trackers=vector<statetracker*>(),
            vector<constraint*> constraints=vector<constraint*>());
        
        //Timestepping
        virtual void setForces();
        virtual void timestep()=0;
        flt dof();
        
        //Stats
        flt potentialenergy();
        flt energy();
        flt temp();
        flt kinetic();
        inline Vec com(){return atoms.com();};
        inline Vec comv(){return atoms.comv();};
        inline Vec angmomentum(const Vec &loc){return atoms.angmomentum(loc);};
        inline Vec angmomentum(){return atoms.angmomentum(com());};
        flt gyradius(); // Radius of gyration
        virtual ~collection(){};
        
        void resetcomv(){atoms.resetcomv();};
        void resetL(){atoms.resetL();};
        void scaleVelocities(flt T);
        
        void addInteraction(interaction* inter){
            interactions.push_back(inter);
            update_trackers();
        };
        void addTracker(statetracker* track){
            trackers.push_back(track);
            update_trackers();
        };
        
        vector<interaction*> getInteractions(){return interactions;};
        
        uint numInteraction(){ return interactions.size();};
};

class StaticCollec : public collection {
    public:
        StaticCollec(vector<atomgroup*> groups,
            vector<interaction*> interactions=vector<interaction*>(),
            vector<statetracker*> trackers=vector<statetracker*>(),
            vector<constraint*> constraints=vector<constraint*>())
                            : collection(groups, interactions, trackers, constraints){};
        virtual void timestep(){};
        void update(){update_trackers(); update_constraints();};
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
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>());
        void changeT(const flt newdt, const flt damp, const flt desiredT){
            dt = newdt; damping = damp; desT = desiredT; setCs();};
        void timestep();
        void seed(uint n){gauss.seed(n);};
        void seed(){gauss.seed();};
};

class collectionVerlet : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        
    public:
        collectionVerlet(const flt dt, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(groups, interactions, trackers, constraints), dt(dt){};
        void timestep();
};

flt solveCubic1(flt b, flt c, flt d){
    // from Wikipedia
    flt determ = (pow(2*pow(b,2) - 9*b*c + 27*d,2) - 4*pow(b*b - 3*c,3));
    if (determ < 0)
        printf("bad determ: %.4f\n", determ);
    flt firstpartundercube = (2*pow(b,3) - 9*b*c + 27*d)/2;
    flt secondpartundercube = sqrt(determ)/2;
    if (firstpartundercube < secondpartundercube) 
        printf("bad pairs under cube: %.4f < %.4f (%.4f)\n", 
                    firstpartundercube, secondpartundercube,
                    firstpartundercube - secondpartundercube);
    flt cuberoot1=cbrt(firstpartundercube + secondpartundercube);
    flt cuberoot2=cbrt(firstpartundercube - secondpartundercube);
    return (-b/3) - (cuberoot1/3) - (cuberoot2/3);
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

flt solveCubic(flt a1, flt a2, flt a3, flt closeto=0){
    // from numerical recipes
    flt Q = (a1*a1 - 3*a2)/9;
    flt Q3 = Q*Q*Q;
    flt R = ((2*a1*a1*a1) - (9*a1*a2) + 27*a3)/54;
    flt R2 = R*R;
    bool determ = (Q3 >= R2);
    if (determ){
        printf("Multiple Answers: %.4f, %.4f\n", Q,R);
        assert(!determ);
        flt theta = acos(R / sqrt(Q3));
        flt sqQ = -2*sqrt(Q);
        flt x1 = sqQ*cos(theta/3) - (a1/3);
        flt x2 = sqQ*cos((theta + (2*M_PI))/3) - (a1/3);
        flt x3 = sqQ*cos((theta + (4*M_PI))/3) - (a1/3);
        flt d1 = fabs(x1 - closeto);
        flt d2 = fabs(x2 - closeto);
        flt d3 = fabs(x3 - closeto);
        //~ printf("pi: %.4f\n", M_2_PI);
        //~ printf("theta %.4f : %.4f, %.4f, %.4f\n", theta, 
                //~ theta/3, (theta + (2*M_PI))/3, (theta + (4*M_PI))/3);
        //~ printf("%.4f (%.4f), %.4f (%.4f), %.4f (%.4f)\n", x1,d1,x2,d2,x3,d3);
        
        //~ if(d1 < d2 and d1 < d3) return x1;
        //~ if(d2 < d1 and d2 < d3) return x2;
        //~ return x3;
        flt x;
        if(d1 < d2 and d1 < d3) x=x1;
        else if(d2 < d1 and d2 < d3) x=x2;
        else x=x3;
        
        #define tol 1e-3
        if(d1 > tol and d2 > tol and d3 > tol){
            printf("Multiple Answers: %.4f, %.4f\n", Q,R);
            //~ printf("pi: %.4f\n", M_2_PI);
            //~ printf("theta %.4f : %.4f, %.4f, %.4f\n", theta, 
                //~ theta/3, (theta + (2*M_PI))/3, (theta + (4*M_PI))/3);
            printf("%.4f (%.4f), %.4f (%.4f), %.4f (%.4f) : %.4f\n", x1,d1,x2,d2,x3,d3, x);
        }
        return x;
    }
    flt R2Q3 = cbrt(sqrt(R2 - Q3) + fabs(R));
    return -(sgn(R)*(R2Q3 + (Q/R2Q3))) - (a1/3);
}

class collectionNoseHoover : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt, Q, T;
        flt xi, lns;
        
    public:
        collectionNoseHoover(const flt dt, const flt Q, const flt T, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(groups, interactions, trackers, constraints), 
            dt(dt), Q(Q), T(T){
                xi = 0; lns = 0; 
            };
        void timestep();
        flt Hamiltonian();
        flt getxi(){return xi;};
        flt getlns(){return lns;};
};

class collectionGaussianT : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt, Q;
        flt xi;
        flt setxi();
        
    public:
        collectionGaussianT(const flt dt, const flt Q, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(groups, interactions, trackers, constraints), 
            dt(dt), Q(Q){};
        void setForces(){setForces(true);};
        void setForces(bool setxi);
        void timestep();
};

#endif
