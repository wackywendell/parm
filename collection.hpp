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
        Box *box;
        vector<atomgroup*> groups;
        vector<interaction*> interactions;
        vector<statetracker*> trackers;
        vector<constraint*> constraints;
        metagroup atoms;
        void update_trackers();
        void update_constraints();
        virtual flt setForcesGetPressure(bool seta=true);
        
    public:
        collection(Box *box, 
            vector<atomgroup*> groups=vector<atomgroup*>(),
            vector<interaction*> interactions=vector<interaction*>(),
            vector<statetracker*> trackers=vector<statetracker*>(),
            vector<constraint*> constraints=vector<constraint*>());
        
        //Timestepping
        virtual void setForces(bool seta=true);
        virtual void timestep()=0;
        flt dof();
        
        //Stats
        flt potentialenergy();
        flt energy();
        virtual flt temp(bool minuscomv=true);
        virtual flt kinetic();
        virtual flt virial();
        virtual flt pressure();
        Box *getbox(){return box;};
        inline Vec com(){return atoms.com();};
        inline Vec comv(){return atoms.comv();};
        #ifdef VEC3D
        inline Vec angmomentum(const Vec &loc){return atoms.angmomentum(loc, box);};
        inline Vec angmomentum(){return atoms.angmomentum(com(), box);};
        #elif defined VEC2D
        inline flt angmomentum(const Vec &loc){return atoms.angmomentum(loc, box);};
        inline flt angmomentum(){return atoms.angmomentum(com(), box);};
        #endif
        flt gyradius(); // Radius of gyration
        virtual ~collection(){};
        
        void resetcomv(){atoms.resetcomv();};
        void resetL(){atoms.resetL(box);};
        void scaleVs(flt scaleby);
        void scaleVelocitiesT(flt T);
        void scaleVelocitiesE(flt E);
        
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
        StaticCollec(Box *box, vector<atomgroup*> groups,
            vector<interaction*> interactions=vector<interaction*>(),
            vector<statetracker*> trackers=vector<statetracker*>(),
            vector<constraint*> constraints=vector<constraint*>())
                            : collection(box, groups, interactions, trackers, constraints){};
        virtual void timestep(){};
        void update(){update_trackers(); update_constraints();};
};

class collectionSol : public collection {
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
        collectionSol(Box *box, const flt dt, const flt damping, const flt desiredT, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>());
        void changeT(const flt newdt, const flt damp, const flt desiredT){
            dt = newdt; damping = damp; desT = desiredT; setCs();};
        void timestep();
        //void seed(uint n){gauss.seed(n);};
        //void seed(){gauss.seed();};
};

class collectionSolHT : public collection {
    /** for use in solution, with damped forces and random forces
      * 
      * Treats atom.f as the "configurational force", and atom.a as
      * the acceleration due to atom.f + random force.
      * Note that d²r/dt² = atom.f + random force (- damping*v), but we
      * don't include that.
      * 
      * 1) find positions: (x0,v0,a0,f0) -> (x, v, a, f)
      *         sets atom.x from previous force, velocity, and position
      *         x = x0 + dt v0 + 1/2 dt^2 a0 - damping*dt²/2m * v0
      * 2) intermediate v: (v0,f0) -> (v1, f0)
      *         v1 = damped(v0) +  dt/2 f0/m
      *                 where damped(vo) = v0*(1-h*damping/2m + (h*damping/m)²/4)
      * 3) setForces(): (x, f0) -> (x,f)
      *         reset and set the forces
      *         f = grad(V(x))
      * 4) Acceleration: (f, a0) -> (f, a1)
      *         a1 = f/m + gaussian
      *         set acceleration given the forces, adding in random pieces, no damping
      * 5) Finish v and a: (v1, f, a1) -> (v, f, a)
      *         v = v1 + dt/2 a1
      *     
      
     **/
    protected:
        flt dt;
        flt damping;
        flt desT; // desired temperature
        gaussVec gauss;
        void setGauss();
    
    public:
        collectionSolHT(Box *box, const flt dt, const flt damping, const flt desiredT, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>());
        void changeT(const flt newdt, const flt damp, const flt desiredT){
            dt = newdt; damping = damp; desT = desiredT; setGauss();};
        void timestep();
        //void seed(uint n){gauss.seed(n);};
        //void seed(){gauss.seed();};
};

class collectionVerlet : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        
    public:
        collectionVerlet(Box *box, const flt dt, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), dt(dt){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionOverdamped : public collection {
    // over-damped simulation, v = gamma * f
    protected:
        flt dt, gamma;
        
    public:
        collectionOverdamped(Box *box, const flt dt, const float gamma,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints),
                dt(dt), gamma(gamma){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionConjGradient : public collection {
    // over-damped simulation, v = gamma * f
    protected:
        flt dt;
        
    public:
        collectionConjGradient(Box *box, const flt dt,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints),
                dt(dt){};
        void timestep();
        void timestepNewton();
        void reset();
        void setdt(flt newdt){dt=newdt;};
};

class collectionConjGradientBox : public collection {
    // Conjugate-Gradient energy minimization, with 
    // H = H0(x₁, x₂, …, L) + P V
    // More specifically, we take the Nose-Hoover NPH hamiltonian,
    // H = ½m V^⅔ Σṡᵢ² + ½Q V̇² + U(V^⅔ ⃗sᵢ…) + P₀ V
    // and E = U(V^⅔ ⃗sᵢ…) + P₀ V
    // We minimize using ⃗sᵢ and V as the two variables
    /// NOTE: this CANNOT be used with neighbor lists, as it modifies
    /// the box size as it goes; either that, or you have to update
    /// the neighbor list more carefully each time.
    protected:
        flt dt;
        flt P0, kappaV;
        flt hV, FV, lastFV, dV;
        flt maxdV;
        
    public:
        collectionConjGradientBox(OriginBox *box, const flt dt,
                const flt P0, const flt kappaV=1,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints),
                dt(dt), P0(P0), kappaV(kappaV), hV(0), FV(0), lastFV(0),
                maxdV(-1){};
        
        flt kinetic();
        
        void timestep();
        void timestepBox();
        void timestepAtoms();
        void reset();
        void resize(flt V);
        void setdt(flt newdt){dt=newdt; reset();};
        void setP(flt P){P0 = P; reset();};
        void setMaxdV(flt diff){maxdV = diff;};
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
    // NVT
    protected:
        flt dt, Q, T;
        flt xi, lns;
        
    public:
        collectionNoseHoover(Box *box, const flt dt, const flt Q, const flt T, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), 
            dt(dt), Q(Q), T(T){
                xi = 0; lns = 0; 
            };
        void setdt(flt newdt){dt=newdt;};
        void setQ(flt newQ){Q=newQ;};
        void resetBath(){xi=0;lns=0;};
        
        void timestep();
        flt Hamiltonian();
        flt getxi(){return xi;};
        flt getlns(){return lns;};
};

class collectionGaussianT : public collection {
    // Gaussian Constraint thermostat
    // NVT
    protected:
        flt dt, Q;
        flt xi;
        flt setxi();
        
    public:
        collectionGaussianT(Box *box, const flt dt, const flt Q, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), 
            dt(dt), Q(Q){};
        void setdt(flt newdt){dt=newdt;};
        void setQ(flt newQ){Q=newQ;};
        void setForces(bool seta=true){setForces(true,true);};
        void setForces(bool seta, bool setxi);
        void timestep();
};

class collectionGear3A : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        
    public:
        collectionGear3A(Box *box, const flt dt, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), dt(dt){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionGear4A : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs;
        void resetbs(){
            uint Natoms = 0;
            vector<atomgroup*>::iterator git;
            for(git = groups.begin(); git<groups.end(); git++){
                Natoms += (*git)->size();
            };
            bs.resize(Natoms, Vec());
        }
        
    public:
        collectionGear4A(Box *box, const flt dt, uint ncorrectionsteps,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, 
                        constraints), dt(dt), ncorrec(ncorrectionsteps){
                resetbs();
            };
        collectionGear4A(Box *box, const flt dt,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
                collection(box, groups, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbs();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionGear5A : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs, cs;
        void resetbcs(){
            uint Natoms = 0;
            vector<atomgroup*>::iterator git;
            for(git = groups.begin(); git<groups.end(); git++){
                Natoms += (*git)->size();
            };
            bs.resize(Natoms, Vec());
            cs.resize(Natoms, Vec());
        }
        
    public:
        collectionGear5A(Box *box, const flt dt, uint ncorrectionsteps,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), 
                        dt(dt), ncorrec(ncorrectionsteps){resetbcs();};
        collectionGear5A(Box *box, const flt dt,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
                collection(box, groups, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbcs();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionGear6A : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs, cs, ds;
        void resetbcds(){
            uint Natoms = 0;
            vector<atomgroup*>::iterator git;
            for(git = groups.begin(); git<groups.end(); git++){
                Natoms += (*git)->size();
            };
            bs.clear(); cs.clear(); ds.clear();
            bs.resize(Natoms, Vec());
            cs.resize(Natoms, Vec());
            ds.resize(Natoms, Vec());
        }
        
    public:
        collectionGear6A(Box *box, const flt dt, uint ncorrectionsteps,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), 
                        dt(dt), ncorrec(ncorrectionsteps){resetbcds();};
        collectionGear6A(Box *box, const flt dt,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
                collection(box, groups, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbcds();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

struct atomRK4 : atom {
    Vec Kxa, Kxb, Kxc, Kxd, Kva, Kvb, Kvc, Kvd;
};

class atomvecRK4 : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms.
    private:
        atomRK4* atoms;
        uint sz;
    public:
        atomvecRK4(vector<flt> masses) : sz(masses.size()){
            atoms = new atomRK4[sz];
            for(uint i=0; i < sz; i++) atoms[i].m = masses[i];
        };
        atomvecRK4(atomgroup &g) : sz(g.size()){
            atoms = new atomRK4[sz];
            for(uint i=0; i < sz; i++){
                (atom &) atoms[i] = g[i];
                //~ if(i > 1) cout << i-1 << ' ' << atoms[i-1].x << ' ' << atoms[i-1].Kvd << '\n';
                //~ cout << i << ' ' << atoms[i].x << ' ' << atoms[i].Kvd << '\n';
            }
        };
        atom& operator[](cuint n){return atoms[n];};
        atom& operator[](cuint n) const {return atoms[n];};
        //~ inline atomRK4& operator[](cuint n){return atoms[n];};
        //~ inline atomRK4& operator[](cuint n) const {return atoms[n];};
        atom* get(cuint n){if(n>=sz) return NULL; return &(atoms[n]);};
        atomRK4* getRK4(cuint n){if(n>=sz) return NULL; return &(atoms[n]);};
        atomid get_id(atom *a);
        inline atomid get_id(uint n) {
            if (n > sz) return atomid(); return atomid(atoms + n,n);};
        //~ inline flt getmass(cuint n) const{return atoms[n].m;};
        //~ inline void setmass(cuint n, flt m){atoms[n].m = m;};
        inline uint size() const {return sz;};
        ~atomvecRK4(){ delete [] atoms;};
};

class collectionRK4 : public collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        vector<atomgroup*> convertRK4vec(vector<atomvecRK4*> rgroups){
            vector<atomgroup*> v;
            vector<atomvecRK4*>::iterator git;
            for(git = rgroups.begin(); git<rgroups.end(); git++){
                v.push_back((atomgroup*) *git);
            }
            return v;
        };
        
    public:
        collectionRK4(Box *box, const flt dt, 
                vector<atomvecRK4*> rgroups=vector<atomvecRK4*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, convertRK4vec(rgroups), interactions, 
                        trackers, constraints), dt(dt){
                setForces();
                update_constraints();
                vector<atomgroup*>::iterator git;
                for(git = groups.begin(); git<groups.end(); git++){
                    atomgroup &m = **git;
                    for(uint i=0; i<m.size(); i++){
                        atomRK4 & a = (atomRK4 &) m[i];
                        a.a = a.f / m.getmass(i);
                    };
                };
            };
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class collectionGear4NPH : public collection {
    // for use in fixed-E, fixed-NPH simulations
    // Nose-Hoover, right?
    protected:
        flt dt;
        flt P, Q; // goal pressure, damping
        flt dV, ddV, dddV; // that's dV²/dt², dV/dt
        uint ncorrec;
        vector<Vec> bs;
        void resetbs(){
            uint Natoms = 0;
            vector<atomgroup*>::iterator git;
            for(git = groups.begin(); git<groups.end(); git++){
                Natoms += (*git)->size();
            };
            bs.resize(Natoms, Vec());
        }
        
    public:
        collectionGear4NPH(OriginBox *box, const flt dt, const flt P,
                const flt Q, uint ncorrectionsteps,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, 
                        constraints), dt(dt), P(P), Q(Q), dV(0), ddV(0), dddV(0), 
                        ncorrec(ncorrectionsteps){
                resetbs();
            };
        collectionGear4NPH(OriginBox *box, const flt dt, const flt P, const flt Q, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
                collection(box, groups, interactions, trackers, constraints),
                        dt(dt), P(P), Q(Q), dV(0), ddV(0), dddV(0), ncorrec(1) {resetbs();};
        void timestep();
        flt kinetic();
        flt temp(bool minuscomv=true);
        flt Hamiltonian(){
            return kinetic() + (Q/2*dV*dV) + potentialenergy() + P*(((OriginBox*)box)->V());
        }
        flt getdV(){return dV;};
        flt getddV(){return ddV;};
        void setdt(flt newdt){dt=newdt;};
};

class xrpsummer : public fpairxFunct {
    private:
        Box *box;
    public:
        flt xsum, rpxsum, vfsum, rfsum;
        xrpsummer(Box *box) : box(box), xsum(0), rpxsum(0), vfsum(0), rfsum(0){};
        virtual void run (forcepairx*);
        inline void reset(){xsum = 0; rpxsum=0; vfsum=0; rfsum=0;};
};

class collectionGear4NPT : public collection {
    // for use in fixed-NPT simulations
    // Gaussian constraint formulation
    public:
        flt dt;
        xrpsummer xrpsums;
        uint ncorrec;
        flt V1, V2, V3, chi, chixi;
        vector<Vec> xs1, xs2, xs3;
        vector<Vec> vs2, vs3;
        void resetbs(){
            uint Natoms = 0;
            vector<atomgroup*>::iterator git;
            for(git = groups.begin(); git<groups.end(); git++){
                Natoms += (*git)->size();
            };
            xs1.resize(Natoms, Vec());
            xs2.resize(Natoms, Vec());
            xs3.resize(Natoms, Vec());
            vs2.resize(Natoms, Vec());
            vs3.resize(Natoms, Vec());
            V1 = V2 = V3 = 0;
        }
        static vector<interaction*> tointerpair(vector<interactionpairsx*>&);
        
    public:
        collectionGear4NPT(OriginBox *box, const flt dt, uint ncorrectionsteps,
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interactionpairsx*> interactions=vector<interactionpairsx*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, tointerpair(interactions), trackers, 
                        constraints), dt(dt), xrpsums(box),
                        ncorrec(ncorrectionsteps){
                resetbs();
            };
        collectionGear4NPT(OriginBox *box, const flt dt, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interactionpairsx*> interactions=vector<interactionpairsx*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, tointerpair(interactions),
                            trackers, constraints),
                    dt(dt), xrpsums(box), ncorrec(1) {
                resetbs();
            };
        void setForces(bool seta=true);
        void timestep();
};



class collectionVerletNPT : public collection {
    // From Toxvaerd 1993, PRE Vol. 47, No. 1
    protected:
        flt dt;
        flt eta, xidot, lastxidot, lastV;
        vector<Vec> vhalf;
        flt P, QP, T, QT, curP;
        void resetvhalf();
        
    public:
        collectionVerletNPT(OriginBox *box, const flt dt, const flt P,
                const flt QP, const flt T, const flt QT, 
                vector<atomgroup*> groups=vector<atomgroup*>(),
                vector<interaction*> interactions=vector<interaction*>(),
                vector<statetracker*> trackers=vector<statetracker*>(),
                vector<constraint*> constraints=vector<constraint*>()) :
            collection(box, groups, interactions, trackers, constraints), 
            dt(dt), eta(0), xidot(0), lastxidot(0), lastV(box->V()), P(P), 
            QP(QP), T(T), QT(QT), curP(0){resetvhalf();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
        
        
        void resetcomv(){collection::resetcomv(); resetvhalf();};
        void resetL(){collection::resetL(); resetvhalf();};
        void scaleVs(flt scaleby){collection::scaleVs(scaleby); resetvhalf();};
        void scaleVelocitiesT(flt T){collection::scaleVelocitiesT(T); resetvhalf();};
        void scaleVelocitiesE(flt E){collection::scaleVelocitiesE(E); resetvhalf();};
        
        flt geteta(){return eta;};
        flt getxidot(){return xidot;};
        flt getP(){return curP;};
        Vec getvhalf(uint n){return vhalf[n];};
};

#endif
