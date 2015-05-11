#ifndef COLLECTION_H
#define COLLECTION_H

#include "vecrand.hpp"
#include "interaction.hpp"
#include "constraints.hpp"
#include <cstdio>
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>

class Collection : public boost::enable_shared_from_this<Collection> {
    /* A group of AtomGroups and interactions, meant to encapsulate an
     * entire simulation.
     *
     * Adds general simulation time-stepping as well as statistical tracking.
     */
    protected:
        sptr<Box> box;
        sptr<AtomGroup> atoms;
        vector<sptr<Interaction> > interactions;
        vector<sptr<StateTracker> > trackers;
        vector<sptr<Constraint> > constraints;

        void update_trackers();
        void update_constraints();
        virtual flt setForcesGetPressure(bool seta=true);

    public:
        Collection(sptr<Box> box,
            sptr<AtomGroup> atoms,
            vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
            vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
            vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >(),
            bool should_initialize=true);

        virtual void initialize();

        //Timestepping
        virtual void setForces(bool seta=true);
        virtual void timestep()=0;
        flt dof();

        //Stats
        flt potentialenergy();
        flt energy();
        virtual flt temp(bool minuscomv=true);
        virtual flt kinetic(){return atoms->kinetic();};
        virtual flt virial();
        virtual flt pressure();
        sptr<Box> getbox(){return box;};
        inline Vec com(){return atoms->com();};
        inline Vec comv(){return atoms->comv();};
        #ifdef VEC3D
        inline Vec angmomentum(const Vec &loc){return atoms->angmomentum(loc, *box);};
        inline Vec angmomentum(){return atoms->angmomentum(com(), *box);};
        #elif defined VEC2D
        inline flt angmomentum(const Vec &loc){return atoms->angmomentum(loc, *box);};
        inline flt angmomentum(){return atoms->angmomentum(com(), *box);};
        #endif
        flt gyradius(); // Radius of gyration
        virtual ~Collection(){};

        void resetcomv(){atoms->resetcomv();};
        void resetL(){atoms->resetL(*box);};
        void scaleVs(flt scaleby);
        void scaleVelocitiesT(flt T, bool minuscomv=true);
        void scaleVelocitiesE(flt E);

        void addInteraction(sptr<Interaction> inter){
            interactions.push_back(inter);
            update_trackers();
        };
        void addTracker(sptr<StateTracker> track){
            trackers.push_back(track);
            update_trackers();
        };
        void addConstraint(sptr<Constraint> c){
            constraints.push_back(c);
            update_trackers();
        };
        void add(sptr<Interaction> a){addInteraction(a);};
        void add(sptr<StateTracker> a){addTracker(a);};
        void add(sptr<Constraint> a){addConstraint(a);};

        vector<sptr<Interaction> > getInteractions(){return interactions;};

        uint numInteraction(){ return (uint) interactions.size();};
};

class StaticCollec : public Collection {
    public:
        StaticCollec(sptr<Box> box, sptr<AtomGroup> atoms,
            vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
            vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
            vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >())
                            : Collection(box, atoms, interactions, trackers, constraints){};
        virtual void timestep(){};
        void update(){update_trackers(); update_constraints();};
};

class CollectionSol : public Collection {
    /** From Allen and Tildesley, p. 263:
     * r(t + dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + dr_G
     * v(t + dt) = c0 v(t) + c1 - c2) dt a(t) + c2 dt a(t + dt) + dv_G
     *
     * where c0 = exp(- h dt)
     *       c1 = (h dt)^-1 (1 - c0)
     *       c2 = (h dt)^-1 (1 - c1)
     * There are expansions for all 3
     *
     * dr_G and dv_G are drawn from correlated Gaussian distributions, with
     *
     * sigma_r^2 = (h dt)^-1 (2 - (h dt)^-1 (3 - 4 exp(-h dt) + exp(-2 h dt))
     * sigma_v^2 = 1 - exp(-2 h dt)
     * c_{rv} = (h dt)^-1 (1 - exp(-h dt))^2 / sigma_r / sigma_v
     *
     * Those are the unitless versions, multiply by dt^2 kB T/m, kB T/m, and dt kB T/m respectively
     * to get the unit-full versions in the book
    **/

    protected:
        BivariateGauss gauss;
        flt dt;
        flt damping;
        flt desT; // desired temperature
        flt sigmar, sigmav, corr; // note that this is sigmar/sqrt(T/m), same for sigmav
                                  // corr is unitless, and correct
        flt c0, c1, c2; // from Allen and Tildesley
        void setCs();

    public:
        CollectionSol(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt damping, const flt desiredT,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >());
        void changeT(const flt damp, const flt desiredT){
            damping = damp; desT = desiredT; setCs();};
        void setdt(const flt newdt){dt = newdt; setCs();};
        void timestep();
        //void seed(uint n){gauss.seed(n);};
        //void seed(){gauss.seed();};
};

class CollectionDamped : public Collection {
    /** Based on CollectionSol, above.
     * From Allen and Tildesley, p. 263:
     * r(t + dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + dr_G
     * v(t + dt) = c0 v(t) + (c1 - c2) dt a(t) + c2 dt a(t + dt) + dv_G
     *
     * where c0 = exp(- h dt)
     *       c1 = (h dt)^-1 (1 - c0)
     *       c2 = (h dt)^-1 (1 - c1)
     * There are expansions for all 3
     *
     * dr_G and dv_G are drawn from correlated Gaussian distributions, which we drop for this one.
    **/

    protected:
        flt dt;
        flt damping;
        flt c0, c1, c2;
        void setCs();

    public:
        CollectionDamped(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt damping,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >());
        void changeDamp(const flt damp){
            damping = damp; setCs();};
        void setdt(const flt newdt){dt = newdt; setCs();};
        void timestep();
};

class CollectionSolHT : public Collection {
    /** for use in solution, with damped forces and random forces
      *
      * Treats Atom.f as the "configurational force", and Atom.a as
      * the acceleration due to Atom.f + random force.
      * Note that d²r/dt² = Atom.f + random force (- damping*v), but we
      * don't include that.
      *
      * 1) find positions: (x0,v0,a0,f0) -> (x, v, a, f)
      *         sets Atom.x from previous force, velocity, and position
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
        GaussVec gauss;
        void setGauss();

    public:
        CollectionSolHT(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt damping, const flt desiredT,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >());
        void changeT(const flt newdt, const flt damp, const flt desiredT){
            dt = newdt; damping = damp; desT = desiredT; setGauss();};
        void timestep();
        //void seed(uint n){gauss.seed(n);};
        //void seed(){gauss.seed();};
};

/**
\example LJatoms.cpp
*/
class CollectionVerlet : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;

    public:
        CollectionVerlet(sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints), dt(dt){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionOverdamped : public Collection {
    // over-damped simulation, v = gamma * f
    protected:
        flt dt, gamma;

    public:
        CollectionOverdamped(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt gamma=1.0,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
                dt(dt), gamma(gamma){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionConjGradient : public Collection {
    // over-damped simulation, v = gamma * f
    protected:
        flt dt;

    public:
        CollectionConjGradient(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
                dt(dt){};
        void timestep();
        void timestepNewton();
        void reset();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionConjGradientBox : public Collection {
    // Conjugate-Gradient energy minimization, with
    // H = H0(x₁, x₂, …, L) + P V
    // More specifically, we take the Nose-Hoover NPH hamiltonian,
    // H = ½m V^⅔ Σṡᵢ² + ½Q V̇² + U(V^⅓ ⃗sᵢ…) + P₀ V
    // and E = U(V^⅓ ⃗sᵢ…) + P₀ V
    // We minimize using ⃗sᵢ and ln V as the two variables
    /// NOTE: this CANNOT be used with neighbor lists, as it modifies
    /// the box size as it goes; either that, or you have to update
    /// the neighbor list more carefully each time.
    protected:
        flt dt;
        flt P0, kappaV;
        flt hV, FV, lastFV, dV;
        flt maxdV;

    public:
        CollectionConjGradientBox(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt P0, const flt kappaV=1.0,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
                dt(dt), P0(P0), kappaV(kappaV), hV(0), FV(0), lastFV(0), dV(0),
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

/**
Conjugate-Gradient energy minimization.

\example packer.cpp
*/
class CollectionNLCG : public Collection {
    // Conjugate-Gradient energy minimization, with
    // H = H0(x₁, x₂, …, L) + P V
    // More specifically, we take the Nose-Hoover NPH hamiltonian,
    // H = ½m V^⅔ Σṡᵢ² + ½Q V̇² + U(V^⅔ ⃗sᵢ…) + P₀ V
    // and E = U(V^⅔ ⃗sᵢ…) + P₀ V
    // We minimize using ⃗sᵢ and κ ln V as the dN+1 variables

    public:
        // Parameters
        flt dt;
        flt seceps;
        uint secmax;
        flt kappa;
        flt alphamax, afrac, dxmax, stepmax, kmax;

        // Goal pressure
        flt P0;

        // To keep between iterations
        flt Knew;
        flt k;
        flt vl, fl, al;

        // For tracking purposes
        flt alpha, beta, betaused, dxsum, alphavmax, maxdV;
        uint sec;

        void stepx(flt dx);
        flt getLsq();
        flt fdota();
        flt fdotf();
        flt fdotv();
        flt vdotv();
        //~ void resizedl(flt dl);

    public:
        CollectionNLCG(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt P0,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >(),
                const flt kappa=10.0, const flt kmax=1000,
                const uint secmax=40, const flt seceps = 1e-20);

        flt kinetic();  // Note: masses are ignored
        flt pressure();
        flt Hamiltonian();
        void setForces(bool seta=true){setForces(seta,true);};
        void setForces(bool seta, bool setV);

        void timestep();
        void descend(); // use steepest descent
        void reset();
        void resize(flt V);

        void setdt(flt newdt){dt=newdt; reset();};
        void setP(flt P){P0 = P; reset();};
        void setkappa(flt k){kappa=k; reset();};
        void setamax(flt a){alphamax=a;};
        void setafrac(flt a){afrac=a;};
        void setdxmax(flt d){dxmax=d;};
        void setstepmax(flt m){stepmax=m;};

};



class CollectionNLCGV : public Collection {
    // Conjugate-Gradient energy minimization, with
    // and E = U(⃗rᵢ…)

    public:
        // Parameters
        flt dt; // initial step attempt
        flt seceps;
        uint secmax;

        // secmax is the max number of iterations within a timestep
        // seceps is the minimum alpha * v before we call it "close enough"
        flt alphamax, afrac, dxmax, stepmax, kmax;
        // alpha is how much larger the next step is, proportionally
        // if alpha < afrac * dxsum, we break, its "close enough"
        // alphamax is the largest proportion
        // dxmax is the maximum "step" you take in one "timestep"; the full
        //    step is v * dxsum
        // stepmax is the largest v*dxsum allowed
        // We reset after kmax iterations

        // To keep between iterations
        flt Knew;
        flt k;
        flt vl, fl, al;

        // For tracking purposes
        flt alpha, beta, betaused, dxsum, alphavmax;
        // alpha is how much bigger one iteration is than the previous
        // beta is how much we use of the previous v
        // betaused is after we take into account 0 <= beta <= 1
        // dxsum is the total of alphas over a timestep; atoms move v*dxsum in one timestep
        // alphavmax = (last alpha) * sqrt(v dot v)
        uint sec;

        void stepx(flt dx);

        flt fdota();
        flt fdotf();
        flt fdotv();
        flt vdotv();
        //~ void resizedl(flt dl);

    public:
        CollectionNLCGV(sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >(),
                const flt kmax=1000, const uint secmax=10,
                const flt seceps = 1e-4);

        flt pressure();

        void reset();
        void descend(); // use steepest descent
        void timestep();

        void setdt(flt newdt){dt=newdt; reset();};
        void setamax(flt a){alphamax=a;};
        void setafrac(flt a){afrac=a;};
        void setdxmax(flt d){dxmax=d;};
        void setstepmax(flt m){stepmax=m;};
};

flt solveCubic1(flt b, flt c, flt d){
    // from Wikipedia
    flt determ = (pow(2*pow(b,2) - 9*b*c + 27*d,2) - 4*pow(b*b - 3*c,3));
    if (determ < 0)
        printf("bad determ: %.4f\n", (double) determ);
    flt firstpartundercube = (2*pow(b,3) - 9*b*c + 27*d)/2;
    flt secondpartundercube = sqrt(determ)/2;
    if (firstpartundercube < secondpartundercube)
        printf("bad pairs under cube: %.4f < %.4f (%.4f)\n",
                    (double) firstpartundercube, (double) secondpartundercube,
                    (double) (firstpartundercube - secondpartundercube));
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
        printf("Multiple Answers: %.4f, %.4f\n", (double) Q,(double) R);
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
            printf("Multiple Answers: %.4f, %.4f\n", (double) Q,(double) R);
            //~ printf("pi: %.4f\n", M_2_PI);
            //~ printf("theta %.4f : %.4f, %.4f, %.4f\n", theta,
                //~ theta/3, (theta + (2*M_PI))/3, (theta + (4*M_PI))/3);
            printf("%.4f (%.4f), %.4f (%.4f), %.4f (%.4f) : %.4f\n",
                (double) x1,(double) d1,(double) x2,(double) d2,(double) x3,(double) d3, (double) x);
        }
        return x;
    }
    flt R2Q3 = cbrt(sqrt(R2 - Q3) + fabs(R));
    return -(sgn(R)*(R2Q3 + (Q/R2Q3))) - (a1/3);
}

class CollectionNoseHoover : public Collection {
    // NVT
    protected:
        flt dt, Q, T;
        flt xi, lns;

    public:
        CollectionNoseHoover(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt Q, const flt T,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
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

class CollectionGaussianT : public Collection {
    // Gaussian Constraint thermostat
    // NVT
    protected:
        flt dt, Q;
        flt xi;
        flt setxi();

    public:
        CollectionGaussianT(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, const flt Q,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            dt(dt), Q(Q), xi(0){};
        void setdt(flt newdt){dt=newdt;};
        void setQ(flt newQ){Q=newQ;};
        void setForces(bool seta=true){setForces(true,true);};
        void setForces(bool seta, bool setxi);
        void timestep();
};

class CollectionGear3A : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;

    public:
        CollectionGear3A(sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints), dt(dt){};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionGear4A : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs;
        void resetbs(){
            bs.resize(atoms->size(), Vec());
        }

    public:
        CollectionGear4A(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, uint ncorrectionsteps,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers,
                        constraints), dt(dt), ncorrec(ncorrectionsteps){
                resetbs();
            };
        CollectionGear4A(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
                Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbs();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionGear5A : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs, cs;
        void resetbcs(){
            uint Natoms = atoms->size();
            bs.resize(Natoms, Vec());
            cs.resize(Natoms, Vec());
        }

    public:
        CollectionGear5A(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, uint ncorrectionsteps,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), ncorrec(ncorrectionsteps){resetbcs();};
        CollectionGear5A(sptr<Box> box, const flt dt,
                sptr<AtomGroup> atoms,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
                Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbcs();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionGear6A : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        uint ncorrec;
        vector<Vec> bs, cs, ds;
        void resetbcds(){
            uint Natoms = atoms->size();
            bs.clear(); cs.clear(); ds.clear();
            bs.resize(Natoms, Vec());
            cs.resize(Natoms, Vec());
            ds.resize(Natoms, Vec());
        }

    public:
        CollectionGear6A(sptr<Box> box, sptr<AtomGroup> atoms,
                const flt dt, uint ncorrectionsteps,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), ncorrec(ncorrectionsteps){resetbcds();};
        CollectionGear6A(sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
                Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), ncorrec(1) {resetbcds();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

struct RK4data {
    Vec Kxa, Kxb, Kxc, Kxd, Kva, Kvb, Kvc, Kvd;
};

class CollectionRK4 : public Collection {
    // for use in fixed-E simulations
    protected:
        flt dt;
        vector<RK4data> data;

    public:
        CollectionRK4(sptr<Box> box, sptr<AtomGroup> ratoms, const flt dt,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, ratoms, interactions,
                        trackers, constraints), dt(dt), data(ratoms->vec().size()){
                setForces();
                update_constraints();
                for(uint i=0; i<atoms->size(); ++i){
                    Atom& a = (*atoms)[i];
                    a.a = a.f / a.m;
                };
            };
        void timestep();
        void setdt(flt newdt){dt=newdt;};
};

class CollectionGear4NPH : public Collection {
    // for use in fixed-E, fixed-NPH simulations
    // Nose-Hoover, right?
    protected:
        flt dt;
        flt P, Q; // goal pressure, damping
        flt dV, ddV, dddV; // that's dV²/dt², dV/dt
        uint ncorrec;
        vector<Vec> bs;
        void resetbs(){
            bs.resize(atoms->size(), Vec());
        }

    public:
        CollectionGear4NPH(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt P,
                const flt Q, uint ncorrectionsteps,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers,
                        constraints), dt(dt), P(P), Q(Q), dV(0), ddV(0), dddV(0),
                        ncorrec(ncorrectionsteps){
                resetbs();
            };
        CollectionGear4NPH(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt P, const flt Q,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
                Collection(box, atoms, interactions, trackers, constraints),
                        dt(dt), P(P), Q(Q), dV(0), ddV(0), dddV(0), ncorrec(1) {resetbs();};
        void timestep();
        flt kinetic();
        flt temp(bool minuscomv=true);
        flt Hamiltonian(){
            return kinetic() + (Q/2*dV*dV) + potentialenergy() + P*(boost::static_pointer_cast<OriginBox>(box)->V());
        }
        flt getdV(){return dV;};
        flt getddV(){return ddV;};
        void setdt(flt newdt){dt=newdt;};
};

class XRPSummer : public FPairXFunct {
    private:
        sptr<Box> box;
    public:
        flt xsum, rpxsum, vfsum, rfsum;
        XRPSummer(sptr<Box> box) : box(box), xsum(0), rpxsum(0), vfsum(0), rfsum(0){};
        virtual void run (ForcePairX*);
        inline void reset(){xsum = 0; rpxsum=0; vfsum=0; rfsum=0;};
};

class CollectionGear4NPT : public Collection {
    // for use in fixed-NPT simulations
    // Gaussian Constraint formulation
    public:
        flt dt;
        XRPSummer xrpsums;
        uint ncorrec;
        flt V1, V2, V3, chi, chixi;
        vector<Vec> xs1, xs2, xs3;
        vector<Vec> vs2, vs3;
        void resetbs(){
            uint Natoms = atoms->size();
            xs1.resize(Natoms, Vec());
            xs2.resize(Natoms, Vec());
            xs3.resize(Natoms, Vec());
            vs2.resize(Natoms, Vec());
            vs3.resize(Natoms, Vec());
            V1 = V2 = V3 = 0;
        }
        static vector<sptr<Interaction> > tointerpair(vector<sptr<InteractionPairsX> >&);

    public:
        CollectionGear4NPT(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, uint ncorrectionsteps,
                vector<sptr<InteractionPairsX> > interactions=vector<sptr<InteractionPairsX> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, tointerpair(interactions), trackers,
                        constraints), dt(dt), xrpsums(box),
                        ncorrec(ncorrectionsteps), chi(0), chixi(0){
                resetbs();
            };
        CollectionGear4NPT(sptr<OriginBox> box, const flt dt,
                sptr<AtomGroup> atoms,
                vector<sptr<InteractionPairsX> > interactions=vector<sptr<InteractionPairsX> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, tointerpair(interactions),
                            trackers, constraints),
                    dt(dt), xrpsums(box), ncorrec(1), chi(0), chixi(0) {
                resetbs();
            };
        void setForces(bool seta=true);
        void timestep();
};



class CollectionVerletNPT : public Collection {
    // From Toxvaerd 1993, PRE Vol. 47, No. 1, http://dx.doi.org/10.1103/PhysRevE.47.343
    // Parameter equivalences (in the form code: paper):
    // dof() or ndof: g
    // QT: g k T t_\eta^2
    // QP: N k T t_\xi^2
    // where N is the number of particles
    public:
        flt dt;
        flt eta, xidot, lastxidot, lastV;
        flt etasum; // for calculating the "hamiltonian"
        vector<Vec> vhalf;
        flt P, QP, T, QT, curP;
        void resetvhalf();

    public:
        CollectionVerletNPT(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt P,
                const flt QP, const flt T, const flt QT,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            dt(dt), eta(0), xidot(0), lastxidot(0), lastV(box->V()), etasum(0), P(P),
            QP(QP), T(T), QT(QT), curP(0){resetvhalf();};
        void timestep();
        void setdt(flt newdt){dt=newdt;};


        void resetcomv(){Collection::resetcomv(); resetvhalf();};
        void resetL(){Collection::resetL(); resetvhalf();};
        void scaleVs(flt scaleby){Collection::scaleVs(scaleby); resetvhalf();};
        void scaleVelocitiesT(flt T){Collection::scaleVelocitiesT(T); resetvhalf();};
        void scaleVelocitiesE(flt E){Collection::scaleVelocitiesE(E); resetvhalf();};

        flt geteta(){return eta;};
        flt getxidot(){return xidot;};
        flt getP(){return curP;};
        Vec getvhalf(uint n){return vhalf[n];};

        // note that this is a constant of the motion, but not a real hamiltonian
        // and also only such at constant T
        flt Hamiltonian(){
			// regular energy
			flt H = kinetic() + potentialenergy();

			if(QT > 0){
				flt gkT = dof()*T;
				H +=gkT*etasum;
				H+= eta*eta*QT/2;
			}
			return H;
		}
};

struct Event {
    flt t; // when it will occur
    AtomID a; // Atom 1
    AtomID b; // Atom 2

    bool operator<(const Event& other ) const {
        if (t < other.t) { return true;};
        if (t > other.t) { return false;};
        if (a < other.a) { return true;};
        if (a > other.a) { return false;};
        if (b < other.b) { return true;};
        return false;
    };

};

flt get_max(vector<flt> v){
    flt mx = 0.0;
    for(vector<flt>::iterator it=v.begin(); it != v.end(); ++it){
        if(*it > mx) mx = *it;
    }
    return mx;
};

/** Collision-Driven Brownian-Dynamics

\example hardspheres.cpp
*/
class CollectionCDBDgrid : public Collection {
    public:
        flt T;
        flt dt, curt;
        long long numevents;
        set<Event> events; // note that this a sorted binary tree
        vector<flt> atomsizes; /// diameters
        flt edge_epsilon;

        void reset_events(bool force=true);
        void line_advance(flt deltat);

        Grid grid;
        flt gridt; // when it was updated

        Event next_event(AtomID a);

    public:
        CollectionCDBDgrid(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt T,
                vector<flt> sizes = vector<flt>(),
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            T(T), dt(dt), curt(0), numevents(0), atomsizes(sizes), edge_epsilon(1e-8),
            grid(box, atoms, get_max(sizes) * (1 + edge_epsilon*10), 2.0),
            gridt(0) {
            assert(atomsizes.size() == atoms->size());
        };
        CollectionCDBDgrid(sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt T,
                flt sizes,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            T(T), dt(dt), curt(0), numevents(0), atomsizes(atoms->size(), sizes),
            edge_epsilon(1e-8),
            grid(box, atoms, sizes * (1 + edge_epsilon*10), 2.0), gridt(0) {
            assert(atomsizes.size() == atoms->size());
        };

        void update_grid(bool force=true);
        Grid &get_grid(){return grid;};
        flt get_epsilon(){return edge_epsilon;};
        void set_epsilon(flt eps){edge_epsilon = eps;};
        void reset_velocities();
        bool take_step(flt tlim=-1); // returns true if it collides, false if it hits the tlim
        void timestep();
        long long events_processed(){return numevents;}; // only counts collision-type events
};

/// Collision-Driven Brownian-Dynamics
class CollectionCDBD : public Collection {
    protected:
        flt T;
        flt dt, curt;
        long long numevents;
        set<Event> events; // note that this a sorted binary tree
        vector<flt> atomsizes; /// diameters

        void reset_events();
        void line_advance(flt deltat);

    public:
        CollectionCDBD(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt T,
                vector<flt> sizes = vector<flt>(),
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            T(T), dt(dt), curt(0), numevents(0), atomsizes(sizes) {
            assert(atomsizes.size() == atoms->size());
        };
        CollectionCDBD(sptr<OriginBox> box, sptr<AtomGroup> atoms,
                const flt dt, const flt T, flt sizes,
                vector<sptr<Interaction> > interactions=vector<sptr<Interaction> >(),
                vector<sptr<StateTracker> > trackers=vector<sptr<StateTracker> >(),
                vector<sptr<Constraint> > constraints=vector<sptr<Constraint> >()) :
            Collection(box, atoms, interactions, trackers, constraints),
            T(T), dt(dt), curt(0), numevents(0), atomsizes(atoms->size(), sizes) {
            assert(atomsizes.size() == atoms->size());
        };

        void reset_velocities();
        bool take_step(flt tlim=-1); // returns true if it collides, false if it hits the tlim
        void timestep();
        long long events_processed(){return numevents;}; // only counts collision-type events
};
#endif
