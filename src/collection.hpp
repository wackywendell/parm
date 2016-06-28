#ifndef COLLECTION_H
#define COLLECTION_H

#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <set>
#include <vector>
#include "constraints.hpp"
#include "interaction.hpp"
#include "vecrand.hpp"

/*!
A "Collection" of atoms, the box, and an integrator. Provides general simulation
time-stepping as
well as statistical tracking.

The most useful method is timestep(), which takes one step forward in time; this
is defined
by each subclass separately, as each subclass uses a different integration
scheme.
*/
class Collection {
   protected:
    sptr<Box> box;
    sptr<AtomGroup> atoms;
    vector<sptr<Interaction> > interactions;
    vector<sptr<StateTracker> > trackers;
    vector<sptr<Constraint> > constraints;

    //! To be called immediately after setting particle positions and
    //! velocities; lets
    //! `StateTracker` instances stay updated automatically
    void update_trackers();

    //! To be called approximately after forces have been set. Constraints will
    //! typically
    //! set forces / velocities in some direction to zero, so
    //! `update_constraints` should be
    //! called after all forces have been set, and any external velocity changes
    //! have been made.
    void update_constraint_positions();
    void update_constraint_velocities();
    void update_constraint_forces();
    virtual flt set_forces_get_pressure(bool constraints_and_a = true);

   public:
    Collection(
        sptr<Box> box, sptr<AtomGroup> atoms,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >(),
        bool should_initialize = true);

    virtual void initialize();

    // Timestepping Methods
    // /////////////////////////////////////////////////////////////////////
    //! Set forces. This should be called at the beginning of the simulation,
    //! and will also be
    //! called by timestep()
    virtual void set_forces(bool constraints_and_a = true);

    //! Take one step forward in time
    virtual void timestep() = 0;

    //! Total degrees of freedom. This takes into account constraints, and
    //! whether the
    //! center-of-mass is free.
    flt degrees_of_freedom();

    // Statistical Methods
    // //////////////////////////////////////////////////////////////////////
    flt potential_energy();
    //! Total energy, including both potential and kinetic_energy
    flt energy();
    virtual flt temp(bool minuscomv = true);
    virtual flt kinetic_energy() { return atoms->kinetic_energy(); };
    virtual flt virial();
    virtual flt pressure();
    sptr<Box> get_box() { return box; };
    inline Vec com() { return atoms->com(); };
    inline Vec com_velocity() { return atoms->com_velocity(); };
#ifdef VEC3D
    //! Shortcut to `AtomGroup` method of the same name
    inline Vec angular_momentum(const Vec &loc) {
        return atoms->angular_momentum(loc);
    };
    //! Shortcut to `AtomGroup` method of the same name
    inline Vec angular_momentum() { return atoms->angular_momentum(com()); };
#elif defined VEC2D
    //! Shortcut to `AtomGroup` method of the same name
    inline flt angular_momentum(const Vec &loc) {
        return atoms->angular_momentum(loc);
    };
    //! Shortcut to `AtomGroup` method of the same name
    inline flt angular_momentum() { return atoms->angular_momentum(com()); };
#endif
    flt gyradius();  // Radius of gyration
    virtual ~Collection(){};

    //! Shortcut to `AtomGroup` method of the same name
    void reset_com_velocity() { atoms->reset_com_velocity(); };
    //! Shortcut to `AtomGroup` method of the same name
    void reset_L() { atoms->reset_L(); };
    //! Scale all velocities by a factor
    void scale_velocities(flt scaleby);
    //! Scale all velocities to get to a specific temperature
    void scale_velocities_to_temp(flt T, bool minuscomv = true);
    //! Scale all velocities to get to a specific total energy
    void scale_velocities_to_energy(flt E);

    virtual void add_interaction(sptr<Interaction> inter) {
        interactions.push_back(inter);
        update_trackers();
    };
    virtual void add_tracker(sptr<StateTracker> track) {
        trackers.push_back(track);
        update_trackers();
    };
    virtual void add_constraint(sptr<Constraint> c) {
        constraints.push_back(c);
        update_trackers();
    };
    void add(sptr<Interaction> a) { add_interaction(a); };
    void add(sptr<StateTracker> a) { add_tracker(a); };
    void add(sptr<Constraint> a) { add_constraint(a); };

    vector<sptr<Interaction> > get_interactions() { return interactions; };
};

//! A "static" Collection, that doesn't move.
class StaticCollec : public Collection {
   public:
    StaticCollec(
        sptr<Box> box, sptr<AtomGroup> atoms,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints){};

    //! Does nothing; a no-op.
    virtual void timestep(){};
    void update() {
        update_trackers();
        update_constraint_positions();
        update_constraint_velocities();
        update_constraint_forces();
    };
};

//! A Collection with a "solvent", using the Langevin equation.
/*! The Langevin Equation ((modified with \f$\vec{f}\f$):

\f$\dot{\vec{p}} = - \xi \vec{p} + \vec{f} + \overset{\circ}{\vec{p}}\f$

Where \f$- \xi \vec{p}\f$ is a drag term, \f$\overset{\circ}{\vec{p}}\f$ is a
"random force"
term, and \f$\vec{f}\f$ is the standard force term. When \f$\vec{f} =
\vec{0}\f$, particles undergo
Brownian motion, with \f$\xi = \frac{k_B T}{m D}\f$.

This particular algorithm is from Allen and Tildesley, p. 263:

\f$
\begin{align*}
\vec{r}\left(t+\delta t\right)) & =\vec{r}\left(t\right)+c_{1}\delta
t\vec{v}\left(t\right)+c_{2}\delta t^{2}\vec{a}\left(t\right)+\delta r^{G}\\
\vec{v}\left(t+\delta t\right) &
=c_{0}v\left(t\right)+\left(c_{1}-c_{2}\right)\delta
t\vec{a}\left(t\right)+c_{2}\delta t\vec{a}\left(t+\delta t\right)+\delta v^{G}
\end{align*}
\f$

where

\f$
\begin{align*}
c_0 & = e^{- \xi \delta t} &
c_1 & = \frac{1 - c_0}{\xi \delta t} &
c_2 & = \frac{1 - c_1}{\xi \delta t}
\end{align*}
\f$

There are expansions for all 3; see Allen and Tildesley page 261.

\f$\delta r^G\f$ and \f$\delta v^G\f$ are drawn from correlated Gaussian
distributions, with

\f$
\begin{align*}
\sigma_{r}^{2} & =\left(\xi dt\right)^{-1}\left[2-\left(\xi
dt\right)^{-1}\left(3-4e^{-\xi dt}+e^{-2\xi dt}\right)\right]\\
\sigma_{v}^{2} & =1-e^{-2\xi dt}\\
c_{rv} & =\frac{\left(1-e^{-\xi dt}\right)^{2}}{\xi dt\sigma_{r}\sigma_{v}}
\end{align*}
\f$

Those are the unitless versions, multiply by \f$\frac{\delta t^2 k_B T}{m}\f$,
\f$\frac{k_B
T}{m}\f$, and \f$\frac{\delta t k_B T}{m}\f$ respectively to get the unit-full
versions in the book

*/
class CollectionSol : public Collection {
   protected:
    //! The random number generator
    BivariateGauss gauss;
    flt dt;
    //! Damping coefficient, \f$\xi\f$
    flt damping;
    flt force_mag;
    //! desired temperature
    flt desT;
    //! note that this is sigmar/sqrt(T/m), same for sigmav
    //! corr is unitless, and correct
    flt sigmar, sigmav, corr;
    flt c0, c1, c2;
    //! Set c0, c1, c2, sigmar, sigmav, corr from desT, dt, and damping
    void set_constants();

   public:
    CollectionSol(
        //! Box for the atoms
        sptr<Box> box,
        //! The atoms
        sptr<AtomGroup> atoms,
        //! The timestep \f$\delta t\f$
        const flt dt,
        //! Damping coefficient, \f$\xi\f$
        const flt damping,
        //! The desired temperature \f$T\f$
        const flt desired_temperature,
        //! The interactions, other than brownian motion. These provide
        //! \f$\vec{f}\f$.
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        //! Trackers, such as a NeighborList
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        //! Constraints
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >());
    //! Change the desired damping coefficient \f$\xi\f$ or temperature \f$T\f$.
    void change_temperature(const flt damp, const flt desired_temperature) {
        damping = damp;
        force_mag = damp;
        desT = desired_temperature;
        set_constants();
    };
    void change_force(const flt damp, const flt fmag,
                      const flt desired_temperature) {
        //! Change the timestep \f$\delta t\f$.
        damping = damp;
        force_mag = fmag;
        desT = desired_temperature;
        set_constants();
    };
    void set_dt(const flt newdt) {
        dt = newdt;
        set_constants();
    };
    void timestep();
    // void seed(uint n){gauss.seed(n);};
    // void seed(){gauss.seed();};
};

//! A damped Collection, equivalent to `CollectionSol` but without the random
//! forces.
/*! Uses the Langevin Equation (modified with \f$\vec{f}\f$):

\f$\dot{\vec{p}} = - \xi \vec{p} + \vec{f} + \overset{\circ}{\vec{p}}\f$

Except that here, we use \f$\overset{\circ}{\vec{p}} = 0\f$.
*/
class CollectionDamped : public Collection {
   protected:
    flt dt;
    flt damping;
    flt c0, c1, c2;
    void set_constants();

   public:
    CollectionDamped(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt, const flt damping,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >());
    void change_damping(const flt damp) {
        damping = damp;
        set_constants();
    };
    void set_dt(const flt newdt) {
        dt = newdt;
        set_constants();
    };
    void timestep();
};

/*! for use in solution, with damped forces and random forces.
  *
  * This uses the same physics as `CollectionSol`, but a different algorithm
 from Honeycutt and Thirumalai.
  *
  * Note that the Honeycutt and Thirumalai algorithm is flawed, and this fixes
 it to some degree.
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
  * 3) set_forces(): (x, f0) -> (x,f)
  *         reset and set the forces
  *         f = grad(V(x))
  * 4) Acceleration: (f, a0) -> (f, a1)
  *         a1 = f/m + gaussian
  *         set acceleration given the forces, adding in random pieces, no
 damping
  * 5) Finish v and a: (v1, f, a1) -> (v, f, a)
  *         v = v1 + dt/2 a1
  *

 **/
class CollectionSolHT : public Collection {
   protected:
    flt dt;
    flt damping;
    flt desT;  // desired temperature
    GaussVec gauss;
    void set_gauss();

   public:
    CollectionSolHT(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt, const flt damping,
        const flt desired_temperature,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >());
    void change_temperature(const flt newdt, const flt damp,
                            const flt desired_temperature) {
        dt = newdt;
        damping = damp;
        desT = desired_temperature;
        set_gauss();
    };
    void timestep();
    // void seed(uint n){gauss.seed(n);};
    // void seed(){gauss.seed();};
};

/**
\example LJatoms.cpp
\example LJ.py
\example polymer.py
*/
class CollectionVerlet : public Collection {
    // for use in fixed-E simulations
   protected:
    flt dt;

   public:
    CollectionVerlet(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints), dt(dt){};
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
};

class CollectionOverdamped : public Collection {
    // over-damped simulation, v = gamma * f
   protected:
    flt dt, gamma;

   public:
    CollectionOverdamped(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        const flt gamma = 1.0,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          gamma(gamma){};
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
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

    virtual void stepx(flt dx);
    flt get_length_squared();
    flt fdota();
    flt fdotf();
    flt fdotv();
    flt vdotv();
    //~ void resizedl(flt dl);

   public:
    CollectionNLCG(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt P0,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >(),
        const flt kappa = 10.0, const flt kmax = 1000, const uint secmax = 40,
        const flt seceps = 1e-20);

    flt kinetic_energy();  // Note: masses are ignored
    flt pressure();
    flt hamiltonian();
    void set_forces(bool constraints_and_a = true) {
        set_forces(constraints_and_a, true);
    };
    void set_forces(bool constraints_and_a, bool setV);

    void timestep();
    void descend();  // use steepest descent
    void reset();

    void set_dt(flt newdt) {
        dt = newdt;
        reset();
    };
    void set_pressure_goal(flt P) {
        P0 = P;
        reset();
    };
    flt get_pressure_goal() { return P0; };
    void set_kappa(flt k) {
        kappa = k;
        reset();
    };
    void set_max_alpha(flt a) { alphamax = a; };
    void set_max_alpha_fraction(flt a) { afrac = a; };
    void set_max_dx(flt d) { dxmax = d; };
    void set_max_step(flt m) { stepmax = m; };
};

class CollectionNLCGFixedL : public CollectionNLCG {
   public:
    CollectionNLCGFixedL(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt P0,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >(),
        const flt kappa = 10.0, const flt kmax = 1000, const uint secmax = 40,
        const flt seceps = 1e-20)
        : CollectionNLCG(box, atoms, dt, P0, interactions, trackers,
                         constraints, kappa, kmax, secmax, seceps){};
    void stepx(flt dx);
};

class CollectionNLCGV : public Collection {
    // Conjugate-Gradient energy minimization, with
    // and E = U(⃗rᵢ…)

   public:
    // Parameters
    flt dt;  // initial step attempt
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
    // dxsum is the total of alphas over a timestep; atoms move v*dxsum in one
    // timestep
    // alphavmax = (last alpha) * sqrt(v dot v)
    uint sec;

    void stepx(flt dx);

    flt fdota();
    flt fdotf();
    flt fdotv();
    flt vdotv();
    //~ void resizedl(flt dl);

   public:
    CollectionNLCGV(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >(),
        const flt kmax = 1000, const uint secmax = 10, const flt seceps = 1e-4);

    flt pressure();

    void reset();
    void descend();  // use steepest descent
    void timestep();

    void set_dt(flt newdt) {
        dt = newdt;
        reset();
    };
    void set_max_alpha(flt a) { alphamax = a; };
    void set_max_alpha_fraction(flt a) { afrac = a; };
    void set_max_dx(flt d) { dxmax = d; };
    void set_max_step(flt m) { stepmax = m; };
};

flt solve_cubic_fast(flt b, flt c, flt d);

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

flt solve_cubic(flt a1, flt a2, flt a3, flt closeto = 0);

class CollectionNoseHoover : public Collection {
    // NVT
   protected:
    flt dt, Q, T;
    flt xi, lns;

   public:
    CollectionNoseHoover(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt, const flt Q,
        const flt T,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          Q(Q),
          T(T) {
        xi = 0;
        lns = 0;
    };
    void set_dt(flt newdt) { dt = newdt; };
    void set_Q(flt newQ) { Q = newQ; };
    void reset_bath() {
        xi = 0;
        lns = 0;
    };

    void timestep();
    flt hamiltonian();
    flt get_xi() { return xi; };
    flt get_lns() { return lns; };
};

class CollectionGaussianT : public Collection {
    // Gaussian Constraint thermostat
    // NVT
   protected:
    flt dt;
    flt xi;
    flt set_xi();

   public:
    CollectionGaussianT(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          xi(0){};
    void set_dt(flt newdt) { dt = newdt; };
    void set_forces(bool constraints_and_a = true) { set_forces(true, true); };
    void set_forces(bool constraints_and_a, bool set_xi);
    void timestep();
};

class CollectionGear3A : public Collection {
    // for use in fixed-E simulations
   protected:
    flt dt;

   public:
    CollectionGear3A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints), dt(dt){};
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
};

class CollectionGear4A : public Collection {
    // for use in fixed-E simulations
   protected:
    flt dt;
    uint ncorrec;
    vector<Vec> bs;
    void resetbs() { bs.resize(atoms->size(), Vec::Zero()); }

   public:
    CollectionGear4A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        uint ncorrectionsteps,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(ncorrectionsteps) {
        resetbs();
    };
    CollectionGear4A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(1) {
        resetbs();
    };
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
};

class CollectionGear5A : public Collection {
    // for use in fixed-E simulations
   protected:
    flt dt;
    uint ncorrec;
    vector<Vec> bs, cs;
    void resetbcs() {
        uint Natoms = atoms->size();
        bs.resize(Natoms, Vec::Zero());
        cs.resize(Natoms, Vec::Zero());
    }

   public:
    CollectionGear5A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        uint ncorrectionsteps,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(ncorrectionsteps) {
        resetbcs();
    };
    CollectionGear5A(
        sptr<Box> box, const flt dt, sptr<AtomGroup> atoms,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(1) {
        resetbcs();
    };
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
};

class CollectionGear6A : public Collection {
    // for use in fixed-E simulations
   protected:
    flt dt;
    uint ncorrec;
    vector<Vec> bs, cs, ds;
    void resetbcds() {
        uint Natoms = atoms->size();
        bs.clear();
        cs.clear();
        ds.clear();
        bs.resize(Natoms, Vec::Zero());
        cs.resize(Natoms, Vec::Zero());
        ds.resize(Natoms, Vec::Zero());
    }

   public:
    CollectionGear6A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        uint ncorrectionsteps,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(ncorrectionsteps) {
        resetbcds();
    };
    CollectionGear6A(
        sptr<Box> box, sptr<AtomGroup> atoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          ncorrec(1) {
        resetbcds();
    };
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
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
    CollectionRK4(
        sptr<Box> box, sptr<AtomGroup> ratoms, const flt dt,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, ratoms, interactions, trackers, constraints),
          dt(dt),
          data(ratoms->vec().size()) {
        set_forces(true);
    };
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };
};

class CollectionGear4NPH : public Collection {
    // for use in fixed-E, fixed-NPH simulations
    // Nose-Hoover, right?
   protected:
    flt dt;
    flt P, Q;           // goal pressure, damping
    flt dV, ddV, dddV;  // that's dV²/dt², dV/dt
    uint ncorrec;
    vector<Vec> bs;
    void resetbs() { bs.resize(atoms->size(), Vec::Zero()); }

   public:
    CollectionGear4NPH(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt P,
        const flt Q, uint ncorrectionsteps,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          P(P),
          Q(Q),
          dV(0),
          ddV(0),
          dddV(0),
          ncorrec(ncorrectionsteps) {
        resetbs();
    };
    CollectionGear4NPH(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt P,
        const flt Q,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          P(P),
          Q(Q),
          dV(0),
          ddV(0),
          dddV(0),
          ncorrec(1) {
        resetbs();
    };
    void timestep();
    flt kinetic_energy();
    flt temp(bool minuscomv = true);
    flt hamiltonian() {
        return kinetic_energy() + (Q / 2 * dV * dV) + potential_energy() +
               P * (boost::static_pointer_cast<OriginBox>(box)->V());
    }
    flt get_dV() { return dV; };
    flt getddV() { return ddV; };
    void set_dt(flt newdt) { dt = newdt; };
};

class XRPSummer : public FPairXFunct {
   private:
    sptr<Box> box;

   public:
    flt xsum, rpxsum, vfsum, rfsum;
    XRPSummer(sptr<Box> box)
        : box(box), xsum(0), rpxsum(0), vfsum(0), rfsum(0){};
    virtual void run(ForcePairX *);
    inline void reset() {
        xsum = 0;
        rpxsum = 0;
        vfsum = 0;
        rfsum = 0;
    };
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
    void resetbs() {
        uint Natoms = atoms->size();
        xs1.resize(Natoms, Vec::Zero());
        xs2.resize(Natoms, Vec::Zero());
        xs3.resize(Natoms, Vec::Zero());
        vs2.resize(Natoms, Vec::Zero());
        vs3.resize(Natoms, Vec::Zero());
        V1 = V2 = V3 = 0;
    }
    static vector<sptr<Interaction> > tointerpair(
        vector<sptr<InteractionPairsX> > &);

   public:
    CollectionGear4NPT(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt,
        uint ncorrectionsteps, vector<sptr<InteractionPairsX> > interactions =
                                   vector<sptr<InteractionPairsX> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, tointerpair(interactions), trackers,
                     constraints),
          dt(dt),
          xrpsums(box),
          ncorrec(ncorrectionsteps),
          chi(0),
          chixi(0) {
        resetbs();
    };
    CollectionGear4NPT(
        sptr<OriginBox> box, const flt dt, sptr<AtomGroup> atoms,
        vector<sptr<InteractionPairsX> > interactions =
            vector<sptr<InteractionPairsX> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, tointerpair(interactions), trackers,
                     constraints),
          dt(dt),
          xrpsums(box),
          ncorrec(1),
          chi(0),
          chixi(0) {
        resetbs();
    };
    void set_forces(bool constraints_and_a = true);
    void timestep();
};

class CollectionVerletNPT : public Collection {
    // From Toxvaerd 1993, PRE Vol. 47, No. 1,
    // http://dx.doi.org/10.1103/PhysRevE.47.343
    // Parameter equivalences (in the form code: paper):
    // degrees_of_freedom() or ndof: g
    // QT: g k T t_\eta^2
    // QP: N k T t_\xi^2
    // where N is the number of particles
   public:
    flt dt;
    flt eta, xidot, lastxidot, lastV;
    flt etasum;  // for calculating the "hamiltonian"
    vector<Vec> vhalf;
    flt P, QP, T, QT, curP;
    void resetvhalf();

   public:
    CollectionVerletNPT(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt P,
        const flt QP, const flt T, const flt QT,
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : Collection(box, atoms, interactions, trackers, constraints),
          dt(dt),
          eta(0),
          xidot(0),
          lastxidot(0),
          lastV(box->V()),
          etasum(0),
          P(P),
          QP(QP),
          T(T),
          QT(QT),
          curP(0) {
        resetvhalf();
    };
    void timestep();
    void set_dt(flt newdt) { dt = newdt; };

    void reset_com_velocity() {
        Collection::reset_com_velocity();
        resetvhalf();
    };
    void reset_L() {
        Collection::reset_L();
        resetvhalf();
    };
    void scale_velocities(flt scaleby) {
        Collection::scale_velocities(scaleby);
        resetvhalf();
    };
    void scale_velocities_to_temp(flt T) {
        Collection::scale_velocities_to_temp(T);
        resetvhalf();
    };
    void scale_velocities_to_energy(flt E) {
        Collection::scale_velocities_to_energy(E);
        resetvhalf();
    };

    flt get_eta() { return eta; };
    flt det_xi_dot() { return xidot; };
    flt get_pressure() { return curP; };
    Vec get_vhalf(uint n) { return vhalf[n]; };

    // note that this is a constant of the motion, but not a real hamiltonian
    // and also only such at constant T
    flt hamiltonian() {
        // regular energy
        flt H = kinetic_energy() + potential_energy();

        if (QT > 0) {
            flt gkT = degrees_of_freedom() * T;
            H += gkT * etasum;
            H += eta * eta * QT / 2;
        }
        return H;
    }
};

struct Event {
    flt t;     // when it will occur
    AtomID a;  // Atom 1
    AtomID b;  // Atom 2

    bool operator<(const Event &other) const {
        if (t < other.t) {
            return true;
        };
        if (t > other.t) {
            return false;
        };
        if (a < other.a) {
            return true;
        };
        if (a > other.a) {
            return false;
        };
        if (b < other.b) {
            return true;
        };
        return false;
    };
};

flt get_max(vector<flt> v);

/// Collision-Driven Dynamics
class CollectionCD : public Collection {
   public:
    flt dt, curt;
    long long numevents;
    set<Event> events;      // note that this a sorted binary tree
    vector<flt> atomsizes;  /// diameters
    vector<sptr<StateTracker> > collision_trackers;

    void reset_events();
    void line_advance(flt deltat);

   public:
    CollectionCD(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt,
        vector<flt> sizes = vector<flt>(),
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >());

    void reset_velocities(flt T);
    bool take_step(flt tlim = -1);  // returns true if it collides, false if it
                                    // hits the tlim
    void timestep();
    virtual void add_tracker(sptr<StateTracker> track);
    long long events_processed() {
        return numevents;
    };  // only counts collision-type events
};

/** Collision-Driven Brownian-Dynamics

\example hardspheres.cpp
*/
class CollectionCDgrid : public Collection {
   public:
    vector<sptr<StateTracker> > collision_trackers;

    flt dt, curt;
    long long numevents;
    set<Event> events;      // note that this a sorted binary tree
    vector<flt> atomsizes;  /// diameters
    flt edge_epsilon;

    void reset_events(bool force = true);
    void line_advance(flt deltat);

    Grid grid;
    flt gridt;  // when it was updated

    Event next_event(AtomID a);

   public:
    CollectionCDgrid(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt,
        vector<flt> sizes = vector<flt>(),
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >());

    virtual void add_tracker(sptr<StateTracker> track);
    void update_grid(bool force = true);
    Grid &get_grid() { return grid; };
    flt get_epsilon() { return edge_epsilon; };
    void set_epsilon(flt eps) { edge_epsilon = eps; };
    void reset_velocities(flt T);
    bool take_step(flt tlim = -1);  // returns true if it collides, false if it
                                    // hits the tlim
    void timestep();
    long long events_processed() {
        return numevents;
    };  // only counts collision-type events
};

/// Collision-Driven Brownian-Dynamics
class CollectionCDBD : public CollectionCD {
   protected:
    flt T;

   public:
    CollectionCDBD(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt T,
        vector<flt> sizes = vector<flt>(),
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : CollectionCD(box, atoms, dt, sizes, interactions, trackers,
                       constraints),
          T(T){};

    void timestep() {
        reset_velocities(T);
        CollectionCD::timestep();
    };
};

/// Collision-Driven Brownian-Dynamics
class CollectionCDBDgrid : public CollectionCDgrid {
   protected:
    flt T;

   public:
    CollectionCDBDgrid(
        sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt dt, const flt T,
        vector<flt> sizes = vector<flt>(),
        vector<sptr<Interaction> > interactions = vector<sptr<Interaction> >(),
        vector<sptr<StateTracker> > trackers = vector<sptr<StateTracker> >(),
        vector<sptr<Constraint> > constraints = vector<sptr<Constraint> >())
        : CollectionCDgrid(box, atoms, dt, sizes, interactions, trackers,
                           constraints),
          T(T){};

    void timestep() {
        reset_velocities(T);
        CollectionCDgrid::timestep();
    };
};
#endif
