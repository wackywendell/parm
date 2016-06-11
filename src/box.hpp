#include "vecrand.hpp"

#ifndef BOX_H
#define BOX_H

#include <boost/shared_ptr.hpp>
#include <vector>

#define sptr boost::shared_ptr

typedef const unsigned int cuint;

inline bool to_buffer(vector<Vec*> arr, double* buffer, size_t sizet) {
    if (sizet < NDIM * arr.size()) {
        return false;
    };
    for (uint i = 0; i < arr.size(); i++)
        for (uint j = 0; j < NDIM; j++) {
            buffer[i * NDIM + j] = (double)((*arr[i])[j]);
        }
    return true;
};

class AtomGroup;

/*!
@defgroup basics Basics

The basic classes that almost all simulations will need.
*/

/*!
@defgroup boxes Boxes

Classes and functions related to the boundary conditions of the simulation.
*/

/*!
@defgroup atoms Atoms

Particle related functions.
*/

/*!
@ingroup basics boxes
@brief The virtual interface for the shape of the space and its boundaries.
*/
class Box {
   public:
    //! Distance between two points, given boundary conditions.
    /*!
    This is the main function that Box exists for.
    */
    virtual Vec diff(Vec r1, Vec r2) = 0;
    //! Volume. Can return NaN.
    virtual flt V() = 0;
    virtual ~Box(){};
};

/***********************************************************************
 * Boxes
 */

//! The modulus function for Vec
/*!
@ingroup boxes
*/
#ifdef VEC3D
inline Vec vec_mod(Vec r1, Vec r2) {
    return Vec(remainder(r1[0], r2[0]), remainder(r1[1], r2[1]),
               remainder(r1[2], r2[2]));
};
#endif
#ifdef VEC2D
inline Vec vec_mod(Vec r1, Vec r2) {
    return Vec(remainder(r1[0], r2[0]), remainder(r1[1], r2[1]));
};
#endif

//! An infinite Box, for use with, e.g., sticky conglomerations or proteins.
/*!
@ingroup basics
@ingroup boxes
*/
class InfiniteBox : public Box {
   public:
    //! Simply `r1-r2`.
    Vec diff(Vec r1, Vec r2) { return r1 - r2; };
    //! Returns NaN.
    flt V() { return NAN; };
};

//! A rectilinear Box, with periodic boundary conditions.
/*!
@ingroup basics
*/
class OriginBox : public Box {
   protected:
    Vec boxsize;

   public:
    OriginBox(Vec size) : boxsize(size){};
    Vec diff(Vec r1, Vec r2) { return vec_mod((r1 - r2), boxsize); };
    virtual Vec diff(Vec r1, Vec r2, boost::array<int, NDIM> boxes) {
        Vec dr = r1 - r2;
        for (uint i = 0; i < NDIM; i++) dr[i] -= boxsize[i] * boxes[i];
        return dr;
    };
    //! The `div` function to go with `diff`.
    virtual boost::array<int, NDIM> box_round(Vec r1, Vec r2) {
        boost::array<int, NDIM> boxes;
        Vec dr = r1 - r2;
        for (uint i = 0; i < NDIM; i++)
            boxes[i] = (int)round(dr[i] / boxsize[i]);
        return boxes;
    };
#ifdef VEC3D
    OriginBox(flt L) : boxsize(L, L, L){};
    flt V() { return boxsize[0] * boxsize[1] * boxsize[2]; };
    flt L() { return (boxsize[0] + boxsize[1] + boxsize[2]) / 3.0; };
#endif
#ifdef VEC2D
    OriginBox(flt L) : boxsize(L, L){};
    flt V() { return boxsize[0] * boxsize[1]; };
    flt L() { return (boxsize[0] + boxsize[1]) / 2.0; };
#endif
    //! Resize by a factor. Does not move atoms.
    flt resize(flt factor);
    //! Resize by a factor, and move atoms in the affine direction.
    flt resize(flt factor, AtomGroup& atoms);
    //! Resize to a specific shape. Does not move atoms.
    flt resize_to(Vec newsize);
    //! Resize to a specific shape, and move atoms in the affine direction.
    flt resize_to(Vec newsize, AtomGroup& atoms);
    //! Resize to a specific volume. Does not move atoms.
    flt resize_to_V(flt newV);
    //! Resize to a specific volume, and move atoms in the affine direction.
    flt resize_to_V(flt newV, AtomGroup& atoms);

    //! Resize to a specific length. Does not move atoms.
    flt resize_to_L(flt newL);
    //! Resize to a specific length, and move atoms in the affine direction.
    flt resize_to_L(flt newL, AtomGroup& atoms);
    //! Get a random point in the box.
    Vec rand_loc() {
        Vec v = rand_vec_boxed();
        for (uint i = 0; i < NDIM; i++) {
            v[i] *= boxsize[i];
        }
        return diff(v, Vec::Zero());
    };
    Vec box_shape() { return boxsize; };

    //! Apply a pure shear of epsilon to reshape box.
    void pure_shear_to(flt epsilon);
    //! Apply a pure shear of epsilon, moving atoms.
    void pure_shear_to(flt epsilon, AtomGroup& atoms);
};

//! Lees-Edwards boundary conditions, with shear in the x-direction, relative to
//! y.
/*!
@ingroup boxes
*/
class LeesEdwardsBox : public OriginBox {
   protected:
    flt gamma;

   public:
    LeesEdwardsBox(Vec size, flt gamma = 0.0) : OriginBox(size), gamma(gamma){};
    LeesEdwardsBox(flt L, flt gamma = 0.0) : OriginBox(L), gamma(gamma){};
    Vec diff(Vec r1, Vec r2);
    virtual Vec diff(Vec r1, Vec r2, boost::array<int, NDIM> boxes);
    virtual boost::array<int, NDIM> box_round(Vec r1, Vec r2);

    //! The current shear amount.
    flt get_gamma() { return gamma; };

    //! Change the shear by dgamma.
    void shear(flt dgamma);
    //! Change the shear by dgamma, and move the atoms as necessary.
    void shear(flt dgamma, AtomGroup& atoms);
    //! Change the shear to gamma.
    void shear_to(flt gamma);
    //! Change the shear to gamma, and move the atoms as necessary.
    void shear_to(flt gamma, AtomGroup& atoms);
    Vec non_affine(Vec v) {
        v[0] -= gamma * v[1];
        return v;
    }

    Vec affine(Vec v) {
        v[0] += gamma * v[1];
        return v;
    }
};

//! A spheocylinder box, also known as a capsule.
/*!
@ingroup boxes
Note that this *does not* keep particles inside the box; use an Interaction like
SCBoxed for that.

This class is useful for functions like V(), dist, edge_dist, inside, rand_loc,
etc.

The spherocylinder has an axis along the x-axis, centered at origin

L is length of the central axis of the cylinder from 1 sphere center to the
other, so L=0 is a sphere.
*/
class SCBox : public Box {
   protected:
    flt L, R;

   public:
    SCBox(flt L, flt R);
    Vec diff(Vec r1, Vec r2) { return r1 - r2; };
    flt V();
    Vec dist(Vec r1);
    Vec edge_dist(Vec r1);
    bool inside(Vec r1, flt buffer = 0.0);
    Vec rand_loc(flt min_dist_to_wall = 0.0);
    flt length() { return L; };
    flt radius() { return R; };
};

////////////////////////////////////////////////////////////////////////////////
//! The basic class for representing each particle.
/*!
@ingroup atoms
Normally instantiated through AtomVec.
*/
struct Atom {
    //! location.
    Vec x;

    //! velocity
    Vec v;

    //! acceleration
    Vec a;

    //! forces
    Vec f;

    //! mass
    flt m;
};

//! A pointer to an Atom.
class AtomRef {
   private:
    Atom* ptr;

   public:
    inline AtomRef() : ptr(NULL){};
    inline AtomRef(Atom* a) : ptr(a){};
    inline Atom& operator*() const { return *ptr; };
    inline Atom* operator->() const { return ptr; }
    inline bool operator==(const AtomRef& other) const {
        return other.ptr == ptr;
    };
    inline bool operator==(const Atom* other) const { return other == ptr; };
    inline bool operator!=(const AtomRef& other) const {
        return other.ptr != ptr;
    };
    inline bool operator<(const AtomRef& other) const {
        return ptr < other.ptr;
    };
    inline bool operator<=(const AtomRef& other) const {
        return ptr <= other.ptr;
    };
    inline bool operator>=(const AtomRef& other) const {
        return ptr >= other.ptr;
    };
    inline bool operator>(const AtomRef& other) const {
        return ptr > other.ptr;
    };
    inline bool is_null() { return ptr == NULL; };
};

//! A pointer to an Atom, that also knows its own index in an AtomVec.
/*!
@ingroup atoms
This is used in many Interaction classes to compair atoms.
*/
class AtomID : public AtomRef {
   private:
    uint num;  // note that these are generally only in reference to
               // a specific AtomGroup
   public:
    //! Constructor not recommended: uses UINT_MAX for number, i.e. unknown
    inline AtomID() : AtomRef(), num(UINT_MAX){};
    //! Recommended constructor
    inline AtomID(Atom* a, uint n) : AtomRef(a), num(n){};
    inline uint n() const { return num; };
};

class IDPair {
   private:
    AtomID id1, id2;

   public:
    IDPair() : id1(), id2(){};
    IDPair(AtomID a, AtomID b) : id1(a), id2(b){};
    inline AtomID first() const { return id1; };
    inline AtomID last() const { return id2; };
};

class AtomGroup;

//! For iterating through an AtomGroup.
class AtomIter {
   private:
    uint i;
    AtomGroup& g;

   public:
    AtomIter(AtomGroup& g, uint i) : i(i), g(g){};
    bool operator!=(const AtomIter& other) const { return i != other.i; };
    Atom& operator*() const;
    inline const AtomIter& operator++() {
        ++i;
        return *this;
    };
};

class AtomVec;

//! a group of atoms, such as all of them (AtomVec), or a smaller group such as
//! a molecule, sidebranch, etc.
/*!
@ingroup atoms
*/
class AtomGroup {
   public:
    // access individual atoms
    virtual AtomVec& vec() = 0;
    virtual Atom& operator[](cuint n) = 0;
    virtual Atom& operator[](cuint n) const = 0;
    virtual Atom& get(cuint n) { return ((*this)[n]); };
    virtual AtomID get_id(cuint n) = 0;
    //! Number of atoms in the group
    virtual uint size() const = 0;
    //! For use in a for loop
    virtual AtomIter begin() { return AtomIter(*this, 0); };
    virtual AtomIter end() { return AtomIter(*this, (uint)size()); };

    //! center of mass
    Vec com() const;
    //! center of mass force (i.e., sum of all forces)
    Vec com_force() const;
    //! center of mass velocity
    Vec com_velocity() const;

    //! Mass of the whole group
    flt mass() const;
    //! Total kinetic_energy energy of the group.
    /*!
    This is normally with reference to a "lab" reference frame (velocity
    (0,0,0)), but
    can optionally take a different origin velocity, e.g. `com_velocity()`.
    */
    flt kinetic_energy(const Vec originvelocity = Vec::Zero()) const;
    //! Total momentum.
    Vec momentum() const;
    //! \f$R_g\f$
    flt gyradius() const;
#ifdef VEC3D
    //! Total torque about a given location
    Vec torque(const Vec loc) const;
    Vec torque() const { return torque(com()); };
    //! Moment of inertia of the atoms about an axis
    flt moment_about(const Vec axis, const Vec loc) const;
    flt moment_about(const Vec axis) const {
        return moment_about(axis, com());
    };
    //! Angular momentum
    Vec angular_momentum(const Vec loc) const;
    Vec angular_momentum() const { return angular_momentum(com()); };
    //! Moment of inertia of the atoms as a whole
    Matrix moment(const Vec loc) const;
    Matrix moment() const { return moment(com()); };
    //! Angular velocity
    Vec omega(const Vec loc) const;
    Vec omega() const { return omega(com()); };
    //! Add a given angular velocity to all atoms, by adding to their velocity
    void add_omega(Vec w, Vec origin);
    void add_omega(Vec w) { return add_omega(w, com()); };
    //! Reset angular momentum to 0
    inline void reset_L() {
        Vec c = com(), w = omega(c);
        if (w.squaredNorm() == 0) return;
        add_omega(-w, c);
    }
#elif defined VEC2D
    //! Total torque about a given location
    flt torque(const Vec loc) const;
    //! Moment of inertia of the atoms as a whole
    flt moment(const Vec loc) const;
    flt moment() const { return moment(com()); };
    //! Angular momentum
    flt angular_momentum(const Vec loc) const;
    flt angular_momentum() const { return angular_momentum(com()); };
    //! Angular velocity
    flt omega(const Vec loc) const {
        return angular_momentum(loc) / moment(loc);
    };
    flt omega() const { return omega(com()); };
    //! Add a given angular velocity to all atoms, by adding to their velocity
    void add_omega(flt w, Vec origin);
    void add_omega(flt w) { add_omega(w, com()); }
    //! Reset angular momentum to 0
    inline void reset_L() {
        Vec c = com();
        flt w = omega(c);
        if (w == 0) return;
        add_omega(-w, c);
    }
#endif

    //! for resetting. Adds a fixed velocity to all atoms
    void add_velocity(Vec v);
    //! Subtracts the center of mass velocity from all atoms
    void reset_com_velocity() { add_velocity(-com_velocity()); };
    //! Randomize velocities, for a specific temperature
    void randomize_velocities(flt T);

    //! for timestepping
    void reset_forces();
    virtual ~AtomGroup(){};
};

inline Atom& AtomIter::operator*() const { return g[i]; };

/*!
@ingroup basics atoms
@brief The main class for representing particles.
*/
class AtomVec : public virtual AtomGroup {
   private:
    Atom* atoms;
    uint sz;
    void zero() {
        for (uint i = 0; i < sz; i++) {
            atoms[i].x = Vec::Zero();
            atoms[i].v = Vec::Zero();
            atoms[i].f = Vec::Zero();
            atoms[i].a = Vec::Zero();
        }
    }

   public:
    AtomVec(vector<double> masses) : sz((uint)masses.size()) {
        atoms = new Atom[sz];
        for (uint i = 0; i < sz; i++) atoms[i].m = masses[i];
        zero();
    };
    AtomVec(uint N, flt mass) : sz(N) {
        atoms = new Atom[sz];
        for (uint i = 0; i < sz; i++) atoms[i].m = mass;
        zero();
    };
    AtomVec(AtomVec& other) : sz(other.size()) {
        atoms = new Atom[sz];
        for (uint i = 0; i < sz; i++) atoms[i] = other.atoms[i];
    };
    AtomVec& vec() { return *this; };
    inline Atom& operator[](cuint n) { return atoms[n]; };
    inline Atom& operator[](cuint n) const { return atoms[n]; };
    inline AtomID get_id(cuint n) {
        if (n > sz) return AtomID();
        return AtomID(atoms + n, n);
    };
    inline uint size() const { return sz; };

    ~AtomVec() { delete[] atoms; };
};

/*!
A class for representing any grouping of atoms that is not the whole set of
atoms, such as a molecule, a side-chain, etc.
@ingroup atoms
*/
class SubGroup : public AtomGroup {
   protected:
    sptr<AtomVec> atoms;
    vector<AtomID> ids;

   public:
    SubGroup(sptr<AtomVec> atoms) : atoms(atoms){};
    AtomVec& vec() { return *atoms; };
    inline Atom& operator[](cuint n) { return *ids[n]; };
    inline Atom& operator[](cuint n) const { return *ids[n]; };
    inline Atom& get(cuint n) { return *ids[n]; };
    inline void add(AtomID a) {
        std::vector<AtomID>::iterator it = std::find(ids.begin(), ids.end(), a);
        if (it != ids.end())
            throw std::invalid_argument(
                "Cannot add AtomID to SubGroup: it already exists.");
        return ids.push_back(a);
    };
    inline AtomID get_id(cuint n) { return ids[n]; };
    inline uint size() const { return (uint)ids.size(); };
};

#endif
