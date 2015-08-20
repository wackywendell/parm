#include "trackers.hpp"

#ifndef INTERACTION_H
#define INTERACTION_H

#define sptr boost::shared_ptr
using namespace boost; // required for SWIG for some reason

/***********************************************************************
 * L-J Note:
 * The convention here is ε(σ⁶/r⁶ - 1)², which has its minimum at r = σ.
 * This is what to use with Alice's parameters.
 */

/*
 * New plan of attack:
 *
 * Interaction<N> - blindly takes N vectors,
 * and calculates energy or forces.
 *
 * AtomGroup(N) (or (N,M) with array?)-
 * has a list of atoms, can tell mass, vel, loc.
 * keeps forces?
 * possibly has M groups of groups?
 * iterator over atoms necessary
 * perhaps class iterator over all of them?
 * perhaps class iterator over all atoms?
 * has "uint size()" function, says number of atoms
 *
 * InteractGroup - has a single Interaction<N>, and pointers to however
 * many AtomGroups it needs, and returns energy and/or sets forces.
 * one example: neighbors; yields pairs of i,i+1.
 * second example: neighbor list, keeps track of nearby atoms.
 * perhaps class iterator over all of them?
 *
 * perhaps split: interactiongroup has a group, and an Interaction, and
 * applies Interaction to each pair in each group.
 * Not incompatible.
 *
 * integrator - keeps a list of Interaction groups, can get energy.
 * Also goes through list to integrate over time.
 *
 * statistic - keeps a list of Interaction groups, derives a statistic
 * from them.
 *
 * statkeeper - keeps a list of statistics, can write to file?
 *
 * constants - keeps track of constants; in map, or as properties?
 * perhaps a struct.
*/

//typedef double flt; defined in vecrand.hpp
//typedef unsigned int uint;
/**
The basic Interaction class, used to represent a potential function. Specific interactions should
derive from this.
*/
class Interaction : public boost::enable_shared_from_this<Interaction> {
    public:
        /**
        Potential energy due to this Interaction.
        */
        virtual flt energy(Box &box)=0;
        virtual void setForces(Box &box)=0;
        /**
        Set forces (`Atom.f`) and return \f$P = \sum_{\left<i,j \right>} \vec r_{ij} \cdot \vec F_{ij}\f$
        at the same time (see `pressure()`).
        */
        virtual flt setForcesGetPressure(Box &box){return NAN;};
        /**
        Partial pressure due to this Interaction.

        \f$P = \sum_{\left<i,j \right>} \vec r_{ij} \cdot \vec F_{ij}\f$, or equivalently
        \f$P = \sum_i \vec r_i \cdot \vec F_i\f$
        Note that the full pressure involves *all* interactions and temperature
        */
        virtual flt pressure(Box &box)=0;
        virtual ~Interaction(){};
};

/***********************************************************************
 * Interaction Basics
 */
class InteractPair {
    public:
        virtual flt energy(const Vec diff)=0;
        virtual Vec forces(const Vec diff)=0;
        virtual ~InteractPair(){};
};

static const flt LJr0 = pow(2.0, 1.0/6.0);
static const flt LJr0sq = pow(2.0, 1.0/3.0);

class LJRepulsive {
    protected:
        flt epsilon;
        flt sigma;
    public:
        LJRepulsive(const flt epsilon, const flt sigma):
            epsilon(epsilon), sigma(sigma){};
        inline static flt energy(const Vec diff, const flt eps, const flt sig){
            flt rsq = diff.squaredNorm()/(sig*sig);
            if(rsq > 1) return 0;
            flt rsix = rsq*rsq*rsq;
            //~ return eps*(4*(1/(rsix*rsix) - 1/rsix) + 1);
            flt mid = (1-1/rsix);
            return eps*(mid*mid);
        };
        inline flt energy(const Vec& diff){return energy(diff, epsilon, sigma);};
        inline static Vec forces(const Vec diff, const flt eps, const flt sig){
            flt dsq = diff.squaredNorm();
            flt rsq = dsq/(sig*sig);
            if(rsq > 1) return Vec::Zero();
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            //~ flt fmagTimesR = eps*(4*(12/(rsix*rsix) - 6/rsix));
            flt fmagTimesR = 12*eps/rsix*(1/rsix - 1);
            //~ cout << "Repulsing " << diff << "with force"
                 //~ << diff * (fmagTimesR / dsq) << '\n';
            //~ cout << "mag: " << diff.norm() << " sig:" << sig << " eps:" << eps
                 //~ << "F: " << (diff * (fmagTimesR / dsq)).norm() << '\n';
            //~ cout << "E: " << energy(diff, sig, eps) << " LJr0: " << LJr0 << ',' << LJr0sq << "\n";
            return diff * (fmagTimesR / dsq);
        };
        inline Vec forces(const Vec& diff){return forces(diff, epsilon, sigma);};
};

class LJAttract {
    protected:
        flt epsilon;
        flt sigma;
    public:
        LJAttract(const flt epsilon, const flt sigma):
            epsilon(epsilon), sigma(sigma){};
        inline static flt energy(const Vec diff, const flt eps, const flt sig){
            flt rsq = diff.squaredNorm()/(sig*sig);
            if(rsq < 1) return -eps;
            flt rsix = rsq*rsq*rsq;
            //~ return eps*(4*(1/(rsix*rsix) - 1/rsix) + 1);
            flt mid = (1-1/rsix);
            return eps*(mid*mid-1);
        };
        inline static flt energy(const flt rsig){
            if(rsig < 1) return -1;
            flt rsq = rsig * rsig;
            flt rsix = rsq*rsq*rsq;
            flt mid = (1-1/rsix);
            return mid*mid-1;
            // Minimum of -1 occurs at r = 1
        };
        inline flt energy(const Vec& diff){return energy(diff, epsilon, sigma);};
        inline static Vec forces(const Vec diff, const flt eps, const flt sig){
            flt dsq = diff.squaredNorm();
            flt rsq = dsq/(sig*sig);
            if(rsq < 1) return Vec::Zero();
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            flt fmagTimesR = 12*eps/rsix*(1/rsix - 1);
            return diff * (fmagTimesR / dsq);
        };
        inline static flt forces(const flt rsig){
            if(rsig < 1) return 0;
            flt rsq = rsig*rsig;
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            flt fmagTimesR = 12/rsix*(1/rsix - 1);
            return fmagTimesR / rsig;
        };
        inline Vec forces(const Vec& diff){return forces(diff, epsilon, sigma);};
};

class LJAttractCut {
    // Purely attractive
    protected:
        flt epsilon;
        flt sigma;
        flt cutR, cutE;
    public:
        inline LJAttractCut(const flt epsilon, const flt sigma, const flt cutsig):
            epsilon(epsilon), sigma(sigma), cutR(cutsig),
            cutE(LJAttract::energy(cutR) * epsilon){};
        inline static flt energy(const Vec diff, const flt eps,
                                    const flt sig, const flt cutsig){
            if(eps == 0) return 0;
            if(diff.squaredNorm() > (cutsig*cutsig*sig*sig)) return 0;
            return (LJAttract::energy(diff, eps, sig) -
                eps*LJAttract::energy(cutsig));
        };
        inline flt energy(const Vec& diff){
            if(epsilon == 0) return 0;
            if(diff.squaredNorm() > (cutR*cutR*sigma*sigma)) return 0;
            return LJAttract::energy(diff, epsilon, sigma) - cutE;
        };
        inline static Vec forces(const Vec diff, const flt eps,
                                    const flt sig, const flt cutsig){
            if(eps == 0) return Vec::Zero();
            flt dsq = diff.squaredNorm();
            flt rsq = dsq/(sig*sig);
            if(rsq < 1 or rsq > (cutsig*cutsig)) return Vec::Zero();
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            flt fmagTimesR = 12*eps/rsix*(1/rsix - 1);
            return diff * (fmagTimesR / dsq);
        };
        inline Vec forces(const Vec& diff){
                        return forces(diff, epsilon, sigma, cutR);};
};

class LJFullCut {
    // uses ε*((1-1/r⁶)² - 1)
    // sigma is thus at the minimum
    // cutR is unitless; that times sigma is the distance we cut the potential
    protected:
        flt epsilon;
        flt sigma;
        flt cutR, cutE;
    public:
        LJFullCut(const flt epsilon, const flt sigma, const flt cutsig):
            epsilon(epsilon), sigma(sigma), cutR(cutsig){
                flt rsix = pow(cutR, 6);
                flt mid = (1-1/rsix);
                cutE = epsilon*(mid*mid-1);
            };
        inline flt energy(const Vec& diff){
            flt rsq = diff.squaredNorm()/(sigma*sigma);
            if(rsq > (cutR*cutR)) return 0;
            flt rsix = rsq*rsq*rsq;
            flt mid = (1-1/rsix);
            return epsilon*(mid*mid-1) - cutE;
        };
        inline static Vec forces(const Vec diff, const flt eps, const flt sig, const flt cutsig){
            flt dsq = diff.squaredNorm();
            flt rsq = dsq/(sig*sig);
            if(rsq > (cutsig*cutsig)) return Vec::Zero();
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            flt fmagTimesR = 12*eps/rsix*(1/rsix - 1);
            return diff * (fmagTimesR / dsq);
        };
        inline static flt forces(const flt rsig, const flt cutsig){
            if(rsig > cutsig) return 0;
            flt rsq = rsig*rsig;
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            flt fmagTimesR = 12/rsix*(1/rsix - 1);
            return fmagTimesR / rsig;
        };
        inline Vec forces(const Vec& diff){return forces(diff, epsilon, sigma, cutR);};
};

class Spring : public InteractPair {
    protected:
        flt springk;
        flt x0;
    public:
        Spring(const flt k, const flt x0) : springk(k),x0(x0){};
        flt energy(const Vec diff);
        Vec forces(const Vec diff);
        ~Spring(){};
};

class BondAngle {
    // two vectors for E and F are from the central one to the other 2
    protected:
        flt springk;
        flt theta0;
        bool usecos;

    public:
        BondAngle(const flt k, const flt theta, const bool cosine=false)
                        :springk(k), theta0(theta), usecos(cosine){};
        inline static flt get_angle(const Vec& r1, const Vec& r2){
            flt costheta = r1.dot(r2) / r1.norm() / r2.norm();
            if(costheta > 1) costheta = 1;
            else if(costheta < -1) costheta = -1;
            return acos(costheta);
        };
        flt energy(const Vec& diff1, const Vec& diff2);
        array<Vec,3> forces(const Vec& diff1, const Vec& diff2);
        ~BondAngle(){};
};

#ifdef VEC3D
struct DihedralDerivs {
    array<Vec,4> derivs;
    flt costheta;
};

class Dihedral {
    // vectors are from 1 to 2, 2 to 3, 3 to 4
    protected:
        vector<flt> coscoeffs;
        vector<flt> sincoeffs;
        bool usepow;
        flt dudcosthetaCOS(const flt costheta) const;
    public:
        Dihedral(const vector<flt> cosvals,
                const vector<flt> sinvals = vector<flt>(),
                bool usepow = true);
        static flt getcos(const Vec& diff1, const Vec& diff2, const Vec& diff3);
        static flt getang(const Vec& diff1, const Vec& diff2, const Vec& diff3);

        /// Returns dr_i / d cos(θ)
        static DihedralDerivs dr_dcostheta(const Vec& diff1, const Vec& diff2, const Vec& diff3);

        flt dudcostheta(const flt theta) const;

        inline flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const {
            return energy(getang(diff1, diff2, diff3));
        };
        flt energy(flt ang) const;
        array<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
};
#endif

class ElectricScreened : public InteractPair {
    protected:
        flt screen;
        flt q1, q2;
        flt cutoff;
        flt cutoffE;
    public:
        ElectricScreened(const flt screenLength, const flt q1,
            const flt q2, const flt cutoff);
        inline flt energy(const Vec r){return energy(r.norm(),q1*q2,screen, 0) - cutoffE;};
        static flt energy(const flt r, const flt qaqb, const flt screen, const flt cutoff=0);
        inline Vec forces(const Vec r){return forces(r,q1*q2,screen, cutoff);};
        static Vec forces(const Vec r, const flt qaqb, const flt screen, const flt cutoff=0);
};

//~ class InteractGroup : public Interaction {
    //~ protected:
        //~ vector<Interaction*> inters;
    //~ public:
        //~ InteractGroup(vector<Interaction*> inters=vector<Interaction*>())
                    //~ : inters(inters){};
        //~ void add(Interaction* a){inters.push_back(a);};
        //~ uint size() const{ return inters.size();};
        //~ flt energy(Box *box);
        //~ void setForces(Box *box);
//~ };

struct FixedForceAtom {
    Vec F;
    AtomID a;
    FixedForceAtom(Vec F, AtomID a) : F(F), a(a) {};
    flt energy(Box &box){return -F.dot(a->x);};
    void setForce(Box &box){a->f += F;};
};

class FixedForce : public Interaction {
    protected:
        vector<FixedForceAtom> atoms;
    public:
        FixedForce(vector<FixedForceAtom> atoms = vector<FixedForceAtom>()) : atoms(atoms){};
        void add(FixedForceAtom a){atoms.push_back(a);};
        void add(Vec F, AtomID a){add(FixedForceAtom(F,a));};
        #ifdef VEC3D
        void add(flt x, flt y, flt z, AtomID a){add(FixedForceAtom(Vec(x,y,z),a));};
        #elif defined VEC2D
        void add(flt x, flt y, AtomID a){add(FixedForceAtom(Vec(x,y),a));};
        #endif
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box &box){
            flt E=0;
            for(vector<FixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<FixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};


struct FixedForceRegionAtom : public AtomID {
    Vec direction;          // will be normalized
    vector<flt> boundaries; // should be 1 smaller than Fs
                            // each boundary is at (direction * b),
                            // where b is in boundaries
    vector<flt> Fs;
    FixedForceRegionAtom(AtomID a, Vec direction, vector<flt> boundaries, vector<flt> Fs);
    flt energy(Box &box);
    void setForce(Box &box);
};

class FixedForceRegion : public Interaction {

    protected:
        vector<FixedForceRegionAtom> atoms;
    public:
        FixedForceRegion(vector<FixedForceRegionAtom> atoms = vector<FixedForceRegionAtom>())
            : atoms(atoms){};

        void add(FixedForceRegionAtom a){atoms.push_back(a);};
        void add(AtomID a, Vec dir, vector<flt> bound, vector<flt> F){
            add(FixedForceRegionAtom(a, dir, bound, F));};
        uint size() const{ return (uint) atoms.size();};

        flt energy(Box &box){
            flt E=0;
            for(vector<FixedForceRegionAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<FixedForceRegionAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};

struct FixedSpringAtom {
    Vec loc;
    flt k;
    bool usecoord[3];
    AtomID a;
    FixedSpringAtom(AtomID a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true) :
            loc(loc), k(k), a(a) {
                usecoord[0] = usex;
                usecoord[1] = usey;
                usecoord[2] = usez;
                };
    flt energy(Box &box){
            Vec diffx = a->x - loc;
            for(uint i=0; i<2; ++i){
                if(!usecoord[i]) diffx[i] = 0;
            }
            return k*diffx.squaredNorm()/2;
        };
    void setForce(Box &box){
        Vec diffx = a->x - loc;
        for(uint i=0; i<2; ++i){
            if(!usecoord[i]) diffx[i] = 0;
        }
        a->f -= diffx * k;
    };
};

class FixedSpring : public Interaction{
    protected:
        vector<FixedSpringAtom> atoms;
    public:
        FixedSpring(vector<FixedSpringAtom> atoms = vector<FixedSpringAtom>()) : atoms(atoms){};
        void add(FixedSpringAtom a){atoms.push_back(a);};
        void add(AtomID a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true){
            add(FixedSpringAtom(a, loc, k, usex, usey, usez));};
        //void add(Vec F, AtomID a){add(FixedForceAtom(F,a));};
        //void add(flt x, flt y, flt z, AtomID a){add(FixedForceAtom(Vec(x,y,z),a));};
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box &box){
            flt E=0;
            for(vector<FixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<FixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};

class COMSpring : public Interaction{
    protected:
        AtomGroup *g1;
        AtomGroup *g2;
        flt k, x0;
        flt m1, m2;
    public:
        COMSpring(AtomGroup *g1, AtomGroup *g2, flt k, flt x0=0) :
            g1(g1), g2(g2), k(k), x0(x0), m1(g1->mass()), m2(g2->mass()){};
        flt energy(Box &box){
            flt dx = (g1->com() - g2->com()).norm() - x0;
            return k/2 * dx * dx;
        };
        void setForces(Box &box){
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.norm();
            flt fmag = -k * (comdist - x0);
            Vec a1 = comvec * (fmag / m1 / comdist);
            for(uint i=0; i < g1->size(); ++i){
                Atom &atm = g1->get(i);
                atm.f += a1 * atm.m;
            }

            Vec a2 = comvec * (-fmag / m2 / comdist);
            for(uint i=0; i < g2->size(); ++i){
                Atom &atm = g2->get(i);
                atm.f += a2 * atm.m;
            }
        };
        flt pressure(Box &box){
            //~ return NAN;
            // I think this is right, but I haven't checked it
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.norm();
            flt fmag = -k * (comdist - x0);

            Vec a12 = comvec * (fmag / m1 / m2 / comdist);

            flt P = 0;
            for(uint i=0; i < g1->size(); ++i){
                for(uint j=0; j < g2->size(); ++j){
                    Atom &atm1 = g1->get(i);
                    Atom &atm2 = g2->get(i);
                    Vec fij = a12 * (atm1.m * atm2.m);
                    Vec rij = box.diff(atm1.x, atm2.x);
                    P += fij.dot(rij);
                }
            }
            return P;
        };
};

//////////////////////////////////////////////////////////////////////////////
// Random Force

enum RandomForceType {
    FIXED, // Always the same magnitude
    UNIFORM, // uniform probability distribution for magnitude, from 0 to force_mag
    GAUSSIAN // Gaussian probability distribution for magnitude
};

struct RandomForceAtom : public AtomRef {
    public:
        flt force_mag;
        flt freq; // in general, how many timesteps on average between kicks
        RandomForceType force_type;
    public:
        RandomForceAtom(AtomID a, flt force_mag, flt freq, RandomForceType force_type=UNIFORM) : 
            AtomRef(a), force_mag(force_mag), freq(freq), force_type(force_type){};
};

class RandomForce : public Interaction {
    public:
        vector<RandomForceAtom> group;
    
    public:
        RandomForce(){};
        RandomForce(AtomGroup& agroup, flt force_mag, flt freq, RandomForceType force_type=UNIFORM){
            for(uint i=0; i<agroup.size(); ++i){
                group.push_back(RandomForceAtom(agroup.get_id(i), force_mag, freq, force_type));
            }
        };
        
        uint size() const{ return (uint) group.size();};
        RandomForceAtom get(uint i) const{ return group[i];};
        
        /// If "replace", a previous pair found will be replaced by the new pair.
        /// If not "replace" and that pair of atoms is already inserted, an error will be thrown.
        bool add(RandomForceAtom a, bool replace=true);
        
        // No potential energy
        flt energy(Box &box){return 0;};
        void setForces(Box &box);
        
        // no pressure from this interaction
        //TODO: maybe this should be average pressure
        flt pressure(Box &box){return 0;};
        
        // no pressure from this interaction
        //TODO: maybe this should actually work; it could
        flt setForcesGetPressure(Box &box){setForces(box); return 0;};
        
};

//////////////////////////////////////////////////////////////////////////////
// Bonds

enum BondDiffType {
    BOXED, // use box.diff(r1, r2)
    UNBOXED, // use r2 - r1
    FIXEDBOX // use r2 - r1 - (original box separation)
};

struct BondGrouping {
    flt k, x0;
    AtomID a1, a2;
    BondDiffType diff_type;
    array<int,NDIM> fixed_box;
    BondGrouping(flt k, flt x0, AtomID a1, AtomID a2,
            BondDiffType diff=UNBOXED, OriginBox *box=NULL);
    Vec diff(Box &box) const;
    int get_fixed(uint i){return fixed_box[i];};
    inline bool same_atoms(BondGrouping &other){
        return ((a1 == other.a1) and (a2 == other.a2)) or ((a1 == other.a2) and (a2 == other.a1));
    };
};

class BondPairs : public Interaction {
    protected:
        bool zeropressure;
        vector<BondGrouping> pairs;
        //inline static Vec diff(Box &box, Vec r1, Vec r2){return r1-r2;};
        //inline static Vec diff(Box &box, Vec r1, Vec r2){return box.diff(r1, r2);};
    public:
        BondPairs(vector<BondGrouping> pairs, bool zeropressure=true);
        BondPairs(bool zeropressure=true);
        /// Add a pair of atoms.
        /// If "replace", a previous pair found will be replaced by the new pair.
        /// If not "replace" and that pair of atoms is already inserted, an error will be thrown.
        bool add(BondGrouping b, bool replace=true);
        inline bool add(flt k, flt x0, AtomID a1, AtomID a2, bool replace=true){
            return add(BondGrouping(k,x0,a1,a2), replace);};
        void add_forced(BondGrouping b){pairs.push_back(b);};
        /// Add a pair of atoms with the current distance.
        inline bool add(flt k, AtomID a1, AtomID a2, bool replace=true){
            flt x0 = (a1->x - a2->x).norm();
            return add(BondGrouping(k,x0,a1,a2), replace);
        };

        uint size() const{ return (uint) pairs.size();};
        BondGrouping get(uint i) const{ return pairs[i];};
        flt mean_dists(Box &box) const;
        flt std_dists(Box &box) const;
        flt energy(Box &box);
        void setForces(Box &box);
        flt pressure(Box &box);
        flt setForcesGetPressure(Box &box);
};

struct AngleGrouping {
    flt k, x0;
    AtomID a1, a2, a3;
    AngleGrouping(flt k, flt x0, AtomID a1, AtomID a2, AtomID a3) :
                k(k),x0(x0), a1(a1), a2(a2), a3(a3){};
    inline bool same_atoms(AngleGrouping &other){
        if(a2 != other.a2) return false;
        return ((a1 == other.a1) and (a3 == other.a3)) or ((a1 == other.a3) and (a3 == other.a1));
    };
};

class AngleTriples : public Interaction {
    protected:
        vector<AngleGrouping> triples;
        inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    public:
        AngleTriples(vector<AngleGrouping> triples = vector<AngleGrouping>());
        /// Add a triple of atoms.
        /// If "replace", a previous triple found will be replaced by the new triple.
        /// If not "replace" and that triple of atoms is already inserted, an error will be thrown.
        bool add(AngleGrouping b, bool replace=true);
        inline bool add(flt k, flt x0, AtomID a1, AtomID a2, AtomID a3, bool replace=true){
            return add(AngleGrouping(k,x0,a1,a2,a3), replace);};
        /// Add a triple of atoms with the current angle.
        bool add(flt k, AtomID a1, AtomID a2, AtomID a3, bool replace=true);
        void add_forced(AngleGrouping b){triples.push_back(b);};

        flt energy(Box &box);
        inline flt pressure(Box &box){return 0;};
        void setForces(Box &box);
        inline flt setForcesGetPressure(Box &box){setForces(box); return 0;};
        uint size() const {return (uint) triples.size();};
        flt mean_dists() const;
        flt std_dists() const;
};

#ifdef VEC3D
struct DihedralGrouping {
    inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    Dihedral dih;
    AtomID a1, a2, a3, a4;
    DihedralGrouping(vector<flt> coscoeffs, vector<flt> sincoeffs,
                AtomID a1, AtomID a2, AtomID a3, AtomID a4, bool usepow=true) :
                dih(coscoeffs, sincoeffs, usepow), a1(a1),
                a2(a2), a3(a3), a4(a4){};
};

class Dihedrals : public Interaction {
    protected:
        vector<DihedralGrouping> groups;
    public:
        Dihedrals(vector<DihedralGrouping> pairs = vector<DihedralGrouping>());
        void add(DihedralGrouping b){groups.push_back(b);};
        inline void add(vector<flt> nums, AtomID a1, AtomID a2, AtomID a3, AtomID a4){
            add(DihedralGrouping(nums, vector<flt>(), a1,a2,a3,a4));};
        inline void add(vector<flt> coscoeffs, vector<flt> sincoeffs,
                            AtomID a1, AtomID a2, AtomID a3, AtomID a4, bool usepow=true){
            add(DihedralGrouping(coscoeffs, sincoeffs,a1,a2,a3,a4, usepow));};
        /** Add 4 atoms with the potential \f$V(\theta) = k (1 + \cos(\theta - \theta_0))\f$. */
        inline void add(flt k, flt theta0, AtomID a1, AtomID a2, AtomID a3, AtomID a4){
            vector<flt> coscoeffs(2,k);
            coscoeffs[1] = -k*cos(theta0);
            vector<flt> sincoeffs(2,0);
            sincoeffs[1] = -k*sin(theta0);
            add(DihedralGrouping(coscoeffs, sincoeffs, a1,a2,a3,a4, true));
        }
        /**
        Add 4 atoms with the potential \f$V(\theta) = k (1 + \cos(\theta - \theta_0))\f$, where
        \f$\theta_0\f$ is set to match the current positions of the atoms.
        */
        inline void add(flt k, AtomID a1, AtomID a2, AtomID a3, AtomID a4){
            Vec r1 = DihedralGrouping::diff(a2->x, a1->x);
            Vec r2 = DihedralGrouping::diff(a3->x, a2->x);
            Vec r3 = DihedralGrouping::diff(a4->x, a3->x);
            add(k,  Dihedral::getang(r1, r2, r3), a1, a2, a3, a4);
        };
        uint size() const{ return uint(groups.size());};
        flt mean_dists() const;
        //~ flt std_dists() const;
        flt energy(Box &box);
        void setForces(Box &box);
        inline flt pressure(Box &box){return 0;};
        inline flt setForcesGetPressure(Box &box){setForces(box); return 0;};
};
#endif


////////////////////////////////////////////////////////////////////////
struct ForcePair {
    Atom a1, a2;
    Vec fij;
};

struct ForcePairX {
    AtomID a1, a2;
    flt xij;
    Vec fij;
};

//
class FPairXFunct {
    public:
        virtual void run(ForcePairX*) = 0;
        virtual ~FPairXFunct() = 0;
};

inline FPairXFunct::~FPairXFunct() { }  // defined even though it's pure virtual; it's faster this way; trust me

class InteractionPairsX : public Interaction {
    public:
        using Interaction::setForces;
        virtual void setForces(Box &box, FPairXFunct*)=0;
        virtual ~InteractionPairsX(){};
};

struct Charged : public AtomID {
    flt q;
    Charged() : AtomID(), q(0){};
    Charged(flt q, AtomID a) : AtomID(a), q(q){};
};

struct ChargePair {
    flt q1q2;
    AtomID atom1, atom2;
    ChargePair(Charged a1, Charged a2) : q1q2(a1.q*a2.q){};
};

////////////////////////////////////////////////////////////////////////
// Repulsive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// Potential is V(r) = ε (σ⁶/r⁶ - 1)²
// cutoff at sigma
struct LJatom : public AtomID {
    flt epsilon, sigma;
    LJatom(){};
    LJatom(flt epsilon, flt sigma, AtomID a) : AtomID(a),
            epsilon(epsilon), sigma(sigma){};
    LJatom(AtomID a, LJatom other) : AtomID(a),
                    epsilon(other.epsilon), sigma(other.sigma){};
    flt maxsize(){return sigma;};
};

struct LJpair {
    flt epsilon, sigma;
    AtomID atom1, atom2;
    LJpair(LJatom LJ1, LJatom LJ2) :
            epsilon(sqrt(LJ1.epsilon * LJ2.epsilon)),
            sigma((LJ1.sigma + LJ2.sigma) / 2),
            atom1(LJ1), atom2(LJ2){};
    inline flt energy(Box &box){
        return LJRepulsive::energy(box.diff(atom1->x,atom2->x), epsilon, sigma);};
    inline Vec forces(Box &box){
        return LJRepulsive::forces(box.diff(atom1->x,atom2->x), epsilon, sigma);};
};

////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// cutoff at some sigcut
// Minimum at σ
struct LJatomcut : public LJatom {
    flt sigcut; // sigma units
    LJatomcut(){};
    LJatomcut(flt epsilon, flt sigma, AtomID a, flt cut) :
            LJatom(epsilon, sigma, a), sigcut(cut){};
    LJatomcut(AtomID a, LJatomcut other) : LJatom(a, other),
        sigcut(other.sigcut){};
    flt maxsize(){return sigma*sigcut;};
};



////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = ε₁₂ and σ = σ₁₂ (both indexed)
// cutoff at some sigcut (and r < σ)

struct LJAtomIndexed : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    vector<flt> sigmas; // for finding epsilons
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two HydroAtoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]; same for sigma
    flt sigcut; // sigma units
    LJAtomIndexed(){};
    LJAtomIndexed(vector<flt> epsilons, vector<flt> sigmas, uint indx, AtomID a, flt cut) :
            AtomID(a), epsilons(epsilons), sigmas(sigmas), indx(indx),
             sigcut(cut){
                 assert(sigmas.size() == epsilons.size());
            };
    LJAtomIndexed(AtomID a, LJAtomIndexed other) : AtomID(a), epsilons(other.epsilons),
        sigmas(other.sigmas), indx(other.indx), sigcut(other.sigcut){};
    flt getEpsilon(LJAtomIndexed &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        //~ assert(indx < other.epsilons.size());
        //~ assert(other.epsilons[indx] == myeps);
        return myeps;
    }
    flt getSigma(LJAtomIndexed &other){
        assert(other.indx < sigmas.size());
        flt myeps = sigmas[other.indx];
        //~ assert(indx < other.sigmas.size());
        //~ assert(other.sigmas[indx] == myeps);
        return myeps;
    }
    flt maxsize(){
        flt sigma = sigmas[0];
        for(uint i=1; i<sigmas.size(); ++i){
            if(sigma < sigmas[i]) sigma = sigmas[i];
        }
        return sigma*sigcut;
    };
};

struct LJCutPair {
    LJFullCut inter;
    AtomID atom1, atom2;
    LJCutPair(LJatomcut a1, LJatomcut a2) :
        inter(sqrt(a1.epsilon * a2.epsilon),
              (a1.sigma + a2.sigma) / 2,
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    LJCutPair(LJAtomIndexed a1, LJAtomIndexed a2) :
        inter(a1.getEpsilon(a2),
              a1.getSigma(a2),
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};

    inline flt energy(Box &box){return inter.energy(box.diff(atom1->x, atom2->x));};
    inline Vec forces(Box &box){return inter.forces(box.diff(atom1->x, atom2->x));};
};

struct LJAttractPair {
    LJAttractCut inter;
    AtomID atom1, atom2;
    LJAttractPair(LJatomcut a1, LJatomcut a2) :
        inter(sqrt(a1.epsilon * a2.epsilon),
              (a1.sigma + a2.sigma) / 2,
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    LJAttractPair(LJAtomIndexed a1, LJAtomIndexed a2) :
        inter(a1.getEpsilon(a2),
              a1.getSigma(a2),
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){return inter.energy(box.diff(atom1->x, atom2->x));};
    inline Vec forces(Box &box){return inter.forces(box.diff(atom1->x, atom2->x));};
};

////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = ε₁₂ (indexed) and σ = (σ₁+σ₂)/2
// cutoff at some sigcut (and r < σ)

struct HydroAtom : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two HydroAtoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]
    flt sigma;
    flt sigcut; // sigma units
    HydroAtom(){};
    HydroAtom(vector<flt> epsilons, uint indx, flt sigma, AtomID a, flt cut) :
            AtomID(a), epsilons(epsilons), indx(indx), sigma(sigma),
             sigcut(cut){};
    HydroAtom(AtomID a, HydroAtom other) : AtomID(a), epsilons(other.epsilons),
        indx(other.indx), sigma(other.sigma), sigcut(other.sigcut){};
    flt getEpsilon(HydroAtom &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        assert(indx < other.epsilons.size());
        assert(other.epsilons[indx] == myeps);
        return myeps;
    }
    flt maxsize(){return sigma*sigcut;};
};

struct HydroPair {
    LJAttractCut inter;
    AtomID atom1, atom2;
    HydroPair(HydroAtom a1, HydroAtom a2) :
        inter(a1.getEpsilon(a2),
              (a1.sigma + a2.sigma) / 2,
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){return inter.energy(box.diff(atom1->x, atom2->x));};
    inline Vec forces(Box &box){return inter.forces(box.diff(atom1->x, atom2->x));};
};


////////////////////////////////////////////////////////////////////////
// Full LJish, with ε = ε₁₂ (indexed) and σ = (σ₁+σ₂)
// Potential is V(r) = ε (σ^n/r^n - 1)² - ε (r₀^-n - 1)²
// cutoff at some sigcut r₀

struct LJishAtom : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    flt repeps, sigma;
    flt exponent; // power
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two atoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]
    // these should be the same
    flt sigcut; // sigma units
    LJishAtom(){};
    LJishAtom(AtomID a, vector<flt> epsilons, flt repeps, flt sigma,
                            flt n, uint indx, flt cut) :
            AtomID(a), epsilons(epsilons), repeps(repeps), sigma(sigma),
            exponent(n), indx(indx), sigcut(cut){
        assert(indx < epsilons.size());
        //~ assert(sigma > 0.01);
    };
    LJishAtom(AtomID a, LJishAtom other) :
        AtomID(a), epsilons(other.epsilons), repeps(other.repeps),
            sigma(other.sigma), exponent(other.exponent), indx(other.indx),
        sigcut(other.sigcut){
            //~ assert(other.sigma > 0.01);
            //~ assert(sigma > 0.01);
        };
    flt getEpsilon(LJishAtom &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        return myeps;
    }
    flt getSigma(LJishAtom &other){
        return (sigma + other.sigma)/2.0;
    }
    flt maxsize(){return sigma*sigcut;};

};

struct LJishPair {
    flt epsilon, repeps, sigma, n, cutR, cutE;
    AtomID atom1, atom2;
    LJishPair(LJishAtom LJ1, LJishAtom LJ2) :
            epsilon(LJ1.getEpsilon(LJ2)),
            repeps(sqrt(LJ1.repeps * LJ2.repeps)),
            sigma(LJ1.getSigma(LJ2)),
            n((LJ1.exponent + LJ2.exponent) / 2),
            cutR(max(LJ1.sigcut,  LJ2.sigcut)),
            atom1(LJ1), atom2(LJ2){
        if(epsilon <= 0){
            cutR = 1;
            cutE = 0;
            epsilon = 0;
            return;
        }
        flt mid = (1-pow(cutR,-n));
        cutE = epsilon*(mid*mid);
    };
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt rsq = rij.squaredNorm()/(sigma*sigma);
        if(rsq > cutR * cutR){
            //~ printf("LJish: dx=%.2f, σ=%.2f, rsq=%.2f, cutR²=%.2f\n", rij.norm(), sigma, rsq, cutR);
            return 0;
        }

        flt mid = (1-pow(rsq,-n/2));
        //~ printf("LJish: rsq=%.2f, mid %.2f, epsilon %.2f, repeps %.2f\n", rsq, mid, epsilon, repeps);
        if (rsq > 1) return epsilon*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
    };
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        flt rsq = dsq/(sigma*sigma);
        if(rsq > cutR * cutR) return Vec::Zero();
        flt rmid = pow(rsq,-n/2); // σ^n/r^n
        flt fmagTimesR = 2*n*rmid*(rmid - 1); // 2 n r^-n(r^-n - 1)
        if (rsq < 1) return rij * (repeps * fmagTimesR / dsq);
        return rij * (epsilon * fmagTimesR / dsq); // r⃗ * 2 n * r^-n(r^-n - 1) / r² = 2 n r^-(n-1) * (r^-n - 1) r̂
    }
};

////////////////////////////////////////////////////////////////////////
// Full LJ, with ε = ε₁₂ (indexed) and σ = (σ₁ + σ₂)/2
// Potential is V(r) = ε (σ⁶/r⁶ - 1)² - ε (r₀⁻⁶ - 1)²
// cutoff at some sigcut r₀
// if ε₁₂ < 0, purely repulsive with ε = -ε₁₂


struct LJAttractRepulseAtom : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    flt sig;
    uint indx;
    flt sigcut; // sigma units
    LJAttractRepulseAtom(){};
    LJAttractRepulseAtom(AtomID a, vector<flt> epsilons, flt sigma, uint indx, flt cut) :
            AtomID(a), epsilons(epsilons), sig(sigma), indx(indx),
             sigcut(cut){
                 assert(indx < epsilons.size());
                 //~ printf("Made Atom  with ε=%.2f of size σ=%.2f\n", epsilons[indx], sig);
                 };
    LJAttractRepulseAtom(AtomID a, LJAttractRepulseAtom other) : AtomID(a), epsilons(other.epsilons),
        sig(other.sig), indx(other.indx), sigcut(other.sigcut){};
    flt getEpsilon(LJAttractRepulseAtom &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        return myeps;
    }
    flt maxsize(){return sig*sigcut;};
};

struct LJAttractRepulsePair {
    flt eps, sig;
    flt cutR, cutE;
    AtomID atom1, atom2;
    LJAttractRepulsePair(LJAttractRepulseAtom a1, LJAttractRepulseAtom a2) :
        eps(a1.getEpsilon(a2)), sig((a1.sig + a2.sig)/2.0),
        cutR(max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){
            if(eps <= 0){
                cutR = 1;
                cutE = 0;
                eps = abs(eps);
                return;
            }
            flt mid = (1-pow(cutR,-6));
            cutE = eps*(mid*mid);
        };
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt rsq = rij.squaredNorm()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                    //~ sqrt(rij.squaredNorm()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-pow(rsq,-3));
        //~ if(eps*(mid*mid) - cutE < 0) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                //~ rij.norm(), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        //~ }
        return eps*(mid*mid) - cutE;
    };
    inline Vec forces(Box &box){
        if(eps == 0) return Vec::Zero();
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec::Zero();
        flt rsix = pow(rsq,-3); // σ⁶/r⁶
        flt fmagTimesR = 12*eps*rsix*(rsix - 1);
        return rij * (fmagTimesR / dsq);
    };
};

////////////////////////////////////////////////////////////////////////
// Full LJ, with ε_A = ε₁₂ (indexed) and σ = (σ₁ + σ₂)/2
// Attractive Potential is V(r) = ε (σ⁶/r⁶ - 1)² - ε (r₀⁻⁶ - 1)²
// cutoff at some sigcut r₀
// For repulsive, ε_R = repeps

struct LJAttractFixedRepulseAtom : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    flt repeps, sig;
    uint indx;
    flt sigcut; // sigma units
    LJAttractFixedRepulseAtom(){};
    LJAttractFixedRepulseAtom(AtomID a, vector<flt> epsilons, flt repeps,
            flt sigma, uint indx, flt cut) :
            AtomID(a), epsilons(epsilons), repeps(repeps), sig(sigma),
            indx(indx), sigcut(cut){
                 assert(indx < epsilons.size());
                 };
    LJAttractFixedRepulseAtom(AtomID a, LJAttractFixedRepulseAtom other) :
        AtomID(a), epsilons(other.epsilons),
        repeps(other.repeps), sig(other.sig), indx(other.indx), sigcut(other.sigcut){};
    flt getEpsilon(LJAttractFixedRepulseAtom &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        return myeps;
    }
    flt maxsize(){return sig*sigcut;};
};

// r < 1        εr (1 - σ⁶/r⁶)² - εa (1 - σ⁶/R⁶)²
// 1 < r < R    εa (1 - σ⁶/r⁶)² - εa (1 - σ⁶/R⁶)²
// R < r        0

// so cutE = εa (1 - σ⁶/R⁶)²
struct LJAttractFixedRepulsePair {
    flt eps, repeps, sig;
    flt cutR, cutE;
    bool attract;
    AtomID atom1, atom2;
    LJAttractFixedRepulsePair(){};
    LJAttractFixedRepulsePair(LJAttractFixedRepulseAtom a1, LJAttractFixedRepulseAtom a2) :
        eps(abs(a1.getEpsilon(a2))), repeps(sqrt(a1.repeps * a2.repeps)),
        sig((a1.sig + a2.sig)/2.0),
        cutR(max(a1.sigcut, a2.sigcut)), attract(a1.getEpsilon(a2) > 0),
        atom1(a1), atom2(a2){
            if(!attract or eps == 0){
                cutR = 1;
                cutE = 0;
                eps = 0;
                attract = false;
                return;
            }
            flt mid = (1-pow(cutR,-6.0));
            cutE = eps*(mid*mid);
        };
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt rsq = rij.squaredNorm()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                    //~ sqrt(rij.squaredNorm()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-pow(rsq,-3)); // # 1 - σ⁶/r⁶
        //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                    //~ sqrt(rij.squaredNorm()), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        if (rsq > 1) return eps*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
        //~ flt E;
        //~ if (rsq > 1) E = eps*(mid*mid) - cutE;
        //~ else E = repeps*(mid*mid) - cutE;
        //~ if(E > 1e4)
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                    //~ sqrt(rij.squaredNorm()), E, eps, repeps, sig, cutR, cutE);
        //~ return E;
    };
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec::Zero();
        flt rsix = pow(rsq,-3); // σ⁶/r⁶
        flt fmagTimesR = 12*rsix*(rsix - 1);
        if (rsq < 1) return rij * (repeps * fmagTimesR / dsq);
        return rij * (eps * fmagTimesR / dsq);
        //~ flt fmag;
        //~ if (rsq < 1) fmag = repeps * fmagTimesR / dsq;
        //~ else fmag = eps * fmagTimesR / dsq;
        //~ if(fmag * rij.norm() > 1e4)
            //~ printf("Distance: %.2f Force: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n",
                    //~ sqrt(rij.squaredNorm()), fmag * rij.norm(), eps, repeps, sig, cutR, cutE);
        //~ return rij * fmag;
    };
};

////////////////////////////////////////////////////////////////////////
// Full LJ, with ε = √(ε₁ ε₂) and σ = (σ₁ + σ₂)/2
// Potential is V(r) = ε (σ⁶/r⁶ - 1)² - ε (r₀⁻⁶ - 1)²
// cutoff at some sigcut r₀
// if ε₁₂ < 0, purely repulsive with ε = -ε₁₂
struct LJDoubleAtom : public LJatom {
    flt epsrep;
    flt sigcut; // sigma units
    LJDoubleAtom(){};
    LJDoubleAtom(flt epsilon, flt epsrep, flt sigma, AtomID a, flt cut) :
            LJatom(epsilon, sigma, a), epsrep(epsrep), sigcut(cut){};
    LJDoubleAtom(AtomID a, LJDoubleAtom other) : LJatom(a, other),
        epsrep(other.epsrep), sigcut(other.sigcut){};
};

struct LJDoublePair : public LJAttractFixedRepulsePair {
    LJDoublePair(LJDoubleAtom a1, LJDoubleAtom a2) {
        eps = sqrt(a1.epsilon * a2.epsilon);
        repeps = sqrt(a1.epsrep * a2.epsrep);
        sig = (a1.sigma + a2.sigma)/2.0;
        cutR = max(a1.sigcut, a2.sigcut);
        attract = a1.epsilon * a2.epsilon > 0;
        atom1 = a1;
        atom2 = a2;
        if(!attract or eps == 0){
            cutR = 1;
            cutE = 0;
            eps = 0;
            attract = false;
            return;
        }
        flt mid = (1-pow(cutR,-6));
        cutE = eps*(mid*mid);
    };
};

struct EisMclachlanAtom : public AtomID {
    flt dist, sigmai; // dist is radius of Atom + radius of water (1.4 Å)
    EisMclachlanAtom(){};
    EisMclachlanAtom(flt dist, flt sigmai, AtomID a) : AtomID(a),
            dist(dist), sigmai(sigmai){};
    EisMclachlanAtom(AtomID a, EisMclachlanAtom other) : AtomID(a),
        dist(other.dist), sigmai(other.sigmai){};
    flt maxsize(){return dist;};
};

struct EisMclachlanPair {
    flt c0,c1,c2; // 1/R term, const term, R term in E (coefficients)
    // E = c₀/R + c₁ + c₂ R
    flt cutoff; // r1 + r2 (includes two water radii)
    AtomID atom1, atom2;
    EisMclachlanPair(EisMclachlanAtom a1, EisMclachlanAtom a2) :
        c0(-M_PI*(a1.sigmai*a2.dist - a2.sigmai*a1.dist)
                * (a1.dist*a1.dist - a2.dist*a2.dist)),
        c1(-2*M_PI*(a1.sigmai*a2.dist*a2.dist + a2.sigmai*a1.dist*a1.dist)),
        c2(M_PI*(a1.sigmai*a2.dist + a2.sigmai*a1.dist)),
        cutoff(a1.dist + a2.dist),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt R = rij.norm();
        if (R > cutoff) return 0;
        return c0/R + c1 + c2*R;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > (cutoff*cutoff)) return Vec::Zero();
        flt R = sqrt(dsq);
        return rij * ((c0/dsq-c2)/R);
    }
};

////////////////////////////////////////////////////////////////////////
//! Hertzian potential, with ε = √(ε₁ ε₂) and σ = (σ₁ + σ₂)/2
//! Potential is V(r) = ε/n (1 - r/σ)^n, with n = 5/2 usually
//! cutoff at r = σ

struct HertzianAtom : public AtomID {
    flt eps, sigma, exponent;
    HertzianAtom(){};
    HertzianAtom(AtomID a, flt eps, flt sigma, flt exponent=2.5) : AtomID(a),
            eps(eps), sigma(sigma), exponent(exponent){};
    HertzianAtom(AtomID a, HertzianAtom other) : AtomID(a),
        eps(other.eps), sigma(other.sigma), exponent(other.exponent){};
    flt maxsize(){return sigma;};
};

inline HertzianAtom hertzd(AtomID a, double eps, double sigma, double exponent=2.5){
    return HertzianAtom(a, eps, sigma, exponent);
};

struct EnergyForce {
    Vec f;
    flt E;
    EnergyForce(Vec f, flt E) : f(f), E(E){};
};

//----------------------------------------------------------------------
// Hertzian potential, as above, with ε = ε₁₂ and σ = σ₁₂ (both indexed)
// exponent is n = (n₁ + n₂)/2


struct HertzianAtomIndexed : public AtomID {
    vector<flt> epsilons; // for finding epsilons
    vector<flt> sigmas; // for finding epsilons
    flt exponent;
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two HertzianAtomIndexed atoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]; same for sigma
    HertzianAtomIndexed(){};
    HertzianAtomIndexed(AtomID a, vector<flt> epsilons, vector<flt> sigmas,
        uint indx, flt exponent=2.5) :
            AtomID(a), epsilons(epsilons), sigmas(sigmas), exponent(exponent),
            indx(indx){
                 assert(sigmas.size() == epsilons.size());
            };
    HertzianAtomIndexed(AtomID a, LJAtomIndexed other) : AtomID(a), epsilons(other.epsilons),
        sigmas(other.sigmas), indx(other.indx){};
    flt getEpsilon(HertzianAtomIndexed &other){
        assert(other.indx < epsilons.size());
        flt myeps = epsilons[other.indx];
        //~ assert(indx < other.epsilons.size());
        //~ assert(other.epsilons[indx] == myeps);
        return myeps;
    }
    flt getSigma(HertzianAtomIndexed &other){
        assert(other.indx < sigmas.size());
        flt myeps = sigmas[other.indx];
        //~ assert(indx < other.sigmas.size());
        //~ assert(other.sigmas[indx] == myeps);
        return myeps;
    }
    flt maxsize(){
        flt sigma = sigmas[0];
        for(uint i=1; i<sigmas.size(); ++i){
            if(sigma < sigmas[i]) sigma = sigmas[i];
        }
        return sigma;
    };
};

struct HertzianPair {
    flt eps, sig, exponent;
    AtomID atom1, atom2;
    HertzianPair(HertzianAtom a1, HertzianAtom a2) :
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        exponent((a1.exponent + a2.exponent)/2.0), atom1(a1), atom2(a2){};
    HertzianPair(HertzianAtomIndexed a1, HertzianAtomIndexed a2) :
        eps(a1.getEpsilon(a2)), sig(a1.getSigma(a2)),
        exponent((a1.exponent + a2.exponent)/2.0), atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig) return 0.0;
        flt R = sqrt(dsq);
        return eps * pow(1.0 - (R/sig), exponent) / exponent;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig) return Vec::Zero();
        flt R = sqrt(dsq);
        return rij * (eps * pow(1.0 - (R/sig), exponent-1) /sig/R);
    }
    inline EnergyForce EnergyForces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig) return EnergyForce(Vec::Zero(),0);
        flt R = sqrt(dsq);

        Vec f = rij * (eps * pow(1.0 - (R/sig), exponent-1) /sig/R);
        flt E = eps * pow(1.0 - (R/sig), exponent) / exponent;
        return EnergyForce(f, E);
    }
    //~ inline flt xrij(Box &box){
        //~ Vec rij = box.diff(atom1->x, atom2->x);
        //~ flt dsq = rij.squaredNorm();
        //~ if(dsq > sig*sig) return 0.0;
        //~ flt R = sqrt(dsq);
        //~ return (R*eps*(exponent-1)/sig/sig) * pow(1.0 - (R/sig), exponent-2);
    //~ }
    inline void fill(Box &box, ForcePairX &fpair){
        fpair.a1 = atom1;
        fpair.a2 = atom2;
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig){
            fpair.fij = Vec::Zero();
            fpair.xij = 0.0;
            return;
        }
        flt R = sqrt(dsq);
        fpair.fij = rij * (eps*pow(1.0 - (R/sig), exponent-1) /sig/R);
        flt Rs = R / sig;
        fpair.xij = -Rs*eps*pow(1.0 - (R/sig), exponent-2)*(1-exponent*Rs);
    }
};

////////////////////////////////////////////////////////////////////////
// Hertzian potential with drag, with ε = √(ε₁ ε₂) and σ = (σ₁ + σ₂)/2
// Potential is V(r) = ε/n (1 - r/σ)^n, with n = 5/2 usually
// cutoff at r = σ
// drag is f = -γv in the normal direction

struct HertzianDragAtom : public AtomID {
    flt eps, sigma, exponent, gamma;
    HertzianDragAtom(){};
    HertzianDragAtom(AtomID a, flt eps, flt sigma, flt gamma, flt exponent=2.5) : AtomID(a),
            eps(eps), sigma(sigma), exponent(exponent), gamma(gamma){};
    HertzianDragAtom(AtomID a, HertzianDragAtom other) : AtomID(a),
        eps(other.eps), sigma(other.sigma), exponent(other.exponent), gamma(other.gamma){};
    flt maxsize(){return sigma;};
};

struct HertzianDragPair {
    flt eps, sig, exponent, gamma;
    AtomID atom1, atom2;
    HertzianDragPair(HertzianDragAtom a1, HertzianDragAtom a2) :
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        exponent((a1.exponent + a2.exponent)/2.0),
        gamma((a1.gamma + a2.gamma)/2), atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig) return 0.0;
        flt R = sqrt(dsq);
        return eps * pow(1.0 - (R/sig), exponent) / exponent;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        Vec vij = atom1->v - atom2->v;
        flt dsq = rij.squaredNorm();
        if(dsq > sig*sig) return Vec::Zero();
        flt R = sqrt(dsq);
        Vec v_perp = rij * (vij.dot(rij)) / dsq;
        return rij * (eps * pow(1.0 - (R/sig), exponent-1) /sig/R) -
				(v_perp * gamma);
    }
};

/***********************************************************************
 * Lois-O'Hern potential (that's what I'm calling it) for approximating
 * sticky spheres. Force is linear repulsive, linear attractive, with
 * a repulsive epsilon, a sigma, and an attractive epsilon and
 * lengthscale.
 * See PRL 100, 028001 (2008)
 * They used C = 10^-2, l = 0. How do you have l=0? Well, that's a
 * discontinuity in the force, but not in the energy, so it works.
 *
 * More specifically:
F=\begin{cases}
Y\left(1-\frac{r}{\sigma}\right) & \frac{r}{\sigma}<1+C\\
\frac{CY}{\ell}\left(\frac{r}{\sigma}-\left(1+C+\ell\right)\right) & 1+C\leq\frac{r}{\sigma}\leq1+C+\ell
\end{cases}

U\left(r\right)=\begin{cases}
-\frac{Y\sigma}{2}\left(C\left(C+\ell\right)-\left(\frac{r}{\sigma}-1\right)^{2}\right) & \frac{r}{\sigma}<1+C\\
-\frac{CY\sigma}{2\ell}\left(C+\ell+1-\frac{r}{\sigma}\right)^{2} & 1+C\leq\frac{r}{\sigma}\leq1+C+\ell
\end{cases}

The minimum is at σ, its normal harmonic from 0 to (1+C)σ, and is then
"rounded" between (1+C)σ < r < (1+C+l)σ.
* With l = 0, its a harmonic potential, depth C²σ/2, minimum at σ, attractive
* out to (1+C)σ.

Parameters:
* C: how wide the "first half" of the well is
* l: how wide the "second half" of the well is
* l/C: how steep the "second half" of the well is; l == C means perfectly harmonic
* C(C+l)σ/2: depth of the well

*/

struct LoisOhernAtom : public AtomID {
    flt eps, sigma, C, l;
    LoisOhernAtom(){};
    LoisOhernAtom(AtomID a, flt eps, flt sigma, flt C, flt l) : AtomID(a),
            eps(eps), sigma(sigma), C(C), l(l){};
    LoisOhernAtom(AtomID a, LoisOhernAtom other) : AtomID(a),
        eps(other.eps), sigma(other.sigma), C(other.C), l(other.l){};
    flt maxsize(){return sigma*(1+C+l);};
};

struct LoisOhernPair {
    flt eps, sig, C, l, sigcut;
    AtomID atom1, atom2;
    LoisOhernPair(LoisOhernAtom a1, LoisOhernAtom a2) :
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        C((a1.C + a2.C)/2.0), l((a1.l + a2.l)/2.0), sigcut(sig*(1+C+l)),
        atom1(a1), atom2(a2){};
    LoisOhernPair(LoisOhernAtom a1, LoisOhernAtom a2, flt eps, flt sig, flt C, flt l) :
        eps(eps), sig(sig), C(C), l(l), sigcut(sig*(1+C+l)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
         // using >= to prevent NaNs when l = 0
        if(dsq >= sigcut*sigcut) return 0.0;
        flt R = sqrt(dsq)/sig;
         // using <= to prevent NaNs when l = 0
        if(R <= 1 + C){
            flt dR = R-1;
            return -eps*sig/2*(C*(C+l) - dR*dR);
        }

        flt dR2 = C+l+1 - R;
        return -C*eps*sig/2/l*dR2*dR2;
    }

    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq >= sigcut*sigcut) return Vec::Zero();
        flt R = sqrt(dsq);
        flt rsig = R/sig;

        if(rsig <= 1 + C){
            flt dR = rsig-1;
            return rij*(-eps*dR/R);
        }

        flt dR2 = rsig-(C+l+1);
        return rij*(C*eps/l*dR2/R);
    }
};



struct LoisOhernPairMinCLs : public LoisOhernPair {
    LoisOhernPairMinCLs(LoisOhernAtom a1, LoisOhernAtom a2) :
        LoisOhernPair(a1, a2, sqrt(a1.eps * a2.eps),
        (a1.sigma + a2.sigma)/2.0,
        (a1.C < a2.C ? a1.C : a2.C),
        (a1.l < a2.l ? a1.l : a2.l)){};
};

////////////////////////////////////////////////////////////////////////
/***********************************************************************
 * LoisLin atoms
 * Harmonic repulsive, fixed-force attractive
 * epsilon is the strength of the harmonic Interaction
 * f is the fixed-force amount (real units)
 * l is the width of the fixed-force region (real units)
 * To combine, we use a geometric average of epsilon, and everything else
 * is averaged.
 */

struct LoisLinAtom : public AtomID {
    flt eps, sigma, f, l;
    LoisLinAtom(){};
    LoisLinAtom(AtomID a, flt eps, flt sigma, flt depth, flt width) : AtomID(a),
            eps(eps), sigma(sigma), f(width > 0 ? depth/width : 0), l(width){};
    LoisLinAtom(AtomID a, LoisLinAtom other) : AtomID(a),
        eps(other.eps), sigma(other.sigma), f(other.f), l(other.l){};
    flt maxsize(){return sigma*(1+l);};
};

struct LoisLinPair {
    flt eps, sig, f, l, sigcut;
    AtomID atom1, atom2;
    LoisLinPair(LoisLinAtom a1, LoisLinAtom a2) :
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        f((a1.f + a2.f)/2.0), l((a1.l + a2.l)/2.0), sigcut(sig +l),
        atom1(a1), atom2(a2){};
    LoisLinPair(LoisLinAtom a1, LoisLinAtom a2, flt eps, flt sig, flt f, flt l) :
        eps(eps), sig(sig), f(f), l(l), sigcut(sig+l),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
         // using >= to prevent NaNs when l = 0
        if(dsq >= sigcut*sigcut) return 0.0;
        flt R = sqrt(dsq);
         // using <= to prevent NaNs when l = 0
        if(R <= sig){
			flt dR = 1.0 - (R/sig);
            return eps/2 * dR*dR - f*l;
        }

        return -f*(sig + l - R);
    }

    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.squaredNorm();
        if(dsq >= sigcut*sigcut) return Vec::Zero();
        flt R = sqrt(dsq);

        if(R <= sig){
			flt dR = 1.0 - (R/sig);
			return rij * (eps * dR /sig/R);
        }

        return -rij*(f/R);
    }
};

struct LoisLinPairMin : public LoisLinPair {
    LoisLinPairMin(LoisLinAtom a1, LoisLinAtom a2) :
        LoisLinPair(a1, a2, sqrt(a1.eps * a2.eps),
        (a1.sigma + a2.sigma)/2.0,
        (a1.f < a2.f ? a1.f : a2.f),
        (a1.l < a2.l ? a1.l : a2.l)){};
};

////////////////////////////////////////////////////////////////////////
class LJsimple : public Interaction {
    protected:
        vector<LJatom> atoms;
        PairList ignorepairs;
        AtomID get_id(Atom* a);

    public:
        LJsimple(flt cutoffdist, vector<LJatom> atms=vector<LJatom>());
         // cutoffdist in sigma units

        inline void add(LJatom a){
            assert(a.n() <= atoms.size());
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            atoms[a.n()] = a;};
        inline void add(AtomID a, flt epsilon, flt sigma){
            add(LJatom(epsilon,sigma,a));};
        inline void ignore(AtomID a, AtomID b){ignorepairs.add_pair(a,b);};
        inline void ignore(Atom* a, Atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint atoms_size() const{return (uint) atoms.size();};
        flt energy(Box &box);
        flt pressure(Box &box);
        void setForces(Box &box);
        //~ ~LJsimple(){};
};

template <class A, class P>
class SCBoxed : public Interaction {
    protected:
        sptr<SCBox> box;
        sptr<AtomVec> atoms;
        vector<A> group; // data on which atoms interact with the wall
    public:
        SCBoxed(sptr<AtomVec> atomv, sptr<SCBox> box)
            : box(box), atoms(atomv){};
        inline void add(A atm){group.push_back(atm);};
        flt energy(Box &box);
        flt pressure(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        inline vector<A> &atom_list(){return group;};
};

template <class A, class P>
class SimpleListed : public Interaction {
    protected:
        vector<A> atoms;

    public:
        SimpleListed(){};
        inline void add(A atm){atoms.push_back(atm);};
        //flt energy(Box &box, IDPair &pair);
        flt energy(Box &box);
        flt pressure(Box &box);
        uint size(){return ((uint) (atoms.size()));};
        //inline flt energy_pair(P pair, Box &box){return pair.energy(box);}; // This may need to be written!
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        //inline Vec forces_pair(P pair, Box &box){return pair.forces(box);}; // This may need to be written!
        inline vector<A> &atom_list(){return atoms;};
        ~SimpleListed(){};
};

template <class A, class P>
class NListed : public Interaction {
    // NeighborList keeps track of Atom* pairs within a particular distance.
    // NListed keeps track of atoms with additional properties
    // that interact through a particular Interaction.
    // class A needs to inherit from AtomID, and also have an initialization
    // method A(AtomID a, A other)
    // you also need to implement a couple methods below

    // Implementation detail:
    // note that NeighborList maintains its own MetaGroup, so that when
    // a member of A is created and passed to this group, the AtomID in that
    // member has an A.n() value referring to its place in the *original*
    // AtomGroup, not this one. This is why we need an A(AtomID, A) method;
    // so we can make a new A with all the same properties as before,
    // but with the n() referring to the NeighborList maintained AtomGroup.
    protected:
        vector<A> atoms; // NOT all full.
                         // This is the length of the AtomVec, and if
                         // Atom i has been added, then atoms[i] is
                         // filled in this vector
        vector<P> pairs;
        sptr<NeighborList> neighbors;
        uint lastupdate;
    public:
        NListed(sptr<AtomVec> vec, sptr<NeighborList> neighbors) :
            atoms(vec->size()), neighbors(neighbors), lastupdate(0){}; //group(vec),
        inline void add(A atm){
            neighbors->add(atm, atm.maxsize());
            // assert(atoms.size() > atm.n());
            atoms[atm.n()] = atm;
        }
        NListed(sptr<Box> box, sptr<AtomVec> atomv, const flt skin) :
            atoms(atomv->size()),
            neighbors(new NeighborList(box, atomv, skin)), lastupdate(0){};
        void update_pairs();
        P getpair(IDPair &pair){
            return P(atoms[pair.first().n()], atoms[pair.last().n()]);}
        A& getatom(uint n){return atoms[n];}
        flt energy(Box &box, IDPair &pair);
        flt energy(Box &box);
        
        /// number of atom pairs with E != 0
        unsigned long long contacts(Box &box);
        /// number of atom pairs with E > 0 
        unsigned long long overlaps(Box &box);
        flt pressure(Box &box);
        inline vector<P> &pairiter(){return pairs;};
        uint size(){return ((uint) (atoms.size()));};
        inline flt energy_pair(P pair, Box &box){return pair.energy(box);}; // This may need to be written!
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        inline Vec forces_pair(P pair, Box &box){return pair.forces(box);}; // This may need to be written!
        inline vector<A> &atom_list(){return atoms;};
        inline sptr<NeighborList> nlist(){return neighbors;};
        //~ flt energy_test(flt dist);
        //~ flt force_test(flt dist);
        ~NListed(){};
};

template <class A, class P>
class NListedVirial : public InteractionPairsX {
    private:
        NListed<A,P> nlisted;
    public:
        NListedVirial(sptr<AtomVec> vec, sptr<NeighborList> neighbors) :
                nlisted(vec, neighbors){};
        void setForces(Box &box){nlisted.setForces(box);};
        void setForces(Box &box, FPairXFunct*);
        virtual inline flt setForcesGetPressure(Box &box){return nlisted.setForcesGetPressure(box);};
        virtual flt setForcesGetEnergy(Box &box);
        virtual inline flt energy(Box &box){return nlisted.energy(box);};
        virtual inline flt pressure(Box &box){return nlisted.pressure(box);};
        inline void add(A atm){nlisted.add(atm);}
        inline sptr<NeighborList> nlist(){return nlisted.nlist();};
};

template <class A, class P>
flt SCBoxed<A, P>::energy(Box &newbox){
    flt E = 0;
    Atom wallatom;
    wallatom.m = NAN;
    AtomID wallid = AtomID(&wallatom, atoms->size());
    typename vector<A>::iterator it;
    for(it = group.begin(); it != group.end(); ++it){
        Vec r0 = (*it)->x;
        Vec edger = box->edgedist(r0);
        wallatom.x = r0 + (edger*2);
        E += P(*it, HertzianAtom(wallid, *it)).energy(newbox) / 2; // no energy from fake Atom
    }
    return E;
};

template <class A, class P>
void SCBoxed<A, P>::setForces(Box &newbox){
    Atom wallatom;
    wallatom.m = NAN;
    AtomID wallid = AtomID(&wallatom, atoms->size());
    typename vector<A>::iterator it;
    for(it = group.begin(); it != group.end(); ++it){
        Vec r0 = (*it)->x;
        Vec edger = box->edgedist(r0);
        wallatom.x = r0 + (edger*2);
        P pair = P(*it, HertzianAtom(wallid, *it));
        Vec f = pair.forces(*box);
        (*it)->f += f;
    }
};

template <class A, class P>
flt SCBoxed<A, P>::setForcesGetPressure(Box &newbox){
    flt p = 0;
    Atom wallatom;
    wallatom.m = NAN;
    AtomID wallid = AtomID(&wallatom, atoms->size());
    typename vector<A>::iterator it;
    for(it = group.begin(); it != group.end(); ++it){
        Vec r0 = (*it)->x;
        Vec dr = box->edgedist(r0) * 2;
        wallatom.x = r0 + dr;
        P pair = P(*it, HertzianAtom(wallid, *it));
        Vec f = pair.forces(*box);
        (*it)->f += f;
        p += f.dot(dr);
    }
    return p;
};

template <class A, class P>
flt SCBoxed<A, P>::pressure(Box &newbox){
    flt p = 0;
    Atom wallatom;
    wallatom.m = NAN;
    AtomID wallid = AtomID(&wallatom, atoms->size());
    typename vector<A>::iterator it;
    for(it = group.begin(); it != group.end(); ++it){
        Vec r0 = (*it)->x;
        Vec dr = box->edgedist(r0) * 2;
        wallatom.x = r0 + dr;
        P pair = P(*it, HertzianAtom(wallid, *it));
        Vec f = pair.forces(*box);
        p += f.dot(dr);
    }
    return p;
};

template <class A, class P>
flt SimpleListed<A, P>::energy(Box &box){
    flt E = 0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); ++it1){
        for(it2 = it1+1; it2 != atoms.end(); ++it2){
            E += P(*it1, *it2).energy(box);
        }
    }
    return E;
};

template <class A, class P>
void SimpleListed<A, P>::setForces(Box &box){
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); ++it1){
        for(it2 = it1+1; it2 != atoms.end(); ++it2){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            (*it1)->f += f;
            (*it2)->f -= f;
        }
    }
};

template <class A, class P>
flt SimpleListed<A, P>::setForcesGetPressure(Box &box){
    flt p=0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); ++it1){
        for(it2 = it1+1; it2 != atoms.end(); ++it2){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            (*it1)->f += f;
            (*it2)->f -= f;
            Vec r = box.diff((*it1)->x, (*it2)->x);
            p += r.dot(f);
        }
    }
    //~ cout << "Set forces, got pressure" << p << '\n';
    return p;
};

template <class A, class P>
flt SimpleListed<A, P>::pressure(Box &box){
    flt p=0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); ++it1){
        for(it2 = it1+1; it2 != atoms.end(); ++it2){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            Vec r = box.diff((*it1)->x, (*it2)->x);
            p += r.dot(f);
        }
    }
    return p;
};

template <class A, class P>
void NListed<A, P>::update_pairs(){
    if(lastupdate == neighbors->which()) return; // already updated

    lastupdate = neighbors->which();
    pairs.clear();
    vector<IDPair>::iterator pairit;
    for(pairit = neighbors->begin(); pairit != neighbors->end(); ++pairit){
        // assert(atoms.size() > pairit->first().n());
        // assert(atoms.size() > pairit->last().n());
        A firstatom = atoms[pairit->first().n()];
        A secondatom = atoms[pairit->last().n()];
        pairs.push_back(P(firstatom, secondatom));
    }
}

template <class A, class P>
flt NListed<A, P>::energy(Box &box, IDPair &pair){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    P Epair = P(atoms[pair.first().n()],atoms[pair.last().n()]);
    //~ Vec dist = box.diff(Epair.atom1->x, Epair.atom2->x);
    return energy_pair(Epair, box);
};

template <class A, class P>
unsigned long long NListed<A, P>::contacts(Box &box){
    unsigned long long Nc = 0;
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        flt E = energy_pair(*it, box);
        if(E != 0.0){Nc++;};
    }
    return Nc;
};

template <class A, class P>
unsigned long long NListed<A, P>::overlaps(Box &box){
    unsigned long long Nc = 0;
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        flt E = energy_pair(*it, box);
        if(E > 0.0){Nc++;};
    }
    return Nc;
};

template <class A, class P>
flt NListed<A, P>::energy(Box &box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt E = 0;
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        //~ Vec dist = box.diff(it->atom1->x, it->atom2->x);
        E += energy_pair(*it, box);
    }
    return E;
};

template <class A, class P>
void NListed<A, P>::setForces(Box &box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        Vec f = forces_pair(*it, box);
        it->atom1->f += f;
        it->atom2->f -= f;
        //~ assert(f.squaredNorm() < 1000000);
    }
};

template <class A, class P>
flt NListedVirial<A, P>::setForcesGetEnergy(Box &box){
    nlisted.update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt E = 0;
    vector<P> & pairs = nlisted.pairiter();
    for(typename vector<P>::iterator it = pairs.begin(); it != pairs.end(); ++it){
        EnergyForce EF = it->EnergyForces(box);
        it->atom1->f += EF.f;
        it->atom2->f -= EF.f;
        E += EF.E;
    }
    return E;
}

//~ template <class A, class P>
//~ flt NListed<A, P>::energy_test(flt dist){
    //~ typename vector<P>::iterator it = pairs.begin();
    //~ Atom a1 = Atom(it->first());
    //~ Atom a2 = Atom(it->second());
    //~ A atm1 = A(it->first(), &a1);
    //~ A atm2 = A(it->second(), &a1);
    //~ a1.x = Vec::Zero();
    //~ a1.v = Vec::Zero();
    //~ a1.a = Vec::Zero();
    //~ a1.f = Vec::Zero();
    //~
    //~ a2.x = Vec::Zero();
    //~ a2.v = Vec::Zero();
    //~ a2.a = Vec::Zero();
    //~ a2.f = Vec::Zero();
    //~ a2.x.setx(dist);
    //~ P Epair = P(atm1,atm2);
    //~ return energy_pair(Epair, infbox);
//~ };

template <class A, class P>
void NListedVirial<A, P>::setForces(Box &box, FPairXFunct* funct){
    nlisted.update_pairs(); // make sure the LJpairs match the neighbor list ones
    ForcePairX myfpair;
    typename vector<P>::iterator it;
    for(it = nlisted.pairiter().begin(); it != nlisted.pairiter().end(); ++it){
        it->fill(box, myfpair);
        //~ assert(!isnan(myfpair.fij[0]));
        funct->run(&myfpair);
        it->atom1->f += myfpair.fij;
        it->atom2->f -= myfpair.fij;
        //~ assert(f.squaredNorm() < 1000000);
    }
};

template <class A, class P>
flt NListed<A, P>::setForcesGetPressure(Box &box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt p=0;
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        Vec f = forces_pair(*it, box);
        it->atom1->f += f;
        it->atom2->f -= f;
        Vec r = box.diff(it->atom1->x, it->atom2->x);
        p += r.dot(f);
    }
    //~ cout << "Set forces, got pressure" << p << '\n';
    return p;
};

template <class A, class P>
flt NListed<A, P>::pressure(Box &box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt p=0;
    //~ printf("Updating pressure...\n");
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); ++it){
        Vec f = forces_pair(*it, box);
        Vec r = box.diff(it->atom1->x, it->atom2->x);
        //Vec r = it->atom1->x - it->atom2->x;
        p += r.dot(f);
        //~ assert(!isnan(p));
    }
    return p;
};

class Charges : public Interaction {
    protected:
        vector<Charged> atoms;
        PairList ignorepairs;
        flt screen;
        flt k;
        AtomID get_id(Atom* a);

    public:
        Charges(flt screenlength, flt k=1, vector<Charged> atms=vector<Charged>());
         // cutoffdist in sigma units

        inline void add(Charged a){atoms.push_back(a); return;};
        inline void add(AtomID a, flt q){add(Charged(q,a));};
        inline void ignore(AtomID a, AtomID b){ignorepairs.add_pair(a,b);};
        inline void ignore(Atom* a, Atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint size() const{return (uint) atoms.size();};
        flt energy(Box &box);
        flt pressure(Box &box);
        void setForces(Box &box);
        //~ ~LJsimple(){};
};

////////////////////////////////////////////////////////////////////////
// A wall

struct WallAtom : AtomRef {
    public:
        flt sigma;
        flt epsilon;
    public:
        WallAtom(AtomID a, flt sigma, flt epsilon=1.0) :
            AtomRef(a), sigma(sigma), epsilon(epsilon){};
};

class SoftWall : public Interaction {
    protected:
        Vec loc;
        Vec norm;
        flt expt;
        vector<WallAtom> group;
        flt lastf;
    public:
        SoftWall(Vec loc, Vec norm, flt expt=2.0) :
            loc(loc), norm(norm.normalized()), expt(expt), lastf(NAN){};
        void add(WallAtom a){group.push_back(a);};
        flt energy(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        flt pressure(Box &box);

        void setLoc(Vec newloc){loc = newloc;};
        Vec getLoc(){return loc;};
        void setNorm(Vec newNorm){norm = newNorm.normalized();};
        Vec getNorm(){return norm;};

        flt get_last_f(){return lastf;};
};

class SoftWallCylinder : public Interaction {
    protected:
        Vec loc;
        Vec axis;
        flt radius;
        flt expt;
        flt lastf;
        vector<WallAtom> group;
    public:
        SoftWallCylinder(Vec loc, Vec axis, flt radius, flt expt=2.0) :
            loc(loc), axis(axis.normalized()), radius(radius), expt(expt), lastf(NAN){};
        void add(WallAtom a){
            if(a.sigma > radius*2)
                throw std::invalid_argument("SoftWallCylinder::add: sigma must be less than cylinder diameter");
            group.push_back(a);
        };
        flt energy(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        flt pressure(Box &box);

        void setLoc(Vec new_loc){loc = new_loc;};
        Vec getLoc(){return loc;};
        void setAxis(Vec new_axis){axis = new_axis.normalized();};
        Vec getAxis(){return axis;};
        flt get_last_f(){return lastf;};
};

#ifdef VEC2D
class WalledBox2D : public OriginBox {
    protected:
        bool xwalls, ywalls;
        vector<SoftWall*> walls;
    public:
        WalledBox2D(Vec size, bool xwalled, bool ywalled, flt expt=2.0) :
            OriginBox(size), xwalls(xwalled), ywalls(ywalled){
                if(xwalls){
                    walls.push_back(new SoftWall(
                        Vec(-size[0]/2.0, 0), Vec(1, 0), expt));
                    walls.push_back(new SoftWall(
                        Vec(size[0]/2.0, 0), Vec(-1, 0), expt));
                }
                if(ywalls){
                    walls.push_back(new SoftWall(
                        Vec(0, -size[0]/2.0), Vec(0, 1), expt));
                    walls.push_back(new SoftWall(
                        Vec(0, size[0]/2.0), Vec(0, -1), expt));
                }
            };
        Vec diff(Vec r1, Vec r2){
            Vec dr = r1 - r2;
            if(!xwalls) dr[0] = remainder(r1[0], boxsize[0]);
            if(!ywalls) dr[1] = remainder(r1[1], boxsize[1]);
            return dr;
        }
        Vec diff(Vec r1, Vec r2, array<int,NDIM> boxes){
            Vec dr = r1 - r2;
            if(!xwalls) dr[0] -= boxes[0]*boxsize[0];
            if(!ywalls) dr[1] -= boxes[1]*boxsize[1];
            return dr;
        }
        flt V(){return boxsize[0] * boxsize[1];};
        flt L(){return (boxsize[0] + boxsize[1])/2.0;};

        flt resize(flt factor){
            boxsize *= factor;
            vector<SoftWall*>::iterator it;
            for(it = walls.begin(); it != walls.end(); ++it){
                (*it)->setLoc((*it)->getLoc() * factor);
            };
            return V();
        }
        flt resizeV(flt newV){
            flt curV = V();
            return resize(pow(newV/curV, OVERNDIM));
        }
        Vec randLoc(flt walldist){
            Vec bxvec = boxsize;
            if(xwalls) bxvec[0] -= walldist/2.0;
            if(ywalls) bxvec[1] -= walldist/2.0;
            Vec v = randVecBoxed();
            for(uint i=0; i<NDIM; ++i){
                v[i] *= bxvec[i];
            }
            return diff(v, Vec::Zero());
        };
        vector<SoftWall*> getWalls() {return walls;};
        ~WalledBox2D(){
            for(vector<SoftWall*>::iterator it = walls.begin(); it != walls.end(); ++it){
                delete *it;
            }
        };
};
#endif

inline flt confineRange(flt minimum, flt val, flt maximum){
    if(val <= minimum) return minimum;
    if(val >= maximum) return maximum;
    return val;
}

class SCAtomVec : public virtual AtomGroup {
    // this is an AtomGroup which actually owns the atoms, which are
    // arranged in pairs.
    private:
        AtomVec atoms;
    public:
        SCAtomVec(vector<double> masses) : atoms((uint) masses.size()*2, 0.0){
            for(uint i=0; i < masses.size(); ++i){
                atoms[2*i].m = masses[i]/2;
                atoms[2*i+1].m = masses[i]/2;
            }
        };
        SCAtomVec(uint N, flt mass) : atoms(N*2, mass/2.0){};
        inline AtomVec &vec(){return atoms;};
        inline Atom& operator[](cuint n){return atoms[n];};
        inline Atom& operator[](cuint n) const {return atoms[n];};
        inline AtomID get_id(cuint n){return atoms.get_id(n);};
        inline IDPair pair(cuint n){return IDPair(atoms.get_id(n*2), atoms.get_id(n*2 + 1));};
        inline uint size() const {return atoms.size();};
        inline uint pairs() const {return atoms.size()/2;};
        ~SCAtomVec(){};
};

struct SpheroCylinderDiff{
    Vec delta, r;
    flt lambda1, lambda2;
};

struct SCPair {
    IDPair &p1;
    IDPair &p2;
    flt l1, l2;
    SCPair(IDPair &p1, IDPair &p2, flt l1, flt l2) :
        p1(p1), p2(p2), l1(l1), l2(l2){};
    SCPair(IDPair &p1, IDPair &p2, flt l) :
        p1(p1), p2(p2), l1(l), l2(l){};
    SCPair(const SCPair &other) : p1(other.p1), p2(other.p2),
            l1(other.l1), l2(other.l2){}
    SpheroCylinderDiff NearestLoc(Box &box);
    void applyForce(Box &box, Vec f, SpheroCylinderDiff diff, flt I1, flt I2);
    inline void applyForce(Box &box, Vec f, SpheroCylinderDiff diff, flt I){
        return applyForce(box, f, diff, I, I);
    }
};

struct SCSpringPair : public SCPair {
    /// Harmonic repulsive interactions between spherocylinders.
    flt eps, sig;

    SCSpringPair(IDPair &p1, IDPair &p2, flt eps, flt sig, flt l1, flt l2) :
        SCPair(p1, p2, l1, l2), eps(eps), sig(sig){};
    SCSpringPair(IDPair &p1, IDPair &p2, flt eps, flt sig, flt l) :
        SCPair(p1, p2, l), eps(eps), sig(sig){};

    inline flt maxdist(){return sig + (l1+l2)/2;};
    inline flt maxdelta(){return sig;};

    flt energy(Box &box, SpheroCylinderDiff diff){
        flt dsq = diff.delta.squaredNorm();
        if(dsq > sig*sig) return 0;
        flt d = sqrt(dsq);
        flt dsig = d-sig;
        return dsig*dsig*eps/2;
    }
    Vec forces(Box &box, SpheroCylinderDiff diff){
        flt dsq = diff.delta.squaredNorm();
        if(dsq > sig*sig) return Vec::Zero();
        flt dmag = sqrt(dsq);
        Vec dhat = diff.delta / dmag;

        return dhat * (eps * (sig - dmag));
    };
};

class SCSpringList : public Interaction {
    /// Harmonic repulsive interactions between spherocylinders.
    private:
        SCAtomVec *scs;
        flt eps, sig;
        vector<flt> ls;
        set<array<uint, 2> > ignore_list;
    public:
        SCSpringList(SCAtomVec *scs, flt eps, flt sig, flt l) :
            scs(scs), eps(eps), sig(sig), ls(scs->pairs(), l){};
        SCSpringList(SCAtomVec *scs, flt eps, flt sig, vector<flt> ls) :
            scs(scs), eps(eps), sig(sig), ls(ls){};
        flt energy(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box){setForces(box); return NAN;};
        flt pressure(Box &box){return NAN;};
        void ignore(uint n1, uint n2){
            if(n1 > n2){uint n3=n1; n1=n2; n2=n3;}
            array<uint, 2> pair;
            pair[0] = n1;
            pair[1] = n2;
            ignore_list.insert(pair);
        }
        void ignore(AtomID a1, AtomID a2){ignore(a1.n() / 2, a2.n() / 2);}

        ~SCSpringList(){};
};

#endif
