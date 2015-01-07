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
 * interaction<N> - blindly takes N vectors, 
 * and calculates energy or forces.
 * 
 * atomgroup(N) (or (N,M) with array?)- 
 * has a list of atoms, can tell mass, vel, loc.
 * keeps forces?
 * possibly has M groups of groups?
 * iterator over atoms necessary
 * perhaps class iterator over all of them?
 * perhaps class iterator over all atoms?
 * has "uint size()" function, says number of atoms
 * 
 * interactgroup - has a single interaction<N>, and pointers to however
 * many atomgroups it needs, and returns energy and/or sets forces.
 * one example: neighbors; yields pairs of i,i+1.
 * second example: neighbor list, keeps track of nearby atoms.
 * perhaps class iterator over all of them?
 * 
 * perhaps split: interactiongroup has a group, and an interaction, and 
 * applies interaction to each pair in each group.
 * Not incompatible.
 * 
 * integrator - keeps a list of interaction groups, can get energy.
 * Also goes through list to integrate over time.
 * 
 * statistic - keeps a list of interaction groups, derives a statistic 
 * from them.
 * 
 * statkeeper - keeps a list of statistics, can write to file?
 * 
 * constants - keeps track of constants; in map, or as properties?
 * perhaps a struct.
*/

//typedef double flt; defined in vecrand.hpp
//typedef unsigned int uint;
class interaction {
    public:
        virtual flt energy(Box &box)=0;
        virtual void setForces(Box &box)=0;
        virtual flt setForcesGetPressure(Box &box){return NAN;};
        virtual flt pressure(Box &box)=0; // sum_{<i,j>} of r_{ij} \cdot F_{ij}, or equivalently // sum of r_{i} \cdot F_{i}
        virtual ~interaction(){};
};

/***********************************************************************
 * Interaction Basics
 */
class interactpair {
    public:
        virtual flt energy(const Vec diff)=0;
        virtual Vec forces(const Vec diff)=0;
        virtual ~interactpair(){};
};

static const flt LJr0 = pow(2.0, 1.0/6.0);
static const flt LJr0sq = pow(2.0, 1.0/3.0);

class LJrepulsive {
    protected:
        flt epsilon;
        flt sigma;
    public:
        LJrepulsive(const flt epsilon, const flt sigma):
            epsilon(epsilon), sigma(sigma){};
        inline static flt energy(const Vec diff, const flt eps, const flt sig){
            flt rsq = diff.sq()/(sig*sig);
            if(rsq > 1) return 0;
            flt rsix = rsq*rsq*rsq;
            //~ return eps*(4*(1/(rsix*rsix) - 1/rsix) + 1);
            flt mid = (1-1/rsix);
            return eps*(mid*mid);
        };
        inline flt energy(const Vec& diff){return energy(diff, epsilon, sigma);};
        inline static Vec forces(const Vec diff, const flt eps, const flt sig){
            flt dsq = diff.sq();
            flt rsq = dsq/(sig*sig);
            if(rsq > 1) return Vec();
            flt rsix = rsq*rsq*rsq; //r^6 / sigma^6
            //~ flt fmagTimesR = eps*(4*(12/(rsix*rsix) - 6/rsix));
            flt fmagTimesR = 12*eps/rsix*(1/rsix - 1);
            //~ cout << "Repulsing " << diff << "with force"
                 //~ << diff * (fmagTimesR / dsq) << '\n';
            //~ cout << "mag: " << diff.mag() << " sig:" << sig << " eps:" << eps
                 //~ << "F: " << (diff * (fmagTimesR / dsq)).mag() << '\n';
            //~ cout << "E: " << energy(diff, sig, eps) << " LJr0: " << LJr0 << ',' << LJr0sq << "\n";
            return diff * (fmagTimesR / dsq);
        };
        inline Vec forces(const Vec& diff){return forces(diff, epsilon, sigma);};
};

class LJattract {
    protected:
        flt epsilon;
        flt sigma;
    public:
        LJattract(const flt epsilon, const flt sigma):
            epsilon(epsilon), sigma(sigma){};
        inline static flt energy(const Vec diff, const flt eps, const flt sig){
            flt rsq = diff.sq()/(sig*sig);
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
            flt dsq = diff.sq();
            flt rsq = dsq/(sig*sig);
            if(rsq < 1) return Vec();
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

class LJattractCut {
    // Purely attractive
    protected:
        flt epsilon;
        flt sigma;
        flt cutR, cutE;
    public:
        inline LJattractCut(const flt epsilon, const flt sigma, const flt cutsig):
            epsilon(epsilon), sigma(sigma), cutR(cutsig), 
            cutE(LJattract::energy(cutR) * epsilon){};
        inline static flt energy(const Vec diff, const flt eps,
                                    const flt sig, const flt cutsig){
            if(eps == 0) return 0;
            if(diff.sq() > (cutsig*cutsig*sig*sig)) return 0;
            return (LJattract::energy(diff, eps, sig) - 
                eps*LJattract::energy(cutsig));
        };
        inline flt energy(const Vec& diff){
            if(epsilon == 0) return 0;
            if(diff.sq() > (cutR*cutR*sigma*sigma)) return 0;
            return LJattract::energy(diff, epsilon, sigma) - cutE;
        };
        inline static Vec forces(const Vec diff, const flt eps, 
                                    const flt sig, const flt cutsig){
            if(eps == 0) return Vec();
            flt dsq = diff.sq();
            flt rsq = dsq/(sig*sig);
            if(rsq < 1 or rsq > (cutsig*cutsig)) return Vec();
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
            flt rsq = diff.sq()/(sigma*sigma);
            if(rsq > (cutR*cutR)) return 0;
            flt rsix = rsq*rsq*rsq;
            flt mid = (1-1/rsix);
            return epsilon*(mid*mid-1) - cutE;
        };
        inline static Vec forces(const Vec diff, const flt eps, const flt sig, const flt cutsig){
            flt dsq = diff.sq();
            flt rsq = dsq/(sig*sig);
            if(rsq > (cutsig*cutsig)) return Vec();
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

class spring : public interactpair {
    protected:
        flt springk;
        flt x0;
    public:
        spring(const flt k, const flt x0) : springk(k),x0(x0){};
        flt energy(const Vec diff);
        Vec forces(const Vec diff);
        ~spring(){};
};

class bondangle {
    // two vectors for E and F are from the central one to the other 2
    protected:
        flt springk;
        flt theta0;
        bool usecos;
        
    public:
        bondangle(const flt k, const flt theta, const bool cosine=false)
                        :springk(k), theta0(theta), usecos(cosine){};
        inline static flt get_angle(const Vec& r1, const Vec& r2){
            flt costheta = r1.dot(r2) / r1.mag() / r2.mag();
            if(costheta > 1) costheta = 1;
            else if(costheta < -1) costheta = -1;
            return acos(costheta);
        };
        flt energy(const Vec& diff1, const Vec& diff2);
        Nvector<Vec,3> forces(const Vec& diff1, const Vec& diff2);
        ~bondangle(){};
};

#ifdef VEC3D
class dihedral {
    // vectors are from 1 to 2, 2 to 3, 3 to 4
    protected:
        vector<flt> coscoeffs;
        vector<flt> sincoeffs;
        bool usepow;
        Nvector<Vec,4> derivs(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
        flt dudcosthetaCOS(const flt costheta) const;
    public:
        dihedral(const vector<flt> cosvals, 
                const vector<flt> sinvals = vector<flt>(),
                bool usepow = true);
        static flt getcos(const Vec& diff1, const Vec& diff2, const Vec& diff3);
        static flt getang(const Vec& diff1, const Vec& diff2, const Vec& diff3);
        
        flt dudcostheta(const flt theta) const;
        
        inline flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const {
            return energy(getang(diff1, diff2, diff3));
        };
        flt energy(flt ang) const;
        Nvector<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
};
#endif

class electricScreened : public interactpair {
    protected:
        flt screen;
        flt q1, q2;
        flt cutoff;
        flt cutoffE;
    public:
        electricScreened(const flt screenLength, const flt q1, 
            const flt q2, const flt cutoff);
        inline flt energy(const Vec r){return energy(r.mag(),q1*q2,screen, 0) - cutoffE;};
        static flt energy(const flt r, const flt qaqb, const flt screen, const flt cutoff=0);
        inline Vec forces(const Vec r){return forces(r,q1*q2,screen, cutoff);};
        static Vec forces(const Vec r, const flt qaqb, const flt screen, const flt cutoff=0);
};

//~ class interactgroup : public interaction {
    //~ protected:
        //~ vector<interaction*> inters;
    //~ public:
        //~ interactgroup(vector<interaction*> inters=vector<interaction*>())
                    //~ : inters(inters){};
        //~ void add(interaction* a){inters.push_back(a);};
        //~ uint size() const{ return inters.size();};
        //~ flt energy(Box *box);
        //~ void setForces(Box *box);
//~ };

struct fixedForceAtom {
    Vec F;
    atomid a;
    fixedForceAtom(Vec F, atomid a) : F(F), a(a) {};
    flt energy(Box &box){return -F.dot(a->x);};
    void setForce(Box &box){a->f += F;};
};

class fixedForce : public interaction {
    protected:
        vector<fixedForceAtom> atoms;
    public:
        fixedForce(vector<fixedForceAtom> atoms = vector<fixedForceAtom>()) : atoms(atoms){};
        void add(fixedForceAtom a){atoms.push_back(a);};
        void add(Vec F, atomid a){add(fixedForceAtom(F,a));};
        #ifdef VEC3D
        void add(flt x, flt y, flt z, atomid a){add(fixedForceAtom(Vec(x,y,z),a));};
        #elif defined VEC2D
        void add(flt x, flt y, atomid a){add(fixedForceAtom(Vec(x,y),a));};
        #endif
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box &box){
            flt E=0;
            for(vector<fixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<fixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};


struct fixedForceRegionAtom : public atomid {
    Vec direction;          // will be normalized
    vector<flt> boundaries; // should be 1 smaller than Fs
                            // each boundary is at (direction * b), 
                            // where b is in boundaries
    vector<flt> Fs;
    fixedForceRegionAtom(atomid a, Vec direction, vector<flt> boundaries, vector<flt> Fs);
    flt energy(Box &box);
    void setForce(Box &box);
};

class fixedForceRegion : public interaction {
    
    protected:
        vector<fixedForceRegionAtom> atoms;
    public:
        fixedForceRegion(vector<fixedForceRegionAtom> atoms = vector<fixedForceRegionAtom>()) 
            : atoms(atoms){};
        
        void add(fixedForceRegionAtom a){atoms.push_back(a);};
        void add(atomid a, Vec dir, vector<flt> bound, vector<flt> F){
            add(fixedForceRegionAtom(a, dir, bound, F));};
        uint size() const{ return (uint) atoms.size();};
        
        flt energy(Box &box){
            flt E=0;
            for(vector<fixedForceRegionAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<fixedForceRegionAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};

struct fixedSpringAtom {
    Vec loc;
    flt k;
    bool usecoord[3];
    atomid a;
    fixedSpringAtom(atomid a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true) : 
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
            return k*diffx.sq()/2;
        };
    void setForce(Box &box){
        Vec diffx = a->x - loc;
        for(uint i=0; i<2; ++i){
            if(!usecoord[i]) diffx[i] = 0;
        }
        a->f -= diffx * k;
    };
};

class fixedSpring : public interaction{
    protected:
        vector<fixedSpringAtom> atoms;
    public:
        fixedSpring(vector<fixedSpringAtom> atoms = vector<fixedSpringAtom>()) : atoms(atoms){};
        void add(fixedSpringAtom a){atoms.push_back(a);};
        void add(atomid a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true){
            add(fixedSpringAtom(a, loc, k, usex, usey, usez));};
        //void add(Vec F, atomid a){add(fixedForceAtom(F,a));};
        //void add(flt x, flt y, flt z, atomid a){add(fixedForceAtom(Vec(x,y,z),a));};
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box &box){
            flt E=0;
            for(vector<fixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                E += it->energy(box);
            return E;
        };
        void setForces(Box &box){
            for(vector<fixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); ++it)
                it->setForce(box);
        };
        flt pressure(Box &box){return NAN;};
};

class COMSpring : public interaction{
    protected:
        atomgroup *g1;
        atomgroup *g2;
        flt k, x0;
        flt m1, m2;
    public:
        COMSpring(atomgroup *g1, atomgroup *g2, flt k, flt x0=0) : 
            g1(g1), g2(g2), k(k), x0(x0), m1(g1->mass()), m2(g2->mass()){};
        flt energy(Box &box){
            flt dx = (g1->com() - g2->com()).mag() - x0;
            return k/2 * dx * dx;
        };
        void setForces(Box &box){
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.mag();
            flt fmag = -k * (comdist - x0);
            Vec a1 = comvec * (fmag / m1 / comdist);
            for(uint i=0; i < g1->size(); ++i){
                atom &atm = g1->get(i);
                atm.f += a1 * atm.m;
            }
            
            Vec a2 = comvec * (-fmag / m2 / comdist);
            for(uint i=0; i < g2->size(); ++i){
                atom &atm = g2->get(i);
                atm.f += a2 * atm.m;
            }
        };
        flt pressure(Box &box){
            //~ return NAN;
            // I think this is right, but I haven't checked it
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.mag();
            flt fmag = -k * (comdist - x0);
            
            Vec a12 = comvec * (fmag / m1 / m2 / comdist);
            
            flt P = 0;
            for(uint i=0; i < g1->size(); ++i){
                for(uint j=0; j < g2->size(); ++j){
                    atom &atm1 = g1->get(i);
                    atom &atm2 = g2->get(i);
                    Vec fij = a12 * (atm1.m * atm2.m);
                    Vec rij = box.diff(atm1.x, atm2.x);
                    P += fij.dot(rij);
                }
            }
            return P;
        };
};

enum BondDiffType {
    BOXED, // use box.diff(r1, r2)
    UNBOXED, // use r2 - r1
    FIXEDBOX // use r2 - r1 - (original box separation)
};

struct bondgrouping {
    flt k, x0;
    atomid a1, a2;
    BondDiffType diff_type;
    array<int,NDIM> fixed_box;
    bondgrouping(flt k, flt x0, atomid a1, atomid a2, 
            BondDiffType diff=UNBOXED, OriginBox *box=NULL);
    Vec diff(Box &box) const;
    int get_fixed(uint i){return fixed_box[i];};
    inline bool same_atoms(bondgrouping &other){
        return ((a1 == other.a1) and (a2 == other.a2)) or ((a1 == other.a2) and (a2 == other.a1));
    };
};

class bondpairs : public interaction {
    protected:
        bool zeropressure;
        vector<bondgrouping> pairs;
        //inline static Vec diff(Box &box, Vec r1, Vec r2){return r1-r2;};
        //inline static Vec diff(Box &box, Vec r1, Vec r2){return box.diff(r1, r2);};
    public:
        bondpairs(vector<bondgrouping> pairs, bool zeropressure=true);
        bondpairs(bool zeropressure=true);
        /// Add a pair of atoms.
        /// If "replace", a previous pair found will be replaced by the new pair.
        /// If not "replace" and that pair of atoms is already inserted, an error will be thrown.
        bool add(bondgrouping b, bool replace=true);
        inline bool add(flt k, flt x0, atomid a1, atomid a2, bool replace=true){
            return add(bondgrouping(k,x0,a1,a2), replace);};
        void add_forced(bondgrouping b){pairs.push_back(b);};
        /// Add a pair of atoms with the current distance.
        inline bool add(flt k, atomid a1, atomid a2, bool replace=true){
            flt x0 = (a1->x - a2->x).mag();
            return add(bondgrouping(k,x0,a1,a2), replace);
        };
        
        uint size() const{ return (uint) pairs.size();};
        bondgrouping get(uint i) const{ return pairs[i];};
        flt mean_dists(Box &box) const;
        flt std_dists(Box &box) const;
        flt energy(Box &box);
        void setForces(Box &box);
        flt pressure(Box &box);
        flt setForcesGetPressure(Box &box);
};

struct anglegrouping {
    flt k, x0;
    atomid a1, a2, a3;
    anglegrouping(flt k, flt x0, atomid a1, atomid a2, atomid a3) : 
                k(k),x0(x0), a1(a1), a2(a2), a3(a3){};
    inline bool same_atoms(anglegrouping &other){
        if(a2 != other.a2) return false;
        return ((a1 == other.a1) and (a3 == other.a3)) or ((a1 == other.a3) and (a3 == other.a1));
    };
};

class angletriples : public interaction {
    protected:
        vector<anglegrouping> triples;
        inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    public:
        angletriples(vector<anglegrouping> triples = vector<anglegrouping>());
        /// Add a triple of atoms.
        /// If "replace", a previous triple found will be replaced by the new triple.
        /// If not "replace" and that triple of atoms is already inserted, an error will be thrown.
        bool add(anglegrouping b, bool replace=true);
        inline bool add(flt k, flt x0, atomid a1, atomid a2, atomid a3, bool replace=true){
            return add(anglegrouping(k,x0,a1,a2,a3), replace);};
        /// Add a triple of atoms with the current angle.
        bool add(flt k, atomid a1, atomid a2, atomid a3, bool replace=true);
        void add_forced(anglegrouping b){triples.push_back(b);};
        
        flt energy(Box &box);
        inline flt pressure(Box &box){return 0;};
        void setForces(Box &box);
        uint size() const {return (uint) triples.size();};
        flt mean_dists() const;
        flt std_dists() const;
};

#ifdef VEC3D
struct dihedralgrouping {
    inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    dihedral dih;
    atomid a1, a2, a3, a4;
    dihedralgrouping(vector<flt> coscoeffs, vector<flt> sincoeffs,
                atomid a1, atomid a2, atomid a3, atomid a4, bool usepow=true) : 
                dih(coscoeffs, sincoeffs, usepow), a1(a1), 
                a2(a2), a3(a3), a4(a4){};
};

class dihedrals : public interaction {
    protected:
        vector<dihedralgrouping> groups;
    public:
        dihedrals(vector<dihedralgrouping> pairs = vector<dihedralgrouping>());
        void add(dihedralgrouping b){groups.push_back(b);};
        inline void add(vector<flt> nums, atomid a1, atomid a2, atomid a3, atomid a4){
            add(dihedralgrouping(nums, vector<flt>(), a1,a2,a3,a4));};
        inline void add(vector<flt> coscoeffs, vector<flt> sincoeffs, 
                            atomid a1, atomid a2, atomid a3, atomid a4, bool usepow=true){
            add(dihedralgrouping(coscoeffs, sincoeffs,a1,a2,a3,a4, usepow));};
        /// Add 4 atoms with the potential $V(\theta) = k (1 + \cos(\theta - \theta_0))$
        inline void add(flt k, flt theta0, atomid a1, atomid a2, atomid a3, atomid a4){
            vector<flt> coscoeffs(2,k);
            coscoeffs[1] = -k*cos(theta0);
            vector<flt> sincoeffs(2,0);
            sincoeffs[1] = -k*sin(theta0);
            add(dihedralgrouping(coscoeffs, sincoeffs, a1,a2,a3,a4, true));
        }
        /// Add 4 atoms with the potential $V(\theta) = k (1 + \cos(\theta - \theta_0))$, where
        /// \theta_0 is set to the current one
        inline void add(flt k, atomid a1, atomid a2, atomid a3, atomid a4){
            Vec r1 = dihedralgrouping::diff(a2->x, a1->x);
            Vec r2 = dihedralgrouping::diff(a3->x, a2->x);
            Vec r3 = dihedralgrouping::diff(a4->x, a3->x);
            add(k,  dihedral::getang(r1, r2, r3), a1, a2, a3, a4);
        };
        uint size() const{ return uint(groups.size());};
        flt mean_dists() const;
        //~ flt std_dists() const;
        flt energy(Box &box);
        void setForces(Box &box);
        inline flt pressure(Box &box){return 0;};
};
//~ class LJgroup : private vector<LJdata>, public atomgroup {
    //~ public:
        //~ LJgroup() : vector<LJdata>(){};
        //~ add(flt sigma, flt epsilon, atomid a){ 
//~ };
#endif


////////////////////////////////////////////////////////////////////////
struct forcepair {
    atom a1, a2;
    Vec fij;
};

struct forcepairx {
    atomid a1, a2;
    flt xij;
    Vec fij;
}; 

//
class fpairxFunct {
    public:
        virtual void run(forcepairx*) = 0;
        virtual ~fpairxFunct() = 0;
};

inline fpairxFunct::~fpairxFunct() { }  // defined even though it's pure virtual; it's faster this way; trust me

class interactionpairsx : public interaction {
    public:
        using interaction::setForces;
        virtual void setForces(Box &box, fpairxFunct*)=0;
        virtual ~interactionpairsx(){};
};

//~ struct atompaircomp {
    //~ bool operator() (const atompair& lhs, const atompair& rhs) const{
        //~ return (lhs[0] == rhs[0] and lhs[1] == rhs[1]);}
//~ };

struct Charged : public atomid {
    flt q;
    Charged() : atomid(), q(0){};
    Charged(flt q, atomid a) : atomid(a), q(q){};
};

struct ChargePair {
    flt q1q2;
    atomid atom1, atom2;
    ChargePair(Charged a1, Charged a2) : q1q2(a1.q*a2.q){};
};

////////////////////////////////////////////////////////////////////////
// Repulsive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// Potential is V(r) = ε (σ⁶/r⁶ - 1)²
// cutoff at sigma
struct LJatom : public atomid {
    flt epsilon, sigma;
    LJatom(){};
    LJatom(flt epsilon, flt sigma, atomid a) : atomid(a), 
            epsilon(epsilon), sigma(sigma){};
    LJatom(atomid a, LJatom other) : atomid(a),
                    epsilon(other.epsilon), sigma(other.sigma){};
    flt maxsize(){return sigma;};
};

struct LJpair {
    flt epsilon, sigma;
    atomid atom1, atom2;
    LJpair(LJatom LJ1, LJatom LJ2) :
            epsilon(sqrt(LJ1.epsilon * LJ2.epsilon)),
            sigma((LJ1.sigma + LJ2.sigma) / 2),
            atom1(LJ1), atom2(LJ2){};
    inline flt energy(Box &box){
        return LJrepulsive::energy(box.diff(atom1->x,atom2->x), epsilon, sigma);};
    inline Vec forces(Box &box){
        return LJrepulsive::forces(box.diff(atom1->x,atom2->x), epsilon, sigma);};
};

////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// cutoff at some sigcut
// Minimum at σ
struct LJatomcut : public LJatom {
    flt sigcut; // sigma units
    LJatomcut(){};
    LJatomcut(flt epsilon, flt sigma, atomid a, flt cut) : 
            LJatom(epsilon, sigma, a), sigcut(cut){};
    LJatomcut(atomid a, LJatomcut other) : LJatom(a, other), 
        sigcut(other.sigcut){};
    flt maxsize(){return sigma*sigcut;};
};



////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = ε₁₂ and σ = σ₁₂ (both indexed)
// cutoff at some sigcut (and r < σ)

struct LJAtomIndexed : public atomid {
    vector<flt> epsilons; // for finding epsilons 
    vector<flt> sigmas; // for finding epsilons 
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two HydroAtoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]; same for sigma
    flt sigcut; // sigma units
    LJAtomIndexed(){};
    LJAtomIndexed(vector<flt> epsilons, vector<flt> sigmas, uint indx, atomid a, flt cut) : 
            atomid(a), epsilons(epsilons), sigmas(sigmas), indx(indx),
             sigcut(cut){
                 assert(sigmas.size() == epsilons.size());
            };
    LJAtomIndexed(atomid a, LJAtomIndexed other) : atomid(a), epsilons(other.epsilons),
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
    atomid atom1, atom2;
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
    LJattractCut inter;
    atomid atom1, atom2;
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

struct HydroAtom : public atomid {
    vector<flt> epsilons; // for finding epsilons 
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two HydroAtoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]
    flt sigma;
    flt sigcut; // sigma units
    HydroAtom(){};
    HydroAtom(vector<flt> epsilons, uint indx, flt sigma, atomid a, flt cut) : 
            atomid(a), epsilons(epsilons), indx(indx), sigma(sigma),
             sigcut(cut){};
    HydroAtom(atomid a, HydroAtom other) : atomid(a), epsilons(other.epsilons),
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
    LJattractCut inter;
    atomid atom1, atom2;
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

struct LJishAtom : public atomid {
    vector<flt> epsilons; // for finding epsilons 
    flt repeps, sigma;
    flt exponent; // power
    uint indx; // for finding this one in other atoms' epsilon lists
    // note that for two atoms a1, a2, the epsilon for the pair
    // is then either a1.epsilons[a2.indx] or a2.epsilons[a1.indx]
    // these should be the same
    flt sigcut; // sigma units
    LJishAtom(){};
    LJishAtom(atomid a, vector<flt> epsilons, flt repeps, flt sigma,
                            flt n, uint indx, flt cut) : 
            atomid(a), epsilons(epsilons), repeps(repeps), sigma(sigma), 
            exponent(n), indx(indx), sigcut(cut){
        assert(indx < epsilons.size());
        //~ assert(sigma > 0.01);
    };
    LJishAtom(atomid a, LJishAtom other) : 
        atomid(a), epsilons(other.epsilons), repeps(other.repeps), 
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
    atomid atom1, atom2;
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
        flt rsq = rij.sq()/(sigma*sigma);
        if(rsq > cutR * cutR){
            //~ printf("LJish: dx=%.2f, σ=%.2f, rsq=%.2f, cutR²=%.2f\n", rij.mag(), sigma, rsq, cutR);
            return 0;
        }
        
        flt mid = (1-pow(rsq,-n/2));
        //~ printf("LJish: rsq=%.2f, mid %.2f, epsilon %.2f, repeps %.2f\n", rsq, mid, epsilon, repeps);
        if (rsq > 1) return epsilon*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
    };
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        flt rsq = dsq/(sigma*sigma);
        if(rsq > cutR * cutR) return Vec();
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


struct LJAttractRepulseAtom : public atomid {
    vector<flt> epsilons; // for finding epsilons 
    flt sig;
    uint indx;
    flt sigcut; // sigma units
    LJAttractRepulseAtom(){};
    LJAttractRepulseAtom(atomid a, vector<flt> epsilons, flt sigma, uint indx, flt cut) : 
            atomid(a), epsilons(epsilons), sig(sigma), indx(indx),
             sigcut(cut){
                 assert(indx < epsilons.size());
                 //~ printf("Made atom  with ε=%.2f of size σ=%.2f\n", epsilons[indx], sig);
                 };
    LJAttractRepulseAtom(atomid a, LJAttractRepulseAtom other) : atomid(a), epsilons(other.epsilons),
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
    atomid atom1, atom2;
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
        flt rsq = rij.sq()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrt(rij.sq()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-pow(rsq,-3));
        //~ if(eps*(mid*mid) - cutE < 0) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                //~ rij.mag(), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        //~ }
        return eps*(mid*mid) - cutE;
    };
    inline Vec forces(Box &box){
        if(eps == 0) return Vec();
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec();
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

struct LJAttractFixedRepulseAtom : public atomid {
    vector<flt> epsilons; // for finding epsilons 
    flt repeps, sig;
    uint indx;
    flt sigcut; // sigma units
    LJAttractFixedRepulseAtom(){};
    LJAttractFixedRepulseAtom(atomid a, vector<flt> epsilons, flt repeps, 
            flt sigma, uint indx, flt cut) : 
            atomid(a), epsilons(epsilons), repeps(repeps), sig(sigma), 
            indx(indx), sigcut(cut){
                 assert(indx < epsilons.size());
                 };
    LJAttractFixedRepulseAtom(atomid a, LJAttractFixedRepulseAtom other) : 
        atomid(a), epsilons(other.epsilons),
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
    atomid atom1, atom2;
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
        flt rsq = rij.sq()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrt(rij.sq()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-pow(rsq,-3)); // # 1 - σ⁶/r⁶
        //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrt(rij.sq()), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        if (rsq > 1) return eps*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
        //~ flt E;
        //~ if (rsq > 1) E = eps*(mid*mid) - cutE;
        //~ else E = repeps*(mid*mid) - cutE;
        //~ if(E > 1e4)
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrt(rij.sq()), E, eps, repeps, sig, cutR, cutE);
        //~ return E;
    };
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec();
        flt rsix = pow(rsq,-3); // σ⁶/r⁶
        flt fmagTimesR = 12*rsix*(rsix - 1);
        if (rsq < 1) return rij * (repeps * fmagTimesR / dsq);
        return rij * (eps * fmagTimesR / dsq);
        //~ flt fmag;
        //~ if (rsq < 1) fmag = repeps * fmagTimesR / dsq;
        //~ else fmag = eps * fmagTimesR / dsq;
        //~ if(fmag * rij.mag() > 1e4)
            //~ printf("Distance: %.2f Force: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrt(rij.sq()), fmag * rij.mag(), eps, repeps, sig, cutR, cutE);
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
    LJDoubleAtom(flt epsilon, flt epsrep, flt sigma, atomid a, flt cut) : 
            LJatom(epsilon, sigma, a), epsrep(epsrep), sigcut(cut){};
    LJDoubleAtom(atomid a, LJDoubleAtom other) : LJatom(a, other), 
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

struct EisMclachlanAtom : public atomid {
    flt dist, sigmai; // dist is radius of atom + radius of water (1.4 Å)
    EisMclachlanAtom(){};
    EisMclachlanAtom(flt dist, flt sigmai, atomid a) : atomid(a),
            dist(dist), sigmai(sigmai){};
    EisMclachlanAtom(atomid a, EisMclachlanAtom other) : atomid(a), 
        dist(other.dist), sigmai(other.sigmai){};
    flt maxsize(){return dist;};
};

struct EisMclachlanPair {
    flt c0,c1,c2; // 1/R term, const term, R term in E (coefficients)
    // E = c₀/R + c₁ + c₂ R
    flt cutoff; // r1 + r2 (includes two water radii)
    atomid atom1, atom2;
    EisMclachlanPair(EisMclachlanAtom a1, EisMclachlanAtom a2) : 
        c0(-M_PI*(a1.sigmai*a2.dist - a2.sigmai*a1.dist)
                * (a1.dist*a1.dist - a2.dist*a2.dist)),
        c1(-2*M_PI*(a1.sigmai*a2.dist*a2.dist + a2.sigmai*a1.dist*a1.dist)),
        c2(M_PI*(a1.sigmai*a2.dist + a2.sigmai*a1.dist)),
        cutoff(a1.dist + a2.dist),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt R = rij.mag();
        if (R > cutoff) return 0;
        return c0/R + c1 + c2*R;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > (cutoff*cutoff)) return Vec();
        flt R = sqrt(dsq);
        return rij * ((c0/dsq-c2)/R);
    }
};

////////////////////////////////////////////////////////////////////////
// Hertzian potential, with ε = √(ε₁ ε₂) and σ = (σ₁ + σ₂)/2
// Potential is V(r) = ε/n (1 - r/σ)^n, with n = 5/2 usually
// cutoff at r = σ

struct HertzianAtom : public atomid {
    flt eps, sigma, exponent;
    HertzianAtom(){};
    HertzianAtom(atomid a, flt eps, flt sigma, flt exponent=2.5) : atomid(a),
            eps(eps), sigma(sigma), exponent(exponent){};
    HertzianAtom(atomid a, HertzianAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), exponent(other.exponent){};
    flt maxsize(){return sigma;};
};

inline HertzianAtom hertzd(atomid a, double eps, double sigma, double exponent=2.5){
    return HertzianAtom(a, eps, sigma, exponent);
};

struct EnergyForce {
    Vec f;
    flt E;
    EnergyForce(Vec f, flt E) : f(f), E(E){};
};

struct HertzianPair {
    flt eps, sig, exponent;
    atomid atom1, atom2;
    HertzianPair(HertzianAtom a1, HertzianAtom a2) : 
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        exponent((a1.exponent + a2.exponent)/2.0), atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > sig*sig) return 0.0;
        flt R = sqrt(dsq);
        return eps * pow(1.0 - (R/sig), exponent) / exponent;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > sig*sig) return Vec();
        flt R = sqrt(dsq);
        return rij * (eps * pow(1.0 - (R/sig), exponent-1) /sig/R);
    }
    inline EnergyForce EnergyForces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > sig*sig) return EnergyForce(Vec(),0);
        flt R = sqrt(dsq);
        
        Vec f = rij * (eps * pow(1.0 - (R/sig), exponent-1) /sig/R);
        flt E = eps * pow(1.0 - (R/sig), exponent) / exponent;
        return EnergyForce(f, E);
    }
    //~ inline flt xrij(Box &box){
        //~ Vec rij = box.diff(atom1->x, atom2->x);
        //~ flt dsq = rij.sq();
        //~ if(dsq > sig*sig) return 0.0;
        //~ flt R = sqrt(dsq);
        //~ return (R*eps*(exponent-1)/sig/sig) * pow(1.0 - (R/sig), exponent-2);
    //~ }
    inline void fill(Box &box, forcepairx &fpair){
        fpair.a1 = atom1;
        fpair.a2 = atom2;
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > sig*sig){
            fpair.fij = Vec();
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

struct HertzianDragAtom : public atomid {
    flt eps, sigma, exponent, gamma;
    HertzianDragAtom(){};
    HertzianDragAtom(atomid a, flt eps, flt sigma, flt gamma, flt exponent=2.5) : atomid(a),
            eps(eps), sigma(sigma), exponent(exponent), gamma(gamma){};
    HertzianDragAtom(atomid a, HertzianDragAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), exponent(other.exponent), gamma(other.gamma){};
    flt maxsize(){return sigma;};
};

struct HertzianDragPair {
    flt eps, sig, exponent, gamma;
    atomid atom1, atom2;
    HertzianDragPair(HertzianDragAtom a1, HertzianDragAtom a2) : 
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        exponent((a1.exponent + a2.exponent)/2.0), 
        gamma((a1.gamma + a2.gamma)/2), atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
        if(dsq > sig*sig) return 0.0;
        flt R = sqrt(dsq);
        return eps * pow(1.0 - (R/sig), exponent) / exponent;
    }
    inline Vec forces(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        Vec vij = atom1->v - atom2->v;
        flt dsq = rij.sq();
        if(dsq > sig*sig) return Vec();
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

struct LoisOhernAtom : public atomid {
    flt eps, sigma, C, l;
    LoisOhernAtom(){};
    LoisOhernAtom(atomid a, flt eps, flt sigma, flt C, flt l) : atomid(a),
            eps(eps), sigma(sigma), C(C), l(l){};
    LoisOhernAtom(atomid a, LoisOhernAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), C(other.C), l(other.l){};
    flt maxsize(){return sigma*(1+C+l);};
};

struct LoisOhernPair {
    flt eps, sig, C, l, sigcut;
    atomid atom1, atom2;
    LoisOhernPair(LoisOhernAtom a1, LoisOhernAtom a2) : 
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        C((a1.C + a2.C)/2.0), l((a1.l + a2.l)/2.0), sigcut(sig*(1+C+l)),
        atom1(a1), atom2(a2){};
    LoisOhernPair(LoisOhernAtom a1, LoisOhernAtom a2, flt eps, flt sig, flt C, flt l) : 
        eps(eps), sig(sig), C(C), l(l), sigcut(sig*(1+C+l)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
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
        flt dsq = rij.sq();
        if(dsq >= sigcut*sigcut) return Vec();
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
 * epsilon is the strength of the harmonic interaction
 * f is the fixed-force amount (real units)
 * l is the width of the fixed-force region (real units)
 * To combine, we use a geometric average of epsilon, and everything else
 * is averaged.
 */

struct LoisLinAtom : public atomid {
    flt eps, sigma, f, l;
    LoisLinAtom(){};
    LoisLinAtom(atomid a, flt eps, flt sigma, flt depth, flt width) : atomid(a),
            eps(eps), sigma(sigma), f(width > 0 ? depth/width : 0), l(width){};
    LoisLinAtom(atomid a, LoisLinAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), f(other.f), l(other.l){};
    flt maxsize(){return sigma*(1+l);};
};

struct LoisLinPair {
    flt eps, sig, f, l, sigcut;
    atomid atom1, atom2;
    LoisLinPair(LoisLinAtom a1, LoisLinAtom a2) : 
        eps(sqrt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        f((a1.f + a2.f)/2.0), l((a1.l + a2.l)/2.0), sigcut(sig +l),
        atom1(a1), atom2(a2){};
    LoisLinPair(LoisLinAtom a1, LoisLinAtom a2, flt eps, flt sig, flt f, flt l) : 
        eps(eps), sig(sig), f(f), l(l), sigcut(sig+l),
        atom1(a1), atom2(a2){};
    inline flt energy(Box &box){
        Vec rij = box.diff(atom1->x, atom2->x);
        flt dsq = rij.sq();
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
        flt dsq = rij.sq();
        if(dsq >= sigcut*sigcut) return Vec();
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
class LJsimple : public interaction {
    protected:
        vector<LJatom> atoms;
        pairlist ignorepairs;
        atomid get_id(atom* a);

    public:
        LJsimple(flt cutoffdist, vector<LJatom> atms=vector<LJatom>());
         // cutoffdist in sigma units
        
        inline void add(LJatom a){
            assert(a.n() <= atoms.size());
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            atoms[a.n()] = a;};
        inline void add(atomid a, flt epsilon, flt sigma){
            add(LJatom(epsilon,sigma,a));};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        inline void ignore(atom* a, atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint atoms_size() const{return (uint) atoms.size();};
        flt energy(Box &box);
        flt pressure(Box &box);
        void setForces(Box &box);
        //~ ~LJsimple(){};
};

template <class A, class P>
class SCboxed : public interaction {
    protected:
        sptr<SCbox> box;
        sptr<atomvec> atoms;
        vector<A> group; // data on which atoms interact with the wall
    public:
        SCboxed(sptr<atomvec> atomv, sptr<SCbox> box)
            : box(box), atoms(atomv){};
        inline void add(A atm){group.push_back(atm);};
        flt energy(Box &box);
        flt pressure(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        inline vector<A> &atom_list(){return group;};
};

template <class A, class P>
class SimpleListed : public interaction {
    protected:
        vector<A> atoms;
    
    public:
        SimpleListed(){};
        inline void add(A atm){atoms.push_back(atm);};
        //flt energy(Box &box, idpair &pair);
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
class NListed : public interaction {
    // neighborlist keeps track of atom* pairs within a particular distance.
    // NListed keeps track of atoms with additional properties
    // that interact through a particular interaction.
    // class A needs to inherit from atomid, and also have an initialization
    // method A(atomid a, A other)
    // you also need to implement a couple methods below
    
    // Implementation detail:
    // note that neighborlist maintains its own metagroup, so that when
    // a member of A is created and passed to this group, the atomid in that
    // member has an A.n() value referring to its place in the *original*
    // atomgroup, not this one. This is why we need an A(atomid, A) method;
    // so we can make a new A with all the same properties as before,
    // but with the n() referring to the neighborlist maintained atomgroup.
    protected:
        vector<A> atoms; // NOT all full. 
                         // This is the length of the atomvec, and if 
                         // atom i has been added, then atoms[i] is 
                         // filled in this vector
        vector<P> pairs;
        sptr<neighborlist> neighbors;
        uint lastupdate;
    public:
        NListed(sptr<atomvec> vec, sptr<neighborlist> neighbors) : 
            atoms(vec->size()), neighbors(neighbors), lastupdate(0){}; //group(vec), 
        inline void add(A atm){
            neighbors->add(atm, atm.maxsize());
            // assert(atoms.size() > atm.n());
            atoms[atm.n()] = atm;
        }
        NListed(sptr<Box> box, sptr<atomvec> atomv, const flt skin) : 
            atoms(atomv->size()),
            neighbors(new neighborlist(box, atomv, skin)), lastupdate(0){};
        void update_pairs();
        P getpair(idpair &pair){
            return P(atoms[pair.first().n()], atoms[pair.last().n()]);}
        A& getatom(uint n){return atoms[n];}
        flt energy(Box &box, idpair &pair);
        flt energy(Box &box);
        flt pressure(Box &box);
        inline vector<P> &pairiter(){return pairs;};
        uint size(){return ((uint) (atoms.size()));};
        inline flt energy_pair(P pair, Box &box){return pair.energy(box);}; // This may need to be written!
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        inline Vec forces_pair(P pair, Box &box){return pair.forces(box);}; // This may need to be written!
        inline vector<A> &atom_list(){return atoms;};
        inline sptr<neighborlist> nlist(){return neighbors;};
        //~ flt energy_test(flt dist);
        //~ flt force_test(flt dist);
        ~NListed(){};
};

template <class A, class P>
class NListedVirial : public interactionpairsx {
    private:
        NListed<A,P> nlisted;
    public:
        NListedVirial(sptr<atomvec> vec, sptr<neighborlist> neighbors) :
                nlisted(vec, neighbors){};
        void setForces(Box &box){nlisted.setForces(box);};
        void setForces(Box &box, fpairxFunct*);
        virtual inline flt setForcesGetPressure(Box &box){return nlisted.setForcesGetPressure(box);};
        virtual flt setForcesGetEnergy(Box &box);
        virtual inline flt energy(Box &box){return nlisted.energy(box);};
        virtual inline flt pressure(Box &box){return nlisted.pressure(box);};
        inline void add(A atm){nlisted.add(atm);}
        inline sptr<neighborlist> nlist(){return nlisted.nlist();};
};

template <class A, class P>
flt SCboxed<A, P>::energy(Box &newbox){
    flt E = 0;
    atom wallatom;
    wallatom.m = NAN;
    atomid wallid = atomid(&wallatom, atoms->size());
    typename vector<A>::iterator it;
    for(it = group.begin(); it != group.end(); ++it){
        Vec r0 = (*it)->x;
        Vec edger = box->edgedist(r0);
        wallatom.x = r0 + (edger*2);
        E += P(*it, HertzianAtom(wallid, *it)).energy(newbox) / 2; // no energy from fake atom
    }
    return E;
};

template <class A, class P>
void SCboxed<A, P>::setForces(Box &newbox){
    atom wallatom;
    wallatom.m = NAN;
    atomid wallid = atomid(&wallatom, atoms->size());
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
flt SCboxed<A, P>::setForcesGetPressure(Box &newbox){
    flt p = 0;
    atom wallatom;
    wallatom.m = NAN;
    atomid wallid = atomid(&wallatom, atoms->size());
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
flt SCboxed<A, P>::pressure(Box &newbox){
    flt p = 0;
    atom wallatom;
    wallatom.m = NAN;
    atomid wallid = atomid(&wallatom, atoms->size());
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
    vector<idpair>::iterator pairit;
    for(pairit = neighbors->begin(); pairit != neighbors->end(); ++pairit){
        // assert(atoms.size() > pairit->first().n());
        // assert(atoms.size() > pairit->last().n());
        A firstatom = atoms[pairit->first().n()];
        A secondatom = atoms[pairit->last().n()];
        pairs.push_back(P(firstatom, secondatom));
    }
}

template <class A, class P>
flt NListed<A, P>::energy(Box &box, idpair &pair){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    P Epair = P(atoms[pair.first().n()],atoms[pair.last().n()]);
    //~ Vec dist = box.diff(Epair.atom1->x, Epair.atom2->x);
    return energy_pair(Epair, box);
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
        //~ assert(f.sq() < 1000000);
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
    //~ atom a1 = atom(it->first());
    //~ atom a2 = atom(it->second());
    //~ A atm1 = A(it->first(), &a1);
    //~ A atm2 = A(it->second(), &a1);
    //~ a1.x = Vec();
    //~ a1.v = Vec();
    //~ a1.a = Vec();
    //~ a1.f = Vec();
    //~ 
    //~ a2.x = Vec();
    //~ a2.v = Vec();
    //~ a2.a = Vec();
    //~ a2.f = Vec();
    //~ a2.x.setx(dist);
    //~ P Epair = P(atm1,atm2);
    //~ return energy_pair(Epair, infbox);
//~ };

template <class A, class P>
void NListedVirial<A, P>::setForces(Box &box, fpairxFunct* funct){
    nlisted.update_pairs(); // make sure the LJpairs match the neighbor list ones
    forcepairx myfpair;
    typename vector<P>::iterator it;
    for(it = nlisted.pairiter().begin(); it != nlisted.pairiter().end(); ++it){
        it->fill(box, myfpair);
        //~ assert(!isnan(myfpair.fij[0]));
        funct->run(&myfpair);
        it->atom1->f += myfpair.fij;
        it->atom2->f -= myfpair.fij;
        //~ assert(f.sq() < 1000000);
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

class Charges : public interaction {
    protected:
        vector<Charged> atoms;
        pairlist ignorepairs;
        flt screen;
        flt k;
        atomid get_id(atom* a);

    public:
        Charges(flt screenlength, flt k=1, vector<Charged> atms=vector<Charged>());
         // cutoffdist in sigma units
        
        inline void add(Charged a){atoms.push_back(a); return;};
        inline void add(atomid a, flt q){add(Charged(q,a));};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        inline void ignore(atom* a, atom* b){
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

struct WallAtom : atomref {
    public:
        flt sigma;
        flt epsilon;
    public:
        WallAtom(atomid a, flt sigma, flt epsilon=1.0) : 
            atomref(a), sigma(sigma), epsilon(epsilon){};
};

class SoftWall : public interaction {
    protected:
        Vec loc;
        Vec norm;
        flt expt;
        vector<WallAtom> group;
        flt lastf;
    public:
        SoftWall(Vec loc, Vec norm, flt expt=2.0) : 
            loc(loc), norm(norm.norm()), expt(expt), lastf(NAN){};
        void add(WallAtom a){group.push_back(a);};
        flt energy(Box &box);
        void setForces(Box &box);
        flt setForcesGetPressure(Box &box);
        flt pressure(Box &box);
        
        void setLoc(Vec newloc){loc = newloc;};
        Vec getLoc(){return loc;};
        void setNorm(Vec newNorm){norm = newNorm.norm();};
        Vec getNorm(){return norm;};
        
        flt get_last_f(){return lastf;};
};

class SoftWallCylinder : public interaction {
    protected:
        Vec loc;
        Vec axis;
        flt radius;
        flt expt;
        flt lastf;
        vector<WallAtom> group;
    public:
        SoftWallCylinder(Vec loc, Vec axis, flt radius, flt expt=2.0) : 
            loc(loc), axis(axis.norm()), radius(radius), expt(expt), lastf(NAN){};
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
        void setAxis(Vec new_axis){axis = new_axis.norm();};
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
            return diff(v, Vec());
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

class SCatomvec : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms, which are
    // arranged in pairs.
    private:
        atomvec atoms;
    public:
        SCatomvec(vector<double> masses) : atoms((uint) masses.size()*2, 0.0){
            for(uint i=0; i < masses.size(); ++i){
                atoms[2*i].m = masses[i]/2;
                atoms[2*i+1].m = masses[i]/2;
            }
        };
        SCatomvec(uint N, flt mass) : atoms(N*2, mass/2.0){};
        inline atomvec &vec(){return atoms;};
        inline atom& operator[](cuint n){return atoms[n];};
        inline atom& operator[](cuint n) const {return atoms[n];};
        inline atomid get_id(cuint n){return atoms.get_id(n);};
        inline idpair pair(cuint n){return idpair(atoms.get_id(n*2), atoms.get_id(n*2 + 1));};
        inline uint size() const {return atoms.size();};
        inline uint pairs() const {return atoms.size()/2;};
        ~SCatomvec(){};
};

struct SpheroCylinderDiff{
    Vec delta, r;
    flt lambda1, lambda2;
};

struct SCPair {
    idpair &p1;
    idpair &p2;
    flt l1, l2;
    SCPair(idpair &p1, idpair &p2, flt l1, flt l2) : 
        p1(p1), p2(p2), l1(l1), l2(l2){};
    SCPair(idpair &p1, idpair &p2, flt l) : 
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
    
    SCSpringPair(idpair &p1, idpair &p2, flt eps, flt sig, flt l1, flt l2) : 
        SCPair(p1, p2, l1, l2), eps(eps), sig(sig){};
    SCSpringPair(idpair &p1, idpair &p2, flt eps, flt sig, flt l) : 
        SCPair(p1, p2, l), eps(eps), sig(sig){};
    
    inline flt maxdist(){return sig + (l1+l2)/2;};
    inline flt maxdelta(){return sig;};
    
    flt energy(Box &box, SpheroCylinderDiff diff){
        flt dsq = diff.delta.sq();
        if(dsq > sig*sig) return 0;
        flt d = sqrt(dsq);
        flt dsig = d-sig;
        return dsig*dsig*eps/2;
    }
    Vec forces(Box &box, SpheroCylinderDiff diff){
        flt dsq = diff.delta.sq();
        if(dsq > sig*sig) return Vec();
        flt dmag = sqrt(dsq);
        Vec dhat = diff.delta / dmag;
        
        return dhat * (eps * (sig - dmag));
    };
};

class SCSpringList : public interaction {
    /// Harmonic repulsive interactions between spherocylinders.
    private:
        SCatomvec *scs;
        flt eps, sig;
        vector<flt> ls;
        set<array<uint, 2> > ignore_list;
    public:
        SCSpringList(SCatomvec *scs, flt eps, flt sig, flt l) : 
            scs(scs), eps(eps), sig(sig), ls(scs->pairs(), l){};
        SCSpringList(SCatomvec *scs, flt eps, flt sig, vector<flt> ls) : 
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
        void ignore(atomid a1, atomid a2){ignore(a1.n() / 2, a2.n() / 2);}
        
        ~SCSpringList(){};
};

#endif
