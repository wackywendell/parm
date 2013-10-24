#include "vec.hpp"

#include <vector>
#include <bitset> // for contact tracking
#include <set>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <climits>

#ifndef INTERACTION_H
#define INTERACTION_H

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
typedef unsigned int uint;
typedef const unsigned int cuint;

class Box {
    public:
        virtual Vec diff(Vec r1, Vec r2)=0;
        virtual flt V()=0;
        virtual ~Box(){};
};

class interaction {
    public:
        virtual flt energy(Box *box)=0;
        virtual void setForces(Box *box)=0;
        virtual flt setForcesGetPressure(Box *box){return NAN;};
        virtual flt pressure(Box *box)=0;
        virtual ~interaction(){};
};

class statetracker {
    public:
        virtual void update(Box *box) = 0;
        virtual ~statetracker(){};
};

/***********************************************************************
 * Boxes
 */

#ifdef VEC3D
#ifdef LONGFLOAT
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remainderl(r1[0], r2[0]), remainderl(r1[1], r2[1]), remainderl(r1[2], r2[2]));
};
#else
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remainder(r1[0], r2[0]), remainder(r1[1], r2[1]), remainder(r1[2], r2[2]));
};
#endif
#endif
#ifdef VEC2D
#ifdef LONGFLOAT
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remainderl(r1[0], r2[0]), remainderl(r1[1], r2[1]));
};
#else
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remainder(r1[0], r2[0]), remainder(r1[1], r2[1]));
};
#endif
#endif


class InfiniteBox : public Box {
    public:
        Vec diff(Vec r1, Vec r2){return r1-r2;};
        flt V(){return NAN;};
};

InfiniteBox infbox;

class OriginBox : public Box {
    private:
        Vec boxsize;
    public:
        OriginBox(Vec size) : boxsize(size){};
        Vec diff(Vec r1, Vec r2){
            return vecmod((r1-r2), boxsize);
        }
        #ifdef VEC3D
        flt V(){return boxsize[0] * boxsize[1] * boxsize[2];};
        flt L(){return (boxsize[0] + boxsize[1] + boxsize[2])/3.0;};
        #endif
        #ifdef VEC2D
        flt V(){return boxsize[0] * boxsize[1];};
        flt L(){return (boxsize[0] + boxsize[1])/2.0;};
        #endif
        flt resize(flt factor){boxsize *= factor; return V();}
        flt resizeV(flt newV){flt curV = V(); boxsize *= powflt(newV/curV, OVERNDIM); return V();}
        Vec randLoc(){
            Vec v = randVecBoxed();
            for(uint i=0; i<NDIM; i++){
                v[i] *= boxsize[i];
            }
            return diff(v, Vec());
        };
};

/***********************************************************************
 * Atoms
 */
struct atom {
    Vec x; // location
    Vec v; // velocity
    Vec a; // acceleration
    Vec f; // forces
    flt m; // mass
};

class atomref {
    private:
        atom *ptr;
    public:
        inline atomref() : ptr(NULL){};
        inline atomref(atom *a) : ptr(a){};
        inline atom& operator *(){return *ptr;};
        inline atom* pointer(){return ptr;};
        inline Vec& x(){return ptr->x;};
        inline Vec& v(){return ptr->v;};
        inline Vec& f(){return ptr->f;};
        inline Vec& a(){return ptr->a;};
        inline flt& m(){return ptr->m;};
        inline bool operator==(const atomref &other) const {return other.ptr == ptr;};
        inline bool operator==(const atom* other) const {return other == ptr;};
        inline bool operator!=(const atomref &other) const {return other.ptr != ptr;};
        inline bool operator<(const atomref &other) const {return ptr < other.ptr;};
        inline bool operator<=(const atomref &other) const {return ptr <= other.ptr;};
        inline bool operator>=(const atomref &other) const {return ptr >= other.ptr;};
        inline bool operator>(const atomref &other) const {return ptr > other.ptr;};
};

class atomid : public atomref {
    private:
        uint num; // note that these are generally only in reference to 
                  // a specific atomgroup
    public:
        inline atomid() : atomref(), num(UINT_MAX){};
        inline atomid(atom *a) : atomref(a), num(UINT_MAX){};
        inline atomid(atom *a, uint n) : atomref(a), num(n){};
        inline uint n() const {return num;};
};

class idpair : public array<atomid, 2> {
    public:
        idpair(atomid a, atomid b){ vals[0] = a; vals[1] = b;};
        inline atomid first() const {return vals[0];};
        inline atomid last() const {return vals[1];};
};

//~ class atompair : public array<atom*, 2> {
    //~ public:
        //~ atompair(atom* a, atom* b){ vals[0] = a; vals[1] = b;};
        //~ inline atom& first() const {return *(vals[0]);};
        //~ inline atom& last() const {return *(vals[1]);};
//~ };

class atomgroup {
    // a group of atoms, such as a molecule, sidebranch, etc.
    public:
        // access individual atoms
        virtual atom& operator[](cuint n)=0;
        virtual atom& operator[](cuint n) const=0;
        virtual atom* get(cuint n){if(n>=size()) return NULL; return &((*this)[n]);};
        virtual atomid get_id(cuint n){return atomid(get(n),n);};
        virtual uint size() const=0;
        virtual flt getmass(const unsigned int n) const {return (*this)[n].m;};
        
        
        Vec com() const; //center of mass
        Vec comv() const; //center of mass velocity
        
        //Stats
        flt mass() const;
        flt kinetic(const Vec &originvelocity=Vec()) const;
        Vec momentum() const;
        #ifdef VEC3D
        flt moment(const Vec &loc, const Vec &axis, Box *box) const;
        Vec angmomentum(const Vec &loc, Box *box) const;
        Matrix<flt> moment(const Vec &loc, Box *box) const;
        Vec omega(const Vec &loc, Box *box) const;
        void addOmega(Vec w, Vec origin, Box *box);
        inline void resetL(Box *box){
            Vec c = com(), w = omega(c, box);
            if (w.sq() == 0) return;
            addOmega(-w, c, box);
        }
        #elif defined VEC2D
        flt moment(const Vec &loc, Box *box) const;
        flt angmomentum(const Vec &loc, Box *box) const;
        flt omega(const Vec &loc, Box *box) const{return angmomentum(loc, box) / moment(loc, box);};
        void addOmega(flt w, Vec origin, Box *box);
        inline void resetL(Box *box){
            Vec c = com();
            flt w = omega(c, box);
            if (w == 0) return;
            addOmega(-w, c, box);
        }
        #endif
        
        
        // for resetting
        void addv(Vec v);
        void resetcomv(){addv(-comv());};
        
        
        // for timestepping
        void resetForces();
        void setAccel();
        virtual ~atomgroup(){};
};

class atomvec : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms.
    private:
        atom* atoms;
        uint sz;
    public:
        atomvec(vector<double> masses) : sz((uint) masses.size()){
            atoms = new atom[sz];
            for(uint i=0; i < sz; i++) atoms[i].m = masses[i];
        };
        atomvec(uint N, flt mass) : sz(N){
            atoms = new atom[sz];
            for(uint i=0; i < sz; i++) atoms[i].m = mass;
        };
        atomvec(atomvec& other) : sz(other.size()){
            atoms = new atom[sz];
            for(uint i=0; i < sz; i++) atoms[i] = other.atoms[i];
        };
        inline atom& operator[](cuint n){return atoms[n];};
        inline atom& operator[](cuint n) const {return atoms[n];};
        atomid get_id(atom *a);
        inline atomid get_id(uint n) {
            if (n > sz) return atomid(); return atomid(atoms + n,n);};
        //~ inline flt getmass(cuint n) const{return atoms[n].m;};
        //~ inline void setmass(cuint n, flt m){atoms[n].m = m;};
        inline uint size() const {return sz;};
        ~atomvec(){ delete [] atoms;};
};

class metagroup : public atomgroup {
    protected:
        vector<atom*> atoms;
    public:
        metagroup(){};
        //metagroup(vector<atom*> atoms) : atoms(atoms){};
        metagroup(vector<atomgroup*>);
        inline atom& operator[](cuint n){return *atoms[n];};
        inline atom& operator[](cuint n) const{return *atoms[n];};
        inline atom* get(cuint n){if(n>=size()) return NULL; return (atoms[n]);};
        inline void add(atom *a){return atoms.push_back(a);};
        atomid get_id(atom *a);
        inline atomid get_id(uint n) {return atomid(atoms[n],n);};
        inline uint size() const {return (uint) atoms.size();};
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

static const flt LJr0 = powflt(2.0, 1.0/6.0);
static const flt LJr0sq = powflt(2.0, 1.0/3.0);

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
    protected:
        flt epsilon;
        flt sigma;
        flt cutR, cutE;
    public:
        LJFullCut(const flt epsilon, const flt sigma, const flt cutsig):
            epsilon(epsilon), sigma(sigma), cutR(cutsig){
                flt rsix = powflt(cutR, 6);
                flt mid = (1-1/rsix);
                cutE = epsilon*(mid*mid-1);
            };
        inline static flt energy(const Vec diff, const flt eps, const flt sig, const flt cutsig){
            flt rsq = diff.sq()/(sig*sig);
            if(rsq > (cutsig*cutsig)) return 0;
            flt rsix = rsq*rsq*rsq;
            //~ return eps*(4*(1/(rsix*rsix) - 1/rsix) + 1);
            flt mid = (1-1/rsix);
            return eps*(mid*mid-1);
        };
        inline static flt energy(const flt rsig, const flt cutsig){
            if(rsig > cutsig) return 0;
            flt rsq = rsig * rsig;
            flt rsix = rsq*rsq*rsq;
            flt mid = (1-1/rsix);
            return mid*mid-1;
            // Minimum of -1 occurs at r = 1
        };
        inline flt energy(const Vec& diff){return energy(diff, epsilon, sigma, cutR);};
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
        inline flt energy(const Vec r){return energy(r.mag(),q1*q2,screen, cutoff);};
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
    atom *a;
    fixedForceAtom(Vec F, atom *a) : F(F), a(a) {};
    flt energy(Box *box){return -F.dot(a->x);};
    void setForce(Box *box){a->f += F;};
};

class fixedForce : public interaction{
    protected:
        vector<fixedForceAtom> atoms;
    public:
        fixedForce(vector<fixedForceAtom> atoms = vector<fixedForceAtom>()) : atoms(atoms){};
        void add(fixedForceAtom a){atoms.push_back(a);};
        void add(Vec F, atom* a){add(fixedForceAtom(F,a));};
        #ifdef VEC3D
        void add(flt x, flt y, flt z, atom* a){add(fixedForceAtom(Vec(x,y,z),a));};
        #elif defined VEC2D
        void add(flt x, flt y, atom* a){add(fixedForceAtom(Vec(x,y),a));};
        #endif
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box *box){
            flt E=0;
            for(vector<fixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); it++)
                E += it->energy(box);
            return E;
        };
        void setForces(Box *box){
            for(vector<fixedForceAtom>::iterator it = atoms.begin(); it < atoms.end(); it++)
                it->setForce(box);
        };
        flt pressure(Box *box){return NAN;};
};


struct fixedSpringAtom {
    Vec loc;
    flt k;
    bool usecoord[3];
    atom *a;
    fixedSpringAtom(atom *a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true) : 
            loc(loc), k(k), a(a) {
                usecoord[0] = usex;
                usecoord[1] = usey;
                usecoord[2] = usez;
                };
    flt energy(Box *box){
            Vec diffx = a->x - loc;
            for(uint i=0; i<2; i++){
                if(!usecoord[i]) diffx[i] = 0;
            }
            return k*diffx.sq()/2;
        };
    void setForce(Box *box){
        Vec diffx = a->x - loc;
        for(uint i=0; i<2; i++){
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
        void add(atom *a, Vec loc, flt k, bool usex=true, bool usey=true, bool usez=true){
            add(fixedSpringAtom(a, loc, k, usex, usey, usez));};
        //void add(Vec F, atom* a){add(fixedForceAtom(F,a));};
        //void add(flt x, flt y, flt z, atom* a){add(fixedForceAtom(Vec(x,y,z),a));};
        uint size() const{ return (uint) atoms.size();};
        flt energy(Box *box){
            flt E=0;
            for(vector<fixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); it++)
                E += it->energy(box);
            return E;
        };
        void setForces(Box *box){
            for(vector<fixedSpringAtom>::iterator it = atoms.begin(); it < atoms.end(); it++)
                it->setForce(box);
        };
        flt pressure(Box *box){return NAN;};
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
        flt energy(Box *box){
            flt dx = (g1->com() - g2->com()).mag() - x0;
            return k/2 * dx * dx;
        };
        void setForces(Box *box){
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.mag();
            flt fmag = -k * (comdist - x0);
            Vec a1 = comvec * (fmag / m1 / comdist);
            for(uint i=0; i < g1->size(); i++){
                atom &atm = *(g1->get(i));
                atm.f += a1 * atm.m;
            }
            
            Vec a2 = comvec * (-fmag / m2 / comdist);
            for(uint i=0; i < g2->size(); i++){
                atom &atm = *(g2->get(i));
                atm.f += a2 * atm.m;
            }
        };
        flt pressure(Box *box){
            //~ return NAN;
            // I think this is right, but I haven't checked it
            Vec comvec = g1->com() - g2->com();
            flt comdist = comvec.mag();
            flt fmag = -k * (comdist - x0);
            
            Vec a12 = comvec * (fmag / m1 / m2 / comdist);
            
            flt P = 0;
            for(uint i=0; i < g1->size(); i++){
                for(uint j=0; j < g2->size(); j++){
                    atom &atm1 = *(g1->get(i));
                    atom &atm2 = *(g2->get(i));
                    Vec fij = a12 * (atm1.m * atm2.m);
                    Vec rij = box->diff(atm1.x, atm2.x);
                    P += fij.dot(rij);
                }
            }
            return P;
        };
};

struct bondgrouping {
    flt k, x0;
    atom *a1, *a2;
    bondgrouping(flt k, flt x0, atom* a1, atom* a2) : 
                k(k),x0(x0), a1(a1), a2(a2){};
};

class bondpairs : public interaction {
    protected:
        vector<bondgrouping> pairs;
        inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    public:
        bondpairs(vector<bondgrouping> pairs = vector<bondgrouping>());
        void add(bondgrouping b){pairs.push_back(b);};
        void add(flt k, flt x0, atom* a1, atom* a2){add(bondgrouping(k,x0,a1,a2));};
        uint size() const{ return (uint) pairs.size();};
        flt mean_dists() const;
        flt std_dists() const;
        flt energy(Box *box);
        void setForces(Box *box);
        flt pressure(Box *box);
};

struct anglegrouping {
    flt k, x0;
    atom *a1, *a2, *a3;
    anglegrouping(flt k, flt x0, atom* a1, atom* a2, atom *a3) : 
                k(k),x0(x0), a1(a1), a2(a2), a3(a3){};
};

class angletriples : public interaction {
    protected:
        vector<anglegrouping> triples;
        inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    public:
        angletriples(vector<anglegrouping> triples = vector<anglegrouping>());
        void add(anglegrouping b){triples.push_back(b);};
        void add(flt k, flt x0, atom* a1, atom* a2, atom* a3){
                                add(anglegrouping(k,x0,a1,a2,a3));};
        inline flt energy(Box *box);
        inline flt pressure(Box *box){return 0;};
        inline void setForces(Box *box);
        uint size() const {return (uint) triples.size();};
        flt mean_dists() const;
        flt std_dists() const;
};

#ifdef VEC3D
struct dihedralgrouping {
    inline static Vec diff(Vec r1, Vec r2){return r1-r2;};
    dihedral dih;
    atom *a1, *a2, *a3, *a4;
    dihedralgrouping(vector<flt> coscoeffs, vector<flt> sincoeffs,
                atom* a1, atom* a2, atom* a3, atom *a4, bool usepow=true) : 
                dih(coscoeffs, sincoeffs, usepow), a1(a1), 
                a2(a2), a3(a3), a4(a4){};
};

class dihedrals : public interaction {
    protected:
        vector<dihedralgrouping> groups;
    public:
        dihedrals(vector<dihedralgrouping> pairs = vector<dihedralgrouping>());
        void add(dihedralgrouping b){groups.push_back(b);};
        inline void add(vector<flt> nums, atom* a1, atom* a2, atom* a3, atom *a4){
            add(dihedralgrouping(nums, vector<flt>(), a1,a2,a3,a4));};
        inline void add(vector<flt> coscoeffs, vector<flt> sincoeffs, 
                            atom* a1, atom* a2, atom* a3, atom *a4, bool usepow=true){
            add(dihedralgrouping(coscoeffs, sincoeffs,a1,a2,a3,a4, usepow));};
        uint size() const{ return uint(groups.size());};
        flt mean_dists() const;
        //~ flt std_dists() const;
        flt energy(Box *box);
        void setForces(Box *box);
        inline flt pressure(Box *box){return 0;};
};
//~ class LJgroup : private vector<LJdata>, public atomgroup {
    //~ public:
        //~ LJgroup() : vector<LJdata>(){};
        //~ add(flt sigma, flt epsilon, atom *a){ 
//~ };
#endif


////////////////////////////////////////////////////////////////////////
struct forcepair {
    atom *a1, *a2;
    Vec fij;
};

struct forcepairx {
    atom *a1, *a2;
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
        virtual void setForces(Box *box, fpairxFunct*)=0;
        virtual ~interactionpairsx(){};
};

//~ struct atompaircomp {
    //~ bool operator() (const atompair& lhs, const atompair& rhs) const{
        //~ return (lhs[0] == rhs[0] and lhs[1] == rhs[1]);}
//~ };

class pairlist {
    protected:
        map<const atomid, set<atomid> > pairs;
    public:
        //~ pairlist(atomgroup *group);
        pairlist(){};
        
        inline void ensure(const atomid a){
            pairs.insert(std::pair<atomid, set<atomid> >(a, set<atomid>()));
        }
        inline void ensure(vector<atomid> ps){
            vector<atomid>::iterator it;
            for(it=ps.begin(); it != ps.end(); it++) ensure(*it);
        }
        inline void ensure(atomgroup &group){
            for(uint i=0; i<group.size(); i++) ensure(group.get_id(i));
        }
        
        inline bool has_pair(atomid a1, atomid a2){
            if(a1 > a2) return pairs[a1].count(a2) > 0;
            else return pairs[a2].count(a1) > 0;
        }
        
        inline void add_pair(atomid a1, atomid a2){
            //~ cout << "pairlist ignore " << a1.n() << '-' << a2.n() << "\n";
            if(a1 > a2){pairs[a1].insert(a2);}
            else{pairs[a2].insert(a1);};
        }
        
        inline void erase_pair(atomid a1, atomid a2){
            if(a1 > a2){pairs[a1].erase(a2);}
            else{pairs[a2].erase(a1);};
        }
        
        inline set<atomid> get_pairs(const atomid a){ensure(a); return pairs[a];};
        
        // for iterating over neighbors
        inline set<atomid>::iterator begin(const atomid a){return pairs[a].begin();};
        inline set<atomid>::iterator end(const atomid a){return pairs[a].end();};
        
        inline uint size() const { uint N=0; for(
            map<const atomid, set<atomid> >::const_iterator it=pairs.begin();
            it != pairs.end(); it++) N+= (uint) it->second.size();
            return N;
        };
        
        void clear();
};

class neighborlist : public statetracker{
    //maintains a Verlet list of "neighbors": molecules within a 
    // 'skin radius' of each other.
    // note that molecules are counted as neighbors if any point within
    // their molecular radius is within a 'skin radius' of any point
    // within another molecule's molecular radius.
    //
    // update(false) should be called frequently; it checks (O(N)) if
    // any two molecules might conceivably overlap by more than a critical
    // distance, and if so, it updates all the neighbor lists.
    
    // the <bool areneighbors()> function returns whether two molecules
    // are neighbors, and the begin(n), end(n) allow for iterating over
    // the neighbor lists
    protected:
        Box *box;
        flt critdist, skinradius;
        metagroup atoms;
        vector<idpair> curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        atomid get_id(atom* a);
        bool ignorechanged; // if true, forces a full check on next update
        //~ bool checkneighbors(const uint n, const uint m) const;
        // this is a full check
    public:
        neighborlist(Box *box, const flt innerradius, const flt outerradius);
        neighborlist(Box *box, atomgroup &vec, const flt innerradius, 
        const flt outerradius, pairlist ignore = pairlist());
        void update(Box *newbox){assert(newbox == box); update_list(false);};
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first
        
        inline uint which(){return updatenum;};
        inline uint numpairs(){return (uint) curpairs.size();};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b); ignorechanged=true;};
        void ignore(atom*, atom*);
        atomid add(atom* a){
            atomid id = atoms.get_id(a);
            if(id != NULL) return id;
            atoms.add(a);
            assert(lastlocs.size() == atoms.size() - 1);
            lastlocs.push_back(a->x);
            id = atoms.get_id(a);
            ignorechanged = true;
            return id;
        }
        
        inline void changesize(flt inner, flt outer){
            critdist = inner; skinradius = outer; update_list(true);};
        inline void changesize(flt ratio){
            skinradius = critdist*ratio; update_list(true);};

        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint size() const{return atoms.size();};
        inline vector<idpair>::iterator begin(){return curpairs.begin();};
        inline vector<idpair>::iterator end(){return curpairs.end();};
        inline idpair get(uint i){return curpairs[i];};
        //~ inline vector<idpair> getpairs(){return vector<idpair>(curpairs);};
        ~neighborlist(){};
};

neighborlist* neighborlistL(Box *box, const double innerradius, const double outerradius){
    return new neighborlist(box, innerradius, outerradius);
};

class ContactTracker : public statetracker{
    protected:
        atomgroup *atoms;
        vector<flt> dists;
        vector<vector<bool> > contacts;
        
        unsigned long long breaks;
        unsigned long long formations;
        unsigned long long incontact;
    public:
        ContactTracker(Box *box, atomgroup *atoms, vector<flt> dists);
        void update(Box *box);
        
        unsigned long long broken(){return breaks;};
        unsigned long long formed(){return formations;};
        unsigned long long number(){return incontact;};
};

ContactTracker* ContactTrackerD(Box *box, atomgroup *atoms, vector<double> dists){
    vector<flt> newdists = vector<flt>();
    for(uint i=0; i<dists.size(); i++){
        newdists.push_back(dists[i]);
    }
    return new ContactTracker(box, atoms, newdists);
}

class EnergyTracker : public statetracker{
    protected:
        atomgroup *atoms;
        vector<interaction*> interactions;
        
        uint N;
        uint nskip, nskipped;
        flt U0;
        flt Es, Us, Ks;
        flt Esq, Usq, Ksq;
    public:
        EnergyTracker(atomgroup *atoms, 
            vector<interaction*> interactions, uint nskip=1)
             : atoms(atoms),
            interactions(interactions), N(0), nskip(max(nskip,1u)), nskipped(0),
            U0(0),Es(0),Us(0),Ks(0), Esq(0), Usq(0), Ksq(0){};
        void update(Box *box);
        void reset(){
            nskipped=0;
            N=0; Es=0; Us=0; Ks=0;
            Esq=0; Usq=0; Ksq=0;
        };
        void setU0(flt newU0){
            U0 = newU0;
            reset();
        };
        void setU0(Box *box);
        flt getU0(){return U0;};
            
        flt E(){return Es/((flt) N);};
        flt U(){return Us/((flt) N);};
        flt K(){return Ks/((flt) N);};
        flt Estd(){return sqrt(Esq/N -Es*Es/N/N);};
        flt Kstd(){return sqrt(Ksq/N -Ks*Ks/N/N);};
        flt Ustd(){return sqrt(Usq/N -Us*Us/N/N);};
        flt Esqmean(){return Esq/N;};
        flt Ksqmean(){return Ksq/N;};
        flt Usqmean(){return Usq/N;};
        //~ flt Ustd(){return sqrt((Usq -(U*U)) / ((flt) N));};
        //~ flt Kstd(){return sqrt((Ksq -(K*K)) / ((flt) N));};
        uint n(){return N;};
};

template <class A, class P>
class SimpleListed : public interaction {
    protected:
        vector<A> atoms;
    
    public:
        SimpleListed(){};
        inline void add(A atm){atoms.push_back(atm);};
        //flt energy(Box *box, idpair &pair);
        flt energy(Box *box);
        flt pressure(Box *box);
        uint size(){return ((uint) (atoms.size()));};
        //inline flt energy_pair(P pair, Box *box){return pair.energy(box);}; // This may need to be written!
        void setForces(Box *box);
        flt setForcesGetPressure(Box *box);
        //inline Vec forces_pair(P pair, Box *box){return pair.forces(box);}; // This may need to be written!
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
        vector<A> atoms;
        vector<P> pairs;
        neighborlist *neighbors;
        uint lastupdate;
    public:
        NListed(neighborlist *neighbors) : neighbors(neighbors){};
        inline void add(A atm){
            atomid id = neighbors->add(atm.pointer());
            A a = A(id, atm);
            assert(a.n() <= atoms.size());
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            atoms[a.n()] = a;};
        void update_pairs();
        P getpair(idpair &pair){
            return P(atoms[pair.first().n()], atoms[pair.last().n()]);}
        A& getatom(uint n){return atoms[n];}
        flt energy(Box *box, idpair &pair);
        flt energy(Box *box);
        flt pressure(Box *box);
        uint size(){return ((uint) (atoms.size()));};
        inline flt energy_pair(P pair, Box *box){return pair.energy(box);}; // This may need to be written!
        void setForces(Box *box);
        flt setForcesGetPressure(Box *box);
        inline Vec forces_pair(P pair, Box *box){return pair.forces(box);}; // This may need to be written!
        inline vector<A> &atom_list(){return atoms;};
        neighborlist *nlist(){return neighbors;};
        //~ flt energy_test(flt dist);
        //~ flt force_test(flt dist);
        ~NListed(){};
};

template <class A, class P>
class NListedVirial : public interactionpairsx, public NListed<A,P> {
    public:
        NListedVirial(neighborlist *neighbors) : NListed<A,P>(neighbors){};
        void setForces(Box *box){NListed<A,P>::setForces(box);};
        void setForces(Box *box, fpairxFunct*);
        virtual inline flt setForcesGetPressure(Box *box){return NListed<A,P>::setForcesGetPressure(box);};
        virtual flt setForcesGetEnergy(Box *box);
        virtual inline flt energy(Box *box){return NListed<A,P>::energy(box);};
        virtual inline flt pressure(Box *box){return NListed<A,P>::pressure(box);};
};

struct Charged : public atomid {
    flt q;
    Charged() : atomid(), q(0){};
    Charged(flt q, atom *a) : atomid(a), q(q){};
};

struct ChargePair {
    flt q1q2;
    atomid atom1, atom2;
    ChargePair(Charged a1, Charged a2) : q1q2(a1.q*a2.q){};
};

////////////////////////////////////////////////////////////////////////
// Repulsive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// cutoff at sigma
struct LJatom : public atomid {
    flt epsilon, sigma;
    LJatom(flt epsilon, flt sigma, atom* a) : atomid(a), 
            epsilon(epsilon), sigma(sigma){};
    LJatom(atomid a, LJatom other) : atomid(a),
                    epsilon(other.epsilon), sigma(other.sigma){};
};

struct LJpair {
    flt epsilon, sigma;
    atomid atom1, atom2;
    LJpair(LJatom LJ1, LJatom LJ2) :
            epsilon(sqrtflt(LJ1.epsilon * LJ2.epsilon)),
            sigma((LJ1.sigma + LJ2.sigma) / 2),
            atom1(LJ1), atom2(LJ2){};
    inline flt energy(Box *box){
        return LJrepulsive::energy(box->diff(atom1.x(),atom2.x()), epsilon, sigma);};
    inline Vec forces(Box *box){
        return LJrepulsive::forces(box->diff(atom1.x(),atom2.x()), epsilon, sigma);};
};

////////////////////////////////////////////////////////////////////////
// Purely attractive LJ, with ε = √(ε₁ ε₂) and σ = (σ₁+σ₂)/2
// cutoff at some sigcut
struct LJatomcut : public LJatom {
    flt sigcut; // sigma units
    LJatomcut(flt epsilon, flt sigma, atom* a, flt cut) : 
            LJatom(epsilon, sigma, a), sigcut(cut){};
    LJatomcut(atomid a, LJatomcut other) : LJatom(a, other), 
        sigcut(other.sigcut){};
};

struct LJAttractPair {
    LJattractCut inter;
    atomid atom1, atom2;
    LJAttractPair(LJatomcut a1, LJatomcut a2) : 
        inter(sqrtflt(a1.epsilon * a2.epsilon),
              (a1.sigma + a2.sigma) / 2, 
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box *box){return inter.energy(box->diff(atom1.x(), atom2.x()));};
    inline Vec forces(Box *box){return inter.forces(box->diff(atom1.x(), atom2.x()));};
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
    HydroAtom(vector<flt> epsilons, uint indx, flt sigma, atom* a, flt cut) : 
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
};

struct HydroPair {
    LJattractCut inter;
    atomid atom1, atom2;
    HydroPair(HydroAtom a1, HydroAtom a2) : 
        inter(a1.getEpsilon(a2),
              (a1.sigma + a2.sigma) / 2, 
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box *box){return inter.energy(box->diff(atom1.x(), atom2.x()));};
    inline Vec forces(Box *box){return inter.forces(box->diff(atom1.x(), atom2.x()));};
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
    LJAtomIndexed(vector<flt> epsilons, vector<flt> sigmas, uint indx, atom *a, flt cut) : 
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
};

struct LJFullPair {
    LJattractCut inter;
    atomid atom1, atom2;
    LJFullPair(LJAtomIndexed a1, LJAtomIndexed a2) : 
        inter(a1.getEpsilon(a2),
              a1.getSigma(a2),
              max(a1.sigcut, a2.sigcut)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box *box){return inter.energy(box->diff(atom1.x(), atom2.x()));};
    inline Vec forces(Box *box){return inter.forces(box->diff(atom1.x(), atom2.x()));};
};


////////////////////////////////////////////////////////////////////////
// Full LJish, with ε = ε₁₂ and σ = σ₁₂ (both indexed)
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
    LJishAtom(atom *a, vector<flt> epsilons, flt repeps, flt sigma,
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
        return sqrtflt(sigma * other.sigma);
    }
};

struct LJishPair {
    flt epsilon, repeps, sigma, n, cutR, cutE;
    atomid atom1, atom2;
    LJishPair(LJishAtom LJ1, LJishAtom LJ2) :
            epsilon(LJ1.getEpsilon(LJ2)),
            repeps(sqrtflt(LJ1.repeps * LJ2.repeps)),
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
        flt mid = (1-powflt(cutR,-n));
        cutE = epsilon*(mid*mid);
    };
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt rsq = rij.sq()/(sigma*sigma);
        if(rsq > cutR * cutR){
            //~ printf("LJish: dx=%.2f, σ=%.2f, rsq=%.2f, cutR²=%.2f\n", rij.mag(), sigma, rsq, cutR);
            return 0;
        }
        
        flt mid = (1-powflt(rsq,-n/2));
        //~ printf("LJish: rsq=%.2f, mid %.2f, epsilon %.2f, repeps %.2f\n", rsq, mid, epsilon, repeps);
        if (rsq > 1) return epsilon*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
    };
    inline Vec forces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        flt rsq = dsq/(sigma*sigma);
        if(rsq > cutR * cutR) return Vec();
        flt rmid = powflt(rsq,-n/2); // σ^n/r^n
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
    LJAttractRepulseAtom(atom* a, vector<flt> epsilons, flt sigma, uint indx, flt cut) : 
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
            flt mid = (1-powflt(cutR,-6));
            cutE = eps*(mid*mid);
        };
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt rsq = rij.sq()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrtflt(rij.sq()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-powflt(rsq,-3));
        //~ if(eps*(mid*mid) - cutE < 0) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                //~ rij.mag(), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        //~ }
        return eps*(mid*mid) - cutE;
    };
    inline Vec forces(Box *box){
        if(eps == 0) return Vec();
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec();
        flt rsix = powflt(rsq,-3); // σ⁶/r⁶
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
    LJAttractFixedRepulseAtom(atom *a, vector<flt> epsilons, flt repeps, 
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
        eps(abs(a1.getEpsilon(a2))), repeps(sqrtflt(a1.repeps * a2.repeps)), 
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
            flt mid = (1-powflt(cutR,-6));
            cutE = eps*(mid*mid);
        };
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt rsq = rij.sq()/(sig*sig);
        if(rsq > cutR*cutR) {
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrtflt(rij.sq()), 0.0, eps, sig, cutR, cutE);
            return 0;
        }
        flt mid = (1-powflt(rsq,-3)); // # 1 - σ⁶/r⁶
        //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrtflt(rij.sq()), eps*(mid*mid) - cutE, eps, sig, cutR, cutE);
        if (rsq > 1) return eps*(mid*mid) - cutE;
        return repeps*(mid*mid) - cutE;
        //~ flt E;
        //~ if (rsq > 1) E = eps*(mid*mid) - cutE;
        //~ else E = repeps*(mid*mid) - cutE;
        //~ if(E > 1e4)
            //~ printf("Distance: %.2f Energy: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrtflt(rij.sq()), E, eps, repeps, sig, cutR, cutE);
        //~ return E;
    };
    inline Vec forces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        flt rsq = dsq/(sig*sig);
        if(rsq > (cutR*cutR)) return Vec();
        flt rsix = powflt(rsq,-3); // σ⁶/r⁶
        flt fmagTimesR = 12*rsix*(rsix - 1);
        if (rsq < 1) return rij * (repeps * fmagTimesR / dsq);
        return rij * (eps * fmagTimesR / dsq);
        //~ flt fmag;
        //~ if (rsq < 1) fmag = repeps * fmagTimesR / dsq;
        //~ else fmag = eps * fmagTimesR / dsq;
        //~ if(fmag * rij.mag() > 1e4)
            //~ printf("Distance: %.2f Force: %.2f (ε: %.2f,%.2f σ: %.2f cut: %.2f cutE: %.2f)\n", 
                    //~ sqrtflt(rij.sq()), fmag * rij.mag(), eps, repeps, sig, cutR, cutE);
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
    LJDoubleAtom(flt epsilon, flt epsrep, flt sigma, atom *a, flt cut) : 
            LJatom(epsilon, sigma, a), epsrep(epsrep), sigcut(cut){};
    LJDoubleAtom(atomid a, LJDoubleAtom other) : LJatom(a, other), 
        epsrep(other.epsrep), sigcut(other.sigcut){};
};

struct LJDoublePair : public LJAttractFixedRepulsePair {
    flt eps, repeps, sig;
    flt cutR, cutE;
    bool attract;
    atomid atom1, atom2;
    LJDoublePair(LJDoubleAtom a1, LJDoubleAtom a2) : 
        eps(sqrtflt(a1.epsilon * a2.epsilon)), repeps(sqrtflt(a1.epsrep * a2.epsrep)), 
        sig((a1.sigma + a2.sigma)/2.0),
        cutR(max(a1.sigcut, a2.sigcut)), attract(a1.epsilon * a2.epsilon > 0),
        atom1(a1), atom2(a2){
            if(!attract or eps == 0){
                cutR = 1;
                cutE = 0;
                eps = 0;
                attract = false;
                return;
            }
            flt mid = (1-powflt(cutR,-6));
            cutE = eps*(mid*mid);
        };
};

struct EisMclachlanAtom : public atomid {
    flt dist, sigmai; // dist is radius of atom + radius of water (1.4 Å)
    EisMclachlanAtom(flt dist, flt sigmai, atom *a) : atomid(a),
            dist(dist), sigmai(sigmai){};
    EisMclachlanAtom(atomid a, EisMclachlanAtom other) : atomid(a), 
        dist(other.dist), sigmai(other.sigmai){};
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
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt R = rij.mag();
        if (R > cutoff) return 0;
        return c0/R + c1 + c2*R;
    }
    inline Vec forces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > (cutoff*cutoff)) return Vec();
        flt R = sqrtflt(dsq);
        return rij * ((c0/dsq-c2)/R);
    }
};

////////////////////////////////////////////////////////////////////////
// Hertzian potential, with ε = √(ε₁ ε₂) and σ = (σ₁ + σ₂)/2
// Potential is V(r) = ε/n (1 - r/σ)^n, with n = 5/2 usually
// cutoff at r = σ

struct HertzianAtom : public atomid {
    flt eps, sigma, exponent;
    HertzianAtom(atom *a, flt eps, flt sigma, flt exponent=2.5) : atomid(a),
            eps(eps), sigma(sigma), exponent(exponent){};
    HertzianAtom(atomid a, HertzianAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), exponent(other.exponent){};
};

HertzianAtom hertzd(atom *a, double eps, double sigma, double exponent=2.5){
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
        eps(sqrtflt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        exponent((a1.exponent + a2.exponent)/2.0), atom1(a1), atom2(a2){};
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sig*sig) return 0.0;
        flt R = sqrtflt(dsq);
        return eps * powflt(1.0 - (R/sig), exponent) / exponent;
    }
    inline Vec forces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sig*sig) return Vec();
        flt R = sqrtflt(dsq);
        return rij * (eps * powflt(1.0 - (R/sig), exponent-1) /sig/R);
    }
    inline EnergyForce EnergyForces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sig*sig) return EnergyForce(Vec(),0);
        flt R = sqrtflt(dsq);
        
        Vec f = rij * (eps * powflt(1.0 - (R/sig), exponent-1) /sig/R);
        flt E = eps * powflt(1.0 - (R/sig), exponent) / exponent;
        return EnergyForce(f, E);
    }
    //~ inline flt xrij(Box *box){
        //~ Vec rij = box->diff(atom1.x(), atom2.x());
        //~ flt dsq = rij.sq();
        //~ if(dsq > sig*sig) return 0.0;
        //~ flt R = sqrtflt(dsq);
        //~ return (R*eps*(exponent-1)/sig/sig) * powflt(1.0 - (R/sig), exponent-2);
    //~ }
    inline void fill(Box *box, forcepairx &fpair){
        fpair.a1 = atom1.pointer();
        fpair.a2 = atom2.pointer();
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sig*sig){
            fpair.fij = Vec();
            fpair.xij = 0.0;
            return;
        }
        flt R = sqrtflt(dsq);
        fpair.fij = rij * (eps*powflt(1.0 - (R/sig), exponent-1) /sig/R);
        flt Rs = R / sig;
        fpair.xij = -Rs*eps*powflt(1.0 - (R/sig), exponent-2)*(1-exponent*Rs);
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
 
*/

struct LoisOhernAtom : public atomid {
    flt eps, sigma, C, l;
    LoisOhernAtom(atom *a, flt eps, flt sigma, flt C, flt l) : atomid(a),
            eps(eps), sigma(sigma), C(C), l(l){};
    LoisOhernAtom(atomid a, LoisOhernAtom other) : atomid(a), 
        eps(other.eps), sigma(other.sigma), C(other.C), l(other.l){};
};

struct LoisOhernPair {
    flt eps, sig, C, l, sigcut;
    atomid atom1, atom2;
    LoisOhernPair(LoisOhernAtom a1, LoisOhernAtom a2) : 
        eps(sqrtflt(a1.eps * a2.eps)), sig((a1.sigma + a2.sigma)/2.0),
        C((a1.C + a2.C)/2.0), l((a1.l + a2.l)/2.0), sigcut(sig*(1+C+l)),
        atom1(a1), atom2(a2){};
    inline flt energy(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sigcut*sigcut) return 0.0;
        flt R = sqrtflt(dsq)/sig;
        if(R < 1 + C){
            flt dR = R-1;
            return -eps*sig/2*(C*(C+l) - dR*dR);
        }
        
        flt dR2 = C+l+1 - R;
        return -C*eps*sig/2/l*dR2*dR2;
    }
    
    inline Vec forces(Box *box){
        Vec rij = box->diff(atom1.x(), atom2.x());
        flt dsq = rij.sq();
        if(dsq > sigcut*sigcut) return Vec();
        flt R = sqrtflt(dsq);
        flt rsig = R/sig;
        
        if(rsig < 1 + C){
            flt dR = rsig-1;
            return rij*(-eps*dR/R);
        }
        
        flt dR2 = rsig-(C+l+1);
        return rij*(C*eps/l*dR2/R);
    }
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
        inline void add(atom *a, flt epsilon, flt sigma){
            add(LJatom(epsilon,sigma,a));};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        inline void ignore(atom* a, atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint atoms_size() const{return (uint) atoms.size();};
        flt energy(Box *box);
        flt pressure(Box *box);
        void setForces(Box *box);
        //~ ~LJsimple(){};
};

template <class A, class P>
flt SimpleListed<A, P>::energy(Box *box){
    flt E = 0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); it1++){
        for(it2 = it1+1; it2 != atoms.end(); it2++){
            E += P(*it1, *it2).energy(box);
        }
    }
    return E;
};

template <class A, class P>
void SimpleListed<A, P>::setForces(Box *box){
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); it1++){
        for(it2 = it1+1; it2 != atoms.end(); it2++){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            it1->f() += f;
            it2->f() -= f;
        }
    }
};

template <class A, class P>
flt SimpleListed<A, P>::setForcesGetPressure(Box *box){
    flt p=0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); it1++){
        for(it2 = it1+1; it2 != atoms.end(); it2++){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            it1->f() += f;
            it2->f() -= f;
            Vec r = box->diff(it1->x(), it2->x());
            p += r.dot(f);
        }
    }
    //~ cout << "Set forces, got pressure" << p << '\n';
    return p;
};

template <class A, class P>
flt SimpleListed<A, P>::pressure(Box *box){
    flt p=0;
    typename vector<A>::iterator it1;
    typename vector<A>::iterator it2;
    for(it1 = atoms.begin(); it1 != atoms.end(); it1++){
        for(it2 = it1+1; it2 != atoms.end(); it2++){
            P pair = P(*it1, *it2);
            Vec f = pair.forces(box);
            Vec r = box->diff(it1->x(), it2->x());
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
    for(pairit = neighbors->begin(); pairit != neighbors->end(); pairit++){
        A firstatom = atoms[pairit->first().n()];
        A secondatom = atoms[pairit->last().n()];
        pairs.push_back(P(firstatom, secondatom));
    }
}

template <class A, class P>
flt NListed<A, P>::energy(Box *box, idpair &pair){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    P Epair = P(atoms[pair.first().n()],atoms[pair.last().n()]);
    //~ Vec dist = box->diff(Epair.atom1.x(), Epair.atom2.x());
    return energy_pair(Epair, box);
};

template <class A, class P>
flt NListed<A, P>::energy(Box *box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt E = 0;
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        //~ Vec dist = box->diff(it->atom1.x(), it->atom2.x());
        E += energy_pair(*it, box);
    }
    return E;
};

template <class A, class P>
void NListed<A, P>::setForces(Box *box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        Vec f = forces_pair(*it, box);
        it->atom1.f() += f;
        it->atom2.f() -= f;
        //~ assert(f.sq() < 1000000);
    }
};

template <class A, class P>
flt NListedVirial<A, P>::setForcesGetEnergy(Box *box){
    NListed<A, P>::update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt E = 0;
    typename vector<P>::iterator it;
    for(it = NListed<A, P>::pairs.begin(); it != NListed<A, P>::pairs.end(); it++){
        EnergyForce EF = it->EnergyForces(box);
        it->atom1.f() += EF.f;
        it->atom2.f() -= EF.f;
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
void NListedVirial<A, P>::setForces(Box *box, fpairxFunct* funct){
    NListed<A, P>::update_pairs(); // make sure the LJpairs match the neighbor list ones
    forcepairx myfpair;
    typename vector<P>::iterator it;
    for(it = NListed<A, P>::pairs.begin(); it != NListed<A, P>::pairs.end(); it++){
        it->fill(box, myfpair);
        //~ assert(!isnan(myfpair.fij[0]));
        funct->run(&myfpair);
        it->atom1.f() += myfpair.fij;
        it->atom2.f() -= myfpair.fij;
        //~ assert(f.sq() < 1000000);
    }
};

template <class A, class P>
flt NListed<A, P>::setForcesGetPressure(Box *box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt p=0;
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        Vec f = forces_pair(*it, box);
        it->atom1.f() += f;
        it->atom2.f() -= f;
        Vec r = box->diff(it->atom1.x(), it->atom2.x());
        p += r.dot(f);
    }
    //~ cout << "Set forces, got pressure" << p << '\n';
    return p;
};

template <class A, class P>
flt NListed<A, P>::pressure(Box *box){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt p=0;
    //~ printf("Updating pressure...\n");
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        Vec f = forces_pair(*it, box);
        Vec r = box->diff(it->atom1.x(), it->atom2.x());
        //Vec r = it->atom1.x() - it->atom2.x();
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
        inline void add(atom *a, flt q){add(Charged(q,a));};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        inline void ignore(atom* a, atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint size() const{return (uint) atoms.size();};
        flt energy(Box *box);
        flt pressure(Box *box);
        void setForces(Box *box);
        //~ ~LJsimple(){};
};


bool toBuffer(vector<Vec*> arr, double* buffer, size_t sizet) {
    if(sizet < NDIM * arr.size()){return false;};
    for(uint i=0; i < arr.size(); i++)
    for(uint j=0; j < NDIM; j++)
    {
        buffer[i*NDIM + j] = (double) ((*arr[i])[j]);
    }
    return true;
};

////////////////////////////////////////////////////////////////////////
// For comparing two jammed structures

/* We have two packings, A and B, and want to know the sequence {A1, A2, A3...}
 * such that particle A1 of packing 1 matches particle 1 of packing B.
 * A jamminglist is a partial list; it has a list {A1 .. An}, with n / N
 * particles assigned, with a total distance² of distsq.
*/ 
class jamminglist {
    public:
        vector<uint> assigned;
        flt distsq;
        
        jamminglist() : assigned(), distsq(0){};
        jamminglist(const jamminglist& other) 
            : assigned(other.assigned), distsq(other.distsq){};
        jamminglist(const jamminglist& other, uint expand, flt addeddist)
            : assigned(other.size() + 1, 0), distsq(other.distsq + addeddist){
            for(uint i=0; i < other.size(); i++){
                assigned[i] = other.assigned[i];
            }
            assigned[assigned.size()-1] = expand;
        }
        inline uint size() const {return (uint) assigned.size();};
        
        bool operator<(const jamminglist& other);
};

class jammingtree {
    private:
        Box *box;
        list<jamminglist> jlists;
        vector<Vec> A;
        vector<Vec> B;
    public:
        jammingtree(Box *box, vector<Vec>& A, vector<Vec>& B)
            : box(box), jlists(), A(A), B(B) {
            jlists.push_back(jamminglist());
            assert(A.size() <= B.size());
        };

        bool expand(){
            jamminglist curjlist = jlists.front();
            vector<uint>& curlist = curjlist.assigned;
            if(curlist.size() >= A.size()){
                //~ cout << "List already too big\n";
                return false;
            }
            
            list<jamminglist> newlists = list<jamminglist>();
            for(uint i=0; i < B.size(); i++){
                vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
                //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
                if (found != curlist.end()){
                    //~ cout << "Found " << i << "\n";
                    //cout << found << '\n';
                    continue;
                }
                flt newdist = box->diff(A[curlist.size()], B[i]).sq();
                jamminglist newjlist = jamminglist(curjlist, i, newdist);
                newlists.push_back(newjlist);
                //~ cout << "Made " << i << "\n";
            }
            
            if(newlists.size() <= 0){
                //~ cout << "No lists made\n";
                return false;
            }
            //~ cout << "Have " << newlists.size() << "\n";
            newlists.sort();
            //~ cout << "Sorted.\n";
            jlists.pop_front();
            //~ cout << "Popped.\n";
            jlists.merge(newlists);
            //~ cout << "Merged to size " << jlists.size() << "best dist now " << jlists.front().distsq << "\n";
            return true;
        }
        bool expand(uint n){
            bool retval=false;
            for(uint i=0; i<n; i++){
                retval = expand();
            }
            return retval;
        }
        list<jamminglist> &mylist(){return jlists;};
        list<jamminglist> copylist(){return jlists;};
        
        jamminglist curbest(){
            jamminglist j = jamminglist(jlists.front());
            //~ cout << "Best size: " << j.size() << " dist: " << j.distsq;
            //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " << j.assigned[j.size()-1] << "]";
            //~ cout << '\n';
            return j;
            //return jamminglist(jlists.front());
            };
        uint size(){return (uint) jlists.size();};
};

#ifdef VEC2D

class jamminglistrot : public jamminglist {
    public:
        uint rotation;
        
        jamminglistrot() : jamminglist(), rotation(0){};
        jamminglistrot(uint rot) : jamminglist(), rotation(rot){};
        jamminglistrot(const jamminglistrot& other) 
            : jamminglist(other), rotation(other.rotation){};
        jamminglistrot(const jamminglistrot& other, uint expand, flt addeddist)
            : jamminglist(other, expand, addeddist), rotation(other.rotation){};
        
        bool operator<(const jamminglistrot& other);
};

// Includes rotations, flips, and translations.
class jammingtree2 {
    protected:
        Box *box;
        list<jamminglistrot> jlists;
        vector<Vec> A;
        vector<vector<Vec> > Bs;
    public:
        // make all 8 possible rotations / flips
        // then subtract off all possible COMVs
        jammingtree2(Box *box, vector<Vec>& A, vector<Vec>& B);
        flt distance(jamminglistrot& jlist);
        list<jamminglistrot> expand(jamminglistrot curjlist);
        
        virtual bool expand();
        
        bool expand(uint n){
            bool retval=false;
            for(uint i=0; i<n; i++){
                retval = expand();
                if(!retval) break;
            }
            return retval;
        }
        bool expandto(flt maxdistsq){
            bool retval = true;
            while((maxdistsq <= 0 or jlists.front().distsq < maxdistsq) and retval){
                retval = expand();
            };
            return retval;
        }
        static Vec straight_diff(Box *bx, vector<Vec>& A, vector<Vec>& B);
        static flt straight_distsq(Box *bx, vector<Vec>& A, vector<Vec>& B);
        
        list<jamminglistrot> &mylist(){return jlists;};
        list<jamminglistrot> copylist(){return jlists;};
        list<jamminglistrot> copylist(uint n){
            list<jamminglistrot>::iterator last = jlists.begin();
            advance(last, n);
            return list<jamminglistrot>(jlists.begin(), last);
        };
        
        
        jamminglistrot curbest(){
            if(jlists.size() <= 0){
                jamminglistrot bad_list = jamminglistrot();
                bad_list.distsq = -1;
                return bad_list;
                }
            jamminglistrot j = jamminglistrot(jlists.front());
            //~ cout << "Best size: " << j.size() << " dist: " << j.distsq;
            //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " << j.assigned[j.size()-1] << "]";
            //~ cout << '\n';
            return j;
            //return jamminglist(jlists.front());
            };
        
        //jamminglistrot operator[](uint i){
        //    assert(i < jlists.size());
        //    return jamminglistrot(jlists[i]);
        //};
        
        uint size(){return (uint) jlists.size();};
        
        vector<Vec> locationsB(jamminglistrot jlist);
        vector<Vec> locationsB(){return locationsB(curbest());};
        vector<Vec> locationsA(jamminglistrot jlist);
        vector<Vec> locationsA(){return locationsA(curbest());};
        virtual ~jammingtree2(){};
};


class jammingtreeBD : public jammingtree2 {
    /* For a bi-disperse packing.
     * 'cutoff' is the number of particles of the first kind; i.e., the
     * A vector should have A[0]..A[cutoff-1] be of particle type 1,
     * and A[cutoff]..A[N-1] of particle type 2.
     * This does much the same as jammingtree2, but doesn't check any 
     * reordering in which particles of one type are relabeled as another.
     * For exampe, with 2+2 particles (cutoff 2), we check
     * [0123],[1023],[0132],[1032]
     * But not
     * [0213],[0231],[0312],[0321],[1203],[1230],[1302],[1320],...
     * This means at most (cutoff! (N-cutoff)!) combinations are checked,
     * and not all N!, which can save a lot of time (as well as
     *  rejecting false combinations).
     */
    protected:
        uint cutoff1,cutoff2;
    public:
        jammingtreeBD(Box *box, vector<Vec>& A, vector<Vec>& B, uint cutoff) :
            jammingtree2(box, A, B), cutoff1(cutoff), cutoff2(cutoff){};
        jammingtreeBD(Box *box, vector<Vec>& A, vector<Vec>& B, 
                    uint cutoffA, uint cutoffB);// :
            //jammingtree2(box, A, B), cutoff1(cutoffA), cutoff2(cutoffB){};
        
        list<jamminglistrot> expand(jamminglistrot curjlist);
        bool expand();
        bool expand(uint n){return jammingtree2::expand(n);};
};
#endif

flt confineRange(flt minimum, flt val, flt maximum){
    if(val <= minimum) return minimum;
    if(val >= maximum) return maximum;
    return val;
}

class atompair : public array<atom, 2> {
    public:
        atompair() : array<atom, 2>(){};
        atompair(flt m) : array<atom, 2>(){
            vals[0].m = m/2;
            vals[1].m = m/2;
        };
        atompair(flt m1, flt m2) : array<atom, 2>(){
            vals[0].m = m1;
            vals[1].m = m2;
        };
        atompair(atom* a, atom* b){ vals[0] = *a; vals[1] = *b;};
        inline atom& first() {return vals[0];};
        inline atom& last() {return vals[1];};
};

class SCatomvec : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms, which are
    // arranged in pairs.
    private:
        atompair* atoms;
        uint sz;
    public:
        SCatomvec(vector<double> masses) : sz((uint) masses.size()){
            atoms = new atompair[sz];
            for(uint i=0; i < sz; i++){
                atoms[i].first().m = masses[i]/2;
                atoms[i].last().m = masses[i]/2;
            }
        };
        SCatomvec(uint N, flt mass) : sz(N){
            atoms = new atompair[sz];
            for(uint i=0; i < sz; i++){
                atoms[i].first().m = mass/2;
                atoms[i].last().m = mass/2;
            }
        };
        SCatomvec(SCatomvec& other) : sz(other.pairs()){
            atoms = new atompair[sz];
            for(uint i=0; i < sz; i++) atoms[i] = other.atoms[i];
        };
        inline atom& operator[](cuint n){return atoms[n / 2][n % 2];};
        inline atom& operator[](cuint n) const {return atoms[n / 2][n % 2];};
        inline atompair& pair(cuint n){return atoms[n];};
        //inline atompair& pair(cuint n) const {return atoms[n / 2];};
        //~ atomid get_id(atom *a){
            //~ uint n = (uint) (a - atoms); WON'T WORK
            //~ if (n >= sz or a < atoms) return atomid();
            //~ return get_id(n);
        //~ }
        //~ inline atomid get_id(uint n) {
            //~ if (n > sz*2) return atomid(); return atomid(&(atoms[n/2][n%2]),n);};
        inline uint size() const {return sz*2;};
        inline uint pairs() const {return sz;};
        ~SCatomvec(){ delete [] atoms;};
};

struct SpheroCylinderDiff{
    Vec delta, r;
    flt lambda1, lambda2;
};

struct SCPair {
    atompair &p1;
    atompair &p2;
    flt l1, l2;
    SCPair(atompair &p1, atompair &p2, flt l1, flt l2) : 
        p1(p1), p2(p2), l1(l1), l2(l2){};
    SCPair(atompair &p1, atompair &p2, flt l) : 
        p1(p1), p2(p2), l1(l), l2(l){};
    SpheroCylinderDiff NearestLoc(Box *box);
    void applyForce(Box *box, Vec f, SpheroCylinderDiff diff, flt I);
};

struct SCSpringPair : public SCPair {
    flt eps, sig;
    
    SCSpringPair(atompair &p1, atompair &p2, flt eps, flt sig, flt l1, flt l2) : 
        SCPair(p1, p2, l1, l2), eps(eps), sig(sig){};
    SCSpringPair(atompair &p1, atompair &p2, flt eps, flt sig, flt l) : 
        SCPair(p1, p2, l), eps(eps), sig(sig){};
    
    inline flt maxdist(){return sig + (l1+l2)/2;};
    inline flt maxdelta(){return sig;};
    
    flt energy(Box *box, SpheroCylinderDiff diff){
        //~ atom &a1 = p1.first();
        //~ atom &a1p = p1.last();
        //~ atom &a2 = p2.first();
        //~ atom &a2p = p2.last();
        //~ Vec r1 = (a1.x + a1p.x)/2, r2 = (a2.x + a2p.x)/2;
        //~ Vec r12 = r2 - r1;
        //~ flt rsq = r12.sq();
        //~ if(rsq > pow(sig+l, 2)) return 0.0;
        
        //~ SpheroCylinderDiff diff = SCNearestLoc(a1.x, a1p.x, a2.x, a2p.x);
        flt dsq = diff.delta.sq();
        if(dsq > sig*sig) return 0;
        flt d = sqrtflt(dsq);
        flt dsig = d-sig;
        return dsig*dsig*eps/2;
    }
    Vec forces(Box *box, SpheroCylinderDiff diff){
        flt dsq = diff.delta.sq();
        if(dsq > sig*sig) return Vec();
        flt dmag = sqrtflt(dsq);
        Vec dhat = diff.delta / dmag;
        
        return dhat * (eps * (sig - dmag));
    };
};

class SCSpringList {
    private:
        SCatomvec *scs;
        flt eps, sig, l;
    public:
        SCSpringList(SCatomvec *scs, flt eps, flt sig, flt l) : 
            scs(scs), eps(eps), sig(sig), l(l){};
};

#endif
