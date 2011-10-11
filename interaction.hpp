#include "vec.hpp"

#include <vector>

#ifndef INTERACTION_H
#define INTERACTION_H

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

typedef double flt;
typedef unsigned int uint;
typedef const unsigned int cuint;
typedef Vector<flt> Vec;

struct atom {
    Vec x; // location
    Vec v; // velocity
    Vec a; // acceleration
    Vec f; // forces
};

class atomgroup {
    // a group of atoms, such as a molecule, sidebranch, etc.
    public:
        // access individual atoms
        virtual atom& operator[](cuint n)=0;
        virtual atom& operator[](cuint n) const=0;
        atom* get(cuint n){if(n>=N()) return NULL; return &((*this)[n]);};
        virtual uint N() const=0;
        virtual flt getmass(const unsigned int n) const = 0;
        
        
        Vec com() const; //center of mass
        Vec comvel() const; //center of mass velocity
        virtual Vec diff(cuint n, cuint m) const;
        virtual Vec diff(cuint n, const Vec r) const;
        
        //Stats
        flt inline mass() const;
        flt kinetic(const Vec &originvelocity=Vec(0,0,0)) const;
        Vec momentum() const;
        Vec angmomentum(const Vec &loc) const;
        flt mominertia(const Vec &loc, const Vec &axis) const;
        
        // for timestepping
        void resetForces();
        void setAccel();
        virtual ~atomgroup(){};
};

class constL : public virtual atomgroup {
    // Keeps track of an atomgroup inside a fixed size box of length L
    private:
        const flt L;
    public:
        constL(const flt Length) : L(Length){};
        Vec diff(cuint n, cuint m) const;
        Vec diff(cuint n, const Vec r) const;
};

class atomvec : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms.
    private:
        atom* atoms;
        vector<flt> ms;
    public:
        atomvec(vector<flt> masses) : ms(masses){atoms = new atom[N()];};
        inline atom& operator[](cuint n){return atoms[n];};
        inline atom& operator[](cuint n) const {return atoms[n];};
        inline flt getmass(cuint n) const{return ms[n];};
        inline void setmass(cuint n, flt m){ms[n] = m;};
        inline uint N() const {return ms.size();};
        ~atomvec(){ delete [] atoms;};
};

class atomvecL : public constL, public atomvec {
    public:
        atomvecL(vector<flt> masses, const flt Length) :
                constL(Length), atomvec(masses){};
        ~atomvecL(){};
};

class interactpair {
    public:
        virtual flt energy(const Vec& diff)=0;
        virtual Vec forces(const Vec& diff)=0;
        virtual ~interactpair(){};
};

class interacttriple {
    public:
        virtual flt energy(const Vec& diff1, const Vec& diff2)=0;
        virtual Nvector<Vec,3> forces(const Vec& diff1, const Vec& diff2)=0;
        virtual ~interacttriple(){};
};

class interactquad {
    public:
        virtual flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const=0;
        virtual Nvector<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const=0;
        virtual ~interactquad(){};
};

template <uint N>
class interactN {
    public:
        virtual flt energy(const Vec* location[N])=0;
        virtual Nvector<Vec, N> forces(const Vec* location[N])=0;
        
        virtual ~interactN<N>(){};
};

class LJforce : public interactpair {
    protected:
        flt epsilon;
        flt sigma;
    public:
        LJforce(const flt epsilon, const flt sigma);
        virtual flt energy(const Vec& diff);
        virtual Vec forces(const Vec& diff);
        ~LJforce(){};
};

class LJcutoff : public LJforce {
    protected:
        flt cutoff;
        flt cutoffenergy;
    public:
        LJcutoff(const flt epsilon, const flt sigma, const flt cutoff);
        void setcut(const flt cutoff);
        virtual flt energy(const Vec& diff);
        virtual Vec forces(const Vec& diff);
        ~LJcutoff(){};
};

class spring : public interactpair {
    protected:
        flt springk;
        flt x0;
    public:
        spring(const flt k, const flt x0) : springk(k),x0(x0){};
        virtual flt energy(const Vec& diff);
        virtual Vec forces(const Vec& diff);
        ~spring(){};
};

class bondangle : public interacttriple {
    // two vectors for E and F are from the central one to the other 2
    protected:
        flt springk;
        flt theta0;
        bool usecos;
    public:
        bondangle(const flt k, const flt theta, const bool cosine=false)
                        :springk(k), theta0(theta), usecos(cosine){};
        virtual flt energy(const Vec& diff1, const Vec& diff2);
        virtual Nvector<Vec,3> forces(const Vec& diff1, const Vec& diff2);
        ~bondangle(){};
};

class dihedral : public interactquad {
    // vectors are from 1 to 2, 2 to 3, 3 to 4
    protected:
        vector<flt> torsions;
        Nvector<Vec,4> derivs(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
        flt dudcostheta(const flt costheta) const;
    public:
        dihedral(const vector<flt> vals);
        virtual flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
        virtual Nvector<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
};

class interaction {
    public:
        virtual flt energy()=0;
        virtual void setForces()=0;
};

class interMolPair : public interaction {
    private:
        vector<atomgroup*> groups;
        interactpair* pair;
    public:
        interMolPair(vector<atomgroup*> groupvec, interactpair* pair);
        flt energy();
        void setForces();
};

class intraMolNNPair : public interaction {
    private:
        vector<atomgroup*> groups;
        interactpair* pair;
    public:
        intraMolNNPair(vector<atomgroup*> groupvec, interactpair* pair);
        flt energy();
        void setForces();
};

class intraMolPairs : public interaction {
    private:
        vector<atomgroup*> groups;
        interactpair* pair;
        cuint skip;
    public:
        intraMolPairs(vector<atomgroup*> groupvec, interactpair* pair, cuint skip);
        flt energy();
        void setForces();
};

class intraMolNNTriple : public interaction {
    private:
        vector<atomgroup*> groups;
        interacttriple* trip;
    public:
        intraMolNNTriple(vector<atomgroup*> groupvec, interacttriple* trip);
        flt energy();
        void setForces();
};

class intraMolNNQuad : public interaction {
    private:
        vector<atomgroup*> groups;
        interactquad* quad;
    public:
        intraMolNNQuad(vector<atomgroup*> groupvec, interactquad* quad);
        flt energy();
        void setForces();
};
#endif
