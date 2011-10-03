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
    Vec x;
    Vec v;
    Vec a;
    Vec f;
};

class atomgroup {
    public:
        virtual atom& operator[](cuint n)=0;
        virtual atom& operator[](cuint n) const=0;
        virtual uint N() const=0;
        Vec com() const; //center of mass
        Vec comvel() const; //center of mass velocity
        virtual flt getmass(const unsigned int n) const = 0;
        virtual Vec diff(cuint n, cuint m) const = 0;
        virtual Vec diff(cuint n, const Vec r) const = 0;
        flt mass() const;
        flt kinetic(const Vec &originvelocity=Vec(0,0,0)) const;
        Vec momentum() const;
        Vec angmomentum(const Vec &loc) const;
        flt mominertia(const Vec &loc, const Vec &axis) const;
        void resetForces();
        void setAccel();
        void vverlet1(const flt dt);
        void vverlet2(const flt dt);
        atom* get(cuint n){if(n>=N()) return NULL; return &((*this)[n]);};
        ~atomgroup(){};
};

class atomgroupL : public atomgroup {
    private:
        const flt L;
    public:
        atomgroupL(const flt Length) : L(Length){};
        Vec diff(cuint n, cuint m) const;
        Vec diff(cuint n, const Vec r) const;
};

class atomvec : public atomgroupL {
    private:
        atom* atoms;
        vector<flt> masses;
    public:
        atomvec(const flt L, vector<flt> masses);
        atom& operator[](cuint n){return atoms[n];};
        atom& operator[](cuint n) const {return atoms[n];};
        flt getmass(cuint n) const{return masses[n];};
        void setmass(cuint n, flt m){masses[n] = m;};
        uint N() const {return masses.size();};
        ~atomvec(){ delete [] atoms;};
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
        interactquad* trip;
    public:
        intraMolNNQuad(vector<atomgroup*> groupvec, interactquad* trip);
        flt energy();
        void setForces();
};
#endif
