#include "vec.hpp"

#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <boost/foreach.hpp>

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

#define foreach BOOST_FOREACH

inline Vec diff(const Vec a, const Vec b){
    return a-b;
}

struct atom {
    Vec x; // location
    Vec v; // velocity
    Vec a; // acceleration
    Vec f; // forces
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
        inline atomid() : atomref(), num(-1){};
        inline atomid(atom *a, uint n) : atomref(a), num(n){};
        inline uint n(){return num;};
};

class idpair : public array<atomid, 2> {
    public:
        idpair(atomid a, atomid b){ vals[0] = a; vals[1] = b;};
        inline atomid first() const {return vals[0];};
        inline atomid last() const {return vals[1];};
};

class atompair : public array<atom*, 2> {
    public:
        atompair(atom* a, atom* b){ vals[0] = a; vals[1] = b;};
        inline atom& first() const {return *(vals[0]);};
        inline atom& last() const {return *(vals[1]);};
};

class atomtriple : public array<atom*, 3> {
    public:
        atomtriple(atom* a, atom* b, atom* c){vals[0]=a; vals[1]=b; vals[2]=c;};
        inline atom& first() const {return *(vals[0]);};
        inline atom& mid() const {return *(vals[1]);};
        inline atom& last() const {return *(vals[2]);};
};

class atomquad : public array<atom*, 4> {
    public:
        atomquad(atom* a, atom* b, atom* c, atom* d){
                    vals[0]=a; vals[1]=b; vals[2]=c; vals[3]=d;};
        inline atom& first() const {return *(vals[0]);};
        inline atom& mid1() const {return *(vals[1]);};
        inline atom& mid2() const {return *(vals[2]);};
        inline atom& last() const {return *(vals[3]);};
};

class statetracker {
    public:
        virtual void update() = 0;
};

class atomgroup {
    // a group of atoms, such as a molecule, sidebranch, etc.
    public:
        // access individual atoms
        virtual atom& operator[](cuint n)=0;
        virtual atom& operator[](cuint n) const=0;
        atom* get(cuint n){if(n>=size()) return NULL; return &((*this)[n]);};
        atomid get_id(cuint n){return atomid(get(n),n);};
        virtual uint size() const=0;
        virtual flt getmass(const unsigned int n) const = 0;
        
        
        Vec com() const; //center of mass
        Vec comvel() const; //center of mass velocity
        
        //Stats
        flt mass() const;
        flt kinetic(const Vec &originvelocity=Vec(0,0,0)) const;
        Vec momentum() const;
        Vec angmomentum(const Vec &loc) const;
        flt mominertia(const Vec &loc, const Vec &axis) const;
        
        // for timestepping
        void resetForces();
        void setAccel();
        virtual ~atomgroup(){};
};

class atomvec : public virtual atomgroup {
    // this is an atomgroup which actually owns the atoms.
    private:
        atom* atoms;
        vector<flt> ms;
    public:
        atomvec(vector<flt> masses) : ms(masses){atoms = new atom[size()];};
        inline atom& operator[](cuint n){return atoms[n];};
        inline atom& operator[](cuint n) const {return atoms[n];};
        inline flt getmass(cuint n) const{return ms[n];};
        inline void setmass(cuint n, flt m){ms[n] = m;};
        inline uint size() const {return ms.size();};
        ~atomvec(){ delete [] atoms;};
};

class metagroup : public atomgroup {
    protected:
        vector<atom*> atoms;
        vector<flt> masses;
    public:
        metagroup(vector<atomgroup*>);
        inline atom& operator[](cuint n){return *atoms[n];};
        inline atom& operator[](cuint n) const{return *atoms[n];};
        inline atom* get(cuint n){if(n>=size()) return NULL; return (atoms[n]);};
        inline uint size() const {return atoms.size();};
        inline flt getmass(const unsigned int n) const{return masses[n];};
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
        virtual flt energy(const flt diff);
        inline flt energy(const Vec& diff){return energy(diff.mag());};
        virtual flt forces(const flt diff);
        inline Vec forces(const Vec& diff){flt m=diff.mag(); return diff *(forces(m)/m);};
        ~LJforce(){};
};

class LJcutoff : public LJforce {
    protected:
        flt cutoff;
        flt cutoffenergy;
    public:
        LJcutoff(const flt epsilon, const flt sigma, const flt cutoff);
        void setcut(const flt cutoff);
        flt energy(const flt diff);
        flt forces(const flt diff);
        ~LJcutoff(){};
};

class LJcutrepulsive : public LJforce {
    protected:
        flt cutoff;
        flt cutoffenergy;
    public:
        LJcutrepulsive(const flt epsilon, const flt sigma, const flt cutoff);
        void setcut(const flt cutoff);
        flt energy(const flt diff);
        flt forces(const flt diff);
        ~LJcutrepulsive(){};
};

class spring : public interactpair {
    protected:
        flt springk;
        flt x0;
    public:
        spring(const flt k, const flt x0) : springk(k),x0(x0){};
        flt energy(const Vec& diff);
        Vec forces(const Vec& diff);
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
        flt energy(const Vec& diff1, const Vec& diff2);
        Nvector<Vec,3> forces(const Vec& diff1, const Vec& diff2);
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
        flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
        Nvector<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
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

class singletpairs : public interaction {
    protected:
        interactpair* inter;
        vector<atompair> atoms;
    public:
        singletpairs(interactpair* inter, vector<atompair> atoms 
                    = vector<atompair>()) : inter(inter), atoms(atoms){};
        void add(atompair as){atoms.push_back(as);};
        void add(atom* a, atom* b){atoms.push_back(atompair(a,b));};
        flt energy();
        void setForces();
        ~singletpairs(){};
};

class interactgroup : public interaction {
    protected:
        vector<interaction*> inters;
    public:
        interactgroup(vector<interaction*> inters=vector<interaction*>())
                    : inters(inters){};
        void add(interaction* a){inters.push_back(a);};
        uint size() const{ return inters.size();};
        flt energy();
        void setForces();
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
    public:
        bondpairs(vector<bondgrouping> pairs = vector<bondgrouping>());
        void add(bondgrouping b){pairs.push_back(b);};
        void add(flt k, flt x0, atom* a1, atom* a2){add(bondgrouping(k,x0,a1,a2));};
        uint size() const{ return pairs.size();};
        flt energy();
        void setForces();
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
    public:
        angletriples(vector<anglegrouping> triples = vector<anglegrouping>());
        void add(anglegrouping b){triples.push_back(b);};
        void add(flt k, flt x0, atom* a1, atom* a2, atom* a3){
                                add(anglegrouping(k,x0,a1,a2,a3));};
        inline flt energy();
        inline void setForces();
};

//~ class LJgroup : private vector<LJdata>, public atomgroup {
    //~ public:
        //~ LJgroup() : vector<LJdata>(){};
        //~ add(flt sigma, flt epsilon, atom *a){ 
//~ };

struct atompaircomp {
    bool operator() (const atompair& lhs, const atompair& rhs) const{
        return (lhs[0] == rhs[0] and lhs[1] == rhs[1]);}
};

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
            typename vector<atomid>::iterator it;
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
        flt critdist, skinradius;
        vector<atomid> atoms;
        vector<idpair> curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        atomid get_id(atom* a);
        //~ bool checkneighbors(const uint n, const uint m) const;
        // this is a full check
    public:
        neighborlist(atomgroup &vec, const flt innerradius, 
        const flt outerradius, pairlist ignore = pairlist());
        void update(){update_list(false);};
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first
        
        inline uint which(){return updatenum;};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        void ignore(atom*, atom*);
        inline vector<idpair>::iterator begin(){return curpairs.begin();};
        inline vector<idpair>::iterator end(){return curpairs.end();};
        
        ~neighborlist(){};
};

struct LJatom : public atomid {
    flt sigma, epsilon;
    LJatom() : sigma(0), epsilon(0){};
    LJatom(flt sigma, flt epsilon, atomid a) : sigma(sigma), epsilon(epsilon),
        atomid(a){};
};

struct LJpair {
    flt sigma, epsilon;
    atomid atom1, atom2;
    LJpair(LJatom LJ1, LJatom LJ2) : sigma((LJ1.sigma + LJ2.sigma) / 2),
            epsilon(sqrt(LJ1.epsilon * LJ2.epsilon)),
            atom1(LJ1), atom2(LJ2){};
};

class LJgroup : public interaction {
    protected:
        LJcutrepulsive LJ; 
        vector<LJatom> atoms;
        vector<LJpair> pairs;
        neighborlist *neighbors;
        uint lastupdate;
    public:
        LJgroup(neighborlist *neighbors, flt cutoffdist); // cutoffdist in sigma units
        inline void add(LJatom a){
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            if (a.n() > atoms.size()) atoms.resize(a.n());
            atoms[a.n()] = a;};
        inline void add(atomid a, flt sigma, flt epsilon){
            add(LJatom(sigma,epsilon,a));};
        void update_pairs();
        flt energy();
        void setForces();
        neighborlist *list(){return neighbors;};
        ~LJgroup(){};
};

#endif
