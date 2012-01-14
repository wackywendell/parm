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
        inline uint n() const {return num;};
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
        Vec comv() const; //center of mass velocity
        
        //Stats
        flt mass() const;
        flt kinetic(const Vec &originvelocity=Vec(0,0,0)) const;
        Vec momentum() const;
        Vec angmomentum(const Vec &loc) const;
        flt mominertia(const Vec &loc, const Vec &axis) const;
        
        // for resetting
        void addv(Vec v);
        void resetcomv(){addv(-comv());};
        inline void addrot(Vec omega, Vec around);
        inline void resetL(){
            Vec c = com();
            Vec L = angmomentum(c);
            flt I = mominertia(c, L);
            addrot(L / (-I), c);
        }
        
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
        inline atomid get_id(atom *a);
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
        LJforce(const flt epsilon, const flt sigma):epsilon(epsilon),sigma(sigma){};
        inline static flt energy(const Vec diff, const flt sigma, const flt epsilon){
            flt rsq = diff.sq()/(sigma*sigma);
            flt rsix = rsq*rsq*rsq;
            return 4*epsilon*(1/(rsix*rsix) - 1/rsix);
        }
        inline flt energy(const Vec& diff){return energy(diff, sigma, epsilon);};
        inline static Vec forces(const Vec diff, const flt sigma, const flt epsilon){
            flt dsq = diff.sq();
            flt rsq = dsq/(sigma*sigma);
            flt rsix = rsq*rsq*rsq;
            flt fmagTimesR = 4*epsilon*(12/(rsix*rsix) - 6/rsix);
            return diff * (fmagTimesR / dsq);
        }
        inline Vec forces(const Vec& diff){return forces(diff, sigma, epsilon);};
        ~LJforce(){};
};

//~ class LJcutoff {
    //~ protected:
        //~ flt sigma;
        //~ flit epsilon;
        //~ flt cutoff;
        //~ flt cutoffenergy;
    //~ public:
        //~ LJcutoff(const flt epsilon, const flt sigma, const flt cutoff);
        //~ void setcut(const flt cutoff);
        //~ flt energy(const flt diff);
        //~ flt forces(const flt diff);
        //~ ~LJcutoff(){};
//~ };

static const flt LJr0 = pow(2.0, 1.0/6.0);
static const flt LJr0sq = pow(2.0, 1.0/3.0);

class LJrepulsive {
    protected:
        flt sigma;
        flt epsilon;
    public:
        LJrepulsive(const flt epsilon, const flt sigma):
            sigma(sigma), epsilon(epsilon){};
        inline static flt energy(const Vec diff, const flt sig, const flt eps){
            flt rsq = diff.sq()/(sig*sig);
            if(rsq > 1) return 0;
            flt rsix = rsq*rsq*rsq;
            //~ return eps*(4*(1/(rsix*rsix) - 1/rsix) + 1);
            flt mid = (1-1/rsix);
            return eps*(mid*mid);
        };
        inline flt energy(const Vec& diff){return energy(diff, sigma, epsilon);};
        inline static Vec forces(const Vec diff, const flt sig, const flt eps){
            flt dsq = diff.sq();
            flt rsq = dsq/(sig*sig);
            if(rsq > 1) return Vec(0,0,0);
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
        inline Vec forces(const Vec& diff){return forces(diff, sigma, epsilon);};
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
        static flt getcos(const Vec& diff1, const Vec& diff2, const Vec& diff3);
        
        flt energy(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
        Nvector<Vec,4> forces(const Vec& diff1, const Vec& diff2, const Vec& diff3) const;
};

class electricScreened : public interactpair {
    protected:
        flt screen;
        flt q1, q2;
        flt cutoff;
        flt cutoffE;
    public:
        electricScreened(const flt screenLength, const flt q1, 
            const flt q2, const flt cutoff);
        inline flt energy(const Vec &r){return energy(r.mag(),q1*q2,screen, cutoff);};
        static flt energy(const flt r, const flt qaqb, const flt screen, const flt cutoff=0);
        inline Vec forces(const Vec &r){return forces(r,q1*q2,screen, cutoff);};
        static Vec forces(const Vec &r, const flt qaqb, const flt screen, const flt cutoff=0);
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
        flt mean_dists() const;
        flt std_dists() const;
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
        flt mean_dists() const;
        flt std_dists() const;
};

struct dihedralgrouping {
    vector<flt> nums;
    atom *a1, *a2, *a3, *a4;
    dihedralgrouping(vector<flt> nums, atom* a1, atom* a2, atom* a3, atom *a4) : 
                nums(nums), a1(a1), a2(a2), a3(a3), a4(a4){};
};

class dihedrals : public interaction {
    protected:
        vector<dihedralgrouping> groups;
    public:
        dihedrals(vector<dihedralgrouping> pairs = vector<dihedralgrouping>());
        void add(dihedralgrouping b){groups.push_back(b);};
        void add(vector<flt> nums, atom* a1, atom* a2, atom* a3, atom *a4){
        add(dihedralgrouping(nums,a1,a2,a3,a4));};
        uint size() const{ return groups.size();};
        flt mean_dists() const;
        //~ flt std_dists() const;
        flt energy();
        void setForces();
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
            it != pairs.end(); it++) N+= it->second.size();
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
        flt critdist, skinradius;
        vector<atomid> atoms;
        vector<idpair> curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        atomid get_id(atom* a);
        bool ignorechanged;
        //~ bool checkneighbors(const uint n, const uint m) const;
        // this is a full check
    public:
        neighborlist(atomgroup &vec, const flt innerradius, 
        const flt outerradius, pairlist ignore = pairlist());
        void update(){update_list(false);};
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first
        
        inline uint which(){return updatenum;};
        inline uint numpairs(){return curpairs.size();};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b); ignorechanged=true;};
        void ignore(atom*, atom*);
        
        inline void changesize(flt inner, flt outer){
            critdist = inner; skinradius = outer; update_list(true);};
        inline void changesize(flt ratio){
            skinradius = critdist*ratio; update_list(true);};

        inline uint ignore_size() const{return ignorepairs.size();};
        inline vector<idpair>::iterator begin(){return curpairs.begin();};
        inline vector<idpair>::iterator end(){return curpairs.end();};
        
        ~neighborlist(){};
};

template <class A, class P>
class NListed : public interaction {
    protected:
        vector<A> atoms;
        vector<P> pairs;
        neighborlist *neighbors;
        uint lastupdate;
    public:
        NListed(neighborlist *neighbors) : neighbors(neighbors){};
        inline void add(A a){
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            if (a.n() > atoms.size()) atoms.resize(a.n());
            atoms[a.n()] = a;};
        void update_pairs();
        flt energy();
        flt energy_pair(P pair); // This needs to be written!
        void setForces();
        Vec forces_pair(P pair);  // This needs to be written!
        neighborlist *list(){return neighbors;};
        ~NListed(){};
};



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

struct LJatom : public atomid {
    flt sigma, epsilon;
    LJatom() : atomid(), sigma(0), epsilon(0){};
    LJatom(flt sigma, flt epsilon, atomid a) : atomid(a), sigma(sigma),
                    epsilon(epsilon){};
};

struct LJpair {
    flt sigma, epsilon;
    atomid atom1, atom2;
    LJpair(LJatom LJ1, LJatom LJ2) : sigma((LJ1.sigma + LJ2.sigma) / 2),
            epsilon(sqrt(LJ1.epsilon * LJ2.epsilon)),
            atom1(LJ1), atom2(LJ2){};
    inline flt energy(){
        return LJrepulsive::energy((atom1.x()-atom2.x()), sigma,
                                                    epsilon);};
    inline Vec force(){
        return LJrepulsive::forces((atom1.x()-atom2.x()), sigma,
                                                    epsilon);};
};


// For LJgroup
template<>
inline flt NListed<LJatom, LJpair>::energy_pair(LJpair p){
    return LJrepulsive::energy(
                p.atom1.x() - p.atom2.x(), p.sigma, p.epsilon);
};
                
template<>
inline Vec NListed<LJatom, LJpair>::forces_pair(LJpair p){
    return LJrepulsive::forces(
                p.atom1.x() - p.atom2.x(), p.sigma, p.epsilon);
};              

class LJsimple : public interaction {
    protected:
        vector<LJatom> atoms;
        pairlist ignorepairs;
        atomid get_id(atom* a);

    public:
        LJsimple(flt cutoffdist, vector<LJatom> atms=vector<LJatom>());
         // cutoffdist in sigma units
        
        inline void add(LJatom a){
            if (a.n() == atoms.size()) {atoms.push_back(a); return;};
            if (a.n() > atoms.size()) atoms.resize(a.n());
            atoms[a.n()] = a;};
        inline void add(atomid a, flt sigma, flt epsilon){
            add(LJatom(sigma,epsilon,a));};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b);};
        inline void ignore(atom* a, atom* b){
            ignore(get_id(a),get_id(b));};
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint atoms_size() const{return atoms.size();};
        flt energy();
        void setForces();
        //~ ~LJsimple(){};
};

template <class A, class P>
void NListed<A, P>::update_pairs(){
    if(lastupdate == neighbors->which()) return; // already updated
    
    lastupdate = neighbors->which();
    pairs.clear();
    vector<idpair>::iterator pairit;
    for(pairit = neighbors->begin(); pairit != neighbors->end(); pairit++){
        A LJ1 = atoms[pairit->first().n()];
        A LJ2 = atoms[pairit->last().n()];
        pairs.push_back(P(LJ1, LJ2));
    }
}

template <class A, class P>
flt NListed<A, P>::energy(){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    flt E = 0;
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        Vec dist = diff(it->atom1.x(), it->atom2.x());
        E += energy_pair(*it);
    }
    return E;
};

template <class A, class P>
void NListed<A, P>::setForces(){
    update_pairs(); // make sure the LJpairs match the neighbor list ones
    typename vector<P>::iterator it;
    for(it = pairs.begin(); it != pairs.end(); it++){
        Vec f = forces_pair(*it);
        it->atom1.f() += f;
        it->atom2.f() -= f;
    }
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
        flt energy();
        void setForces();
        //~ ~LJsimple(){};
};

#endif
