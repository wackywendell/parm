#include "vecrand.hpp"

#ifndef BOX_H
#define BOX_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>

#define sptr boost::shared_ptr

using namespace boost; // required for SWIG for some reason

typedef const unsigned int cuint;

inline bool toBuffer(vector<Vec*> arr, double* buffer, size_t sizet) {
    if(sizet < NDIM * arr.size()){return false;};
    for(uint i=0; i < arr.size(); i++)
    for(uint j=0; j < NDIM; j++)
    {
        buffer[i*NDIM + j] = (double) ((*arr[i])[j]);
    }
    return true;
};

class atomgroup;

class Box {
    public:
        virtual Vec diff(Vec r1, Vec r2)=0;
        virtual flt V()=0;
        virtual ~Box(){};
};

/***********************************************************************
 * Boxes
 */

#ifdef VEC3D
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remflt(r1[0], r2[0]), remflt(r1[1], r2[1]), remflt(r1[2], r2[2]));
};
#endif
#ifdef VEC2D
inline Vec vecmod(Vec r1, Vec r2){
    return Vec(remflt(r1[0], r2[0]), remflt(r1[1], r2[1]));
};
#endif


class InfiniteBox : public Box {
    public:
        Vec diff(Vec r1, Vec r2){return r1-r2;};
        flt V(){return NAN;};
};

class OriginBox : public Box {
    protected:
        Vec boxsize;
    public:
        OriginBox(Vec size) : boxsize(size){};
        Vec diff(Vec r1, Vec r2){
            return vecmod((r1-r2), boxsize);
        };
        virtual Vec diff(Vec r1, Vec r2, array<int,NDIM> boxes){
            Vec dr = r1 - r2;
            for(uint i=0; i<NDIM; i++) dr[i] -= boxsize[i] * boxes[i];
            return dr;
        };
        virtual array<int,NDIM> box_round(Vec r1, Vec r2){
            array<int,NDIM> boxes;
            Vec dr = r1 - r2;
            for(uint i=0; i<NDIM; i++) boxes[i] = (int) roundflt(dr[i] / boxsize[i]);
            return boxes;
        };
        #ifdef VEC3D
        OriginBox(flt L) : boxsize(L,L,L){};
        flt V(){return boxsize[0] * boxsize[1] * boxsize[2];};
        flt L(){return (boxsize[0] + boxsize[1] + boxsize[2])/3.0;};
        #endif
        #ifdef VEC2D
        OriginBox(flt L) : boxsize(L,L){};
        flt V(){return boxsize[0] * boxsize[1];};
        flt L(){return (boxsize[0] + boxsize[1])/2.0;};
        #endif
        flt resize(flt factor){boxsize *= factor; return V();}
        flt resize(Vec newsize){boxsize = newsize; return V();}
        flt resizeV(flt newV){
            flt curV = V();
            boxsize *= powflt(newV/curV, OVERNDIM);
            return V();
        }
        flt resizeL(flt newL){
            flt curL = powflt(V(), OVERNDIM);
            boxsize *= newL/curL;
            return V();
        }
        Vec randLoc(){
            Vec v = randVecBoxed();
            for(uint i=0; i<NDIM; i++){
                v[i] *= boxsize[i];
            }
            return diff(v, Vec());
        };
        Vec boxshape(){return boxsize;};
};

class LeesEdwardsBox : public OriginBox {
    // Uses shear in the x-direction, relative to y
    protected:
        flt gamma;
    public:
        LeesEdwardsBox(Vec size, flt gamma=0.0) : OriginBox(size), gamma(gamma){};
        Vec diff(Vec r1, Vec r2);
        virtual Vec diff(Vec r1, Vec r2, array<int,NDIM> boxes);
        virtual array<int,NDIM> box_round(Vec r1, Vec r2);
        
        flt get_gamma(){return gamma;};
        
        void shear(flt dgamma, atomgroup &atoms);
        Vec nonaffine(Vec v){
            v[0] -= gamma * v[1];
            return v;
        }
        
        Vec affine(Vec v){
            v[0] += gamma * v[1];
            return v;
        }
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
        inline atom& operator *() const {return *ptr;};
        inline atom *operator->() const{ return ptr;}
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

class idpair : public Array<atomid, 2> {
    public:
        idpair(){vals[0] = atomid(); vals[1] = atomid();};
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

class atomgroup;

class AtomIter{
    private:
        uint i;
        atomgroup &g;
    public:
        AtomIter (atomgroup& g, uint i): i(i), g(g){};
        bool operator!=(const AtomIter& other) const{return i != other.i;};
        atom& operator* () const;
        inline const AtomIter& operator++(){++i; return *this;};
};

class atomvec;

class atomgroup {
    // a group of atoms, such as a molecule, sidebranch, etc.
    public:
        // access individual atoms
        virtual atomvec& vec()=0;
        virtual atom& operator[](cuint n)=0;
        virtual atom& operator[](cuint n) const=0;
        virtual atom& get(cuint n){return ((*this)[n]);};
        virtual atomid get_id(cuint n);
        virtual uint size() const=0;
        virtual flt getmass(const unsigned int n) const {return (*this)[n].m;};
        virtual AtomIter begin(){return AtomIter(*this, 0);};
        virtual AtomIter end(){return AtomIter(*this, size());};
        
        Vec com() const; //center of mass
        Vec comv() const; //center of mass velocity
        
        //Stats
        flt mass() const;
        flt kinetic(const Vec &originvelocity=Vec()) const;
        Vec momentum() const;
        #ifdef VEC3D
        flt moment(const Vec &loc, const Vec &axis, Box &box) const;
        Vec angmomentum(const Vec &loc, Box &box) const;
        Matrix<flt> moment(const Vec &loc, Box &box) const;
        Vec omega(const Vec &loc, Box &box) const;
        void addOmega(Vec w, Vec origin, Box &box);
        inline void resetL(Box &box){
            Vec c = com(), w = omega(c, box);
            if (w.sq() == 0) return;
            addOmega(-w, c, box);
        }
        #elif defined VEC2D
        flt moment(const Vec &loc, Box &box) const;
        flt angmomentum(const Vec &loc, Box &box) const;
        flt omega(const Vec &loc, Box &box) const{return angmomentum(loc, box) / moment(loc, box);};
        void addOmega(flt w, Vec origin, Box &box);
        inline void resetL(Box &box){
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

inline atom& AtomIter::operator*() const{return g[i];};

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
        atomvec& vec(){return *this;};
        inline atom& operator[](cuint n){return atoms[n];};
        inline atom& operator[](cuint n) const {return atoms[n];};
        //atomid get_id(atom *a);
        inline atomid get_id(uint n) {
            if (n > sz) return atomid(); return atomid(atoms + n,n);};
        //~ inline flt getmass(cuint n) const{return atoms[n].m;};
        //~ inline void setmass(cuint n, flt m){atoms[n].m = m;};
        inline uint size() const {return sz;};
        
        ~atomvec(){ delete [] atoms;};
};

class subgroup : public atomgroup {
    protected:
        sptr<atomvec> atoms;
        vector<atomid> ids;
    public:
        subgroup(){};
        //metagroup(vector<atom*> atoms) : atoms(atoms){};
        subgroup(sptr<atomvec> atoms) : atoms(atoms){};
        atomvec& vec(){return *atoms;};
        inline atom& operator[](cuint n){return *ids[n];};
        inline atom& operator[](cuint n) const{return *ids[n];};
        inline atom& get(cuint n){return *ids[n];};
        inline void add(atomid a){return ids.push_back(a);};
        inline atomid get_id(uint n) {return ids[n];};
        inline uint size() const {return (uint) ids.size();};
};

#endif
