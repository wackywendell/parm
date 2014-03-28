#ifndef VECRAND_H
#define VECRAND_H
#ifdef VEC2D
#define NDIM 2
#else
#ifndef VEC3D
#define VEC3D
#endif
#define NDIM 3
#endif
#include <iostream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include "vec.hpp"

#ifdef LONGFLOAT
typedef long double flt;
long double toLD(double e){return (long double) e;};
double fromLD(long double e){return (double) e;};
vector<flt> LDVector(vector<double> dists){
    vector<flt> newdists = vector<flt>();
    for(uint i=0; i<dists.size(); i++){
        newdists.push_back((flt) dists[i]);
    }
    return newdists;
};

inline bool isinfflt(long double n){return isinfl(n);};
inline bool isnanflt(long double n){return isnanl(n);};
inline flt powflt(flt n, flt m){return powl(n,m);};
inline flt sqrtflt(flt n){return sqrtl(n);};
inline flt cbrtflt(flt n){return cbrtl(n);};
inline flt expm1flt(flt n){return expm1l(n);};
inline flt copysignflt(flt n, flt m){return copysignl(n,m);};
inline flt remflt(flt n, flt m){return remainderl(n,m);};
inline flt roundflt(flt n){return roundl(n);};
inline flt ceilflt(flt n){return ceill(n);};
inline flt floorflt(flt n){return floorl(n);};
#else
typedef double flt;
inline bool isinfflt(double n){return isinf(n);};
inline bool isnanflt(double n){return isnan(n);};
inline flt powflt(flt n, flt m){return pow(n,m);};
inline flt sqrtflt(flt n){return sqrt(n);};
inline flt cbrtflt(flt n){return cbrt(n);};
inline flt expm1flt(flt n){return expm1(n);};
inline flt copysignflt(flt n, flt m){return copysign(n,m);};
inline flt remflt(flt n, flt m){return remainder(n,m);};
inline flt roundflt(flt n){return round(n);};
inline flt ceilflt(flt n){return ceil(n);};
inline flt floorflt(flt n){return floor(n);};
#endif

#ifdef VEC2D
typedef Vector2<flt> Vec;
#else
typedef Vector3<flt> Vec;
#endif
typedef Numvector<flt, 2> Pair;
typedef Nvector<Vec, 2> VecPair;
using namespace std;

const flt OVERNDIM = ((flt) 1.0)/NDIM;

inline Vector2<flt> vec(double x, double y){
    return Vector2<flt>(x,y);
};

inline Vector3<flt> vec(double x, double y, double z){
    return Vector3<flt>(x,y,z);
};

inline uint vecsize(){return sizeof(Vec);}

typedef boost::mt19937 engine;
typedef boost::normal_distribution<flt> normdistribution;
typedef boost::uniform_01<flt, flt> lindistribution;
typedef boost::variate_generator<engine&, normdistribution > normgenerator;
typedef boost::variate_generator<engine&, lindistribution > lingenerator;

flt rand01();
Vec randVec();
Vec randVecBoxed();

unsigned int seed(unsigned int n);
unsigned int seed();

class gaussVec {
    protected:
        normdistribution distro;
        normgenerator gauss;
    public:
        gaussVec(flt sigma);
        void set(flt sigma){distro = normdistribution(0,sigma);};
#ifdef VEC2D
        Vec generate(){return Vec(gauss(),gauss());};
#else
        Vec generate(){return Vec(gauss(),gauss(),gauss());};
#endif
};

class bivariateGauss {
    protected:
        normdistribution distro;
        normgenerator gauss;
        flt x11;
        flt x21;
        flt x22;
        
    public:
        bivariateGauss(const flt s1=1, const flt s2=1, const flt corr=0);
        void set(const flt s1, const flt s2, const flt corr);
        Pair generate();
#ifdef VEC2D
        Vec genVec(){return Vec(gauss(), gauss());};
#else
        Vec genVec(){return Vec(gauss(), gauss(), gauss());};
#endif
        VecPair genVecs();
};

#endif
