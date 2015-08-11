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
#include <vector>
#include <cmath>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <Eigen/Dense>

using namespace std;


#ifdef LONGFLOAT
    typedef long double flt;

    // on MAC OS, g++ is short for clang, and they already define the following functions pretty much
    // the same way I do
    #ifndef __APPLE__
        inline flt cbrt(flt n){return cbrtl(n);};
        inline flt expm1(flt n){return expm1l(n);};
        inline flt copysign(flt n, flt m){return copysignl(n,m);};
        inline flt remainder(flt n, flt m){return remainderl(n,m);};
        inline flt round(flt n){return roundl(n);};
    #endif
#else
    typedef double flt;
#endif

typedef Eigen::Matrix<flt, NDIM, 1> Vec;
typedef Eigen::Matrix<flt, 2, 1> Vec2;
typedef Eigen::Matrix<flt, 3, 1> Vec3;
typedef Eigen::Matrix<flt, 3, 3> Matrix;
typedef Eigen::Matrix<flt, NDIM, 2> VecPair;
using namespace std;

const flt OVERNDIM = ((flt) 1.0)/NDIM;

inline Vec2 vec(double x, double y){
    return Vec(x,y);
};

inline Vec3 vec(double x, double y, double z){
    return Vec(x,y,z);
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
#ifdef VEC3D
Vec randVecSphere(flt radius=1);
#endif

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

long double toLD(double e);
double fromLD(long double e);
vector<long double> LDVector(vector<double> dists);

#endif
