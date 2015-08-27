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
#include <boost/array.hpp>
#include <Eigen/Dense>

using namespace std;
typedef unsigned int uint;

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

typedef std::complex<flt> cmplx; // need the std:: for SWIG complex.i, not sure why

typedef Eigen::Matrix<flt, NDIM, 1> Vec;
typedef Eigen::Matrix<flt, 2, 1> Vec2;
typedef Eigen::Matrix<flt, 3, 1> Vec3;
typedef Eigen::Matrix<flt, 3, 3> Matrix;
typedef Eigen::Matrix<flt, NDIM, 2> VecPair;
using namespace std;

const flt OVERNDIM = ((flt) 1.0)/NDIM;

#ifdef VEC3D
inline Vec vec(){return Vec(0,0,0);};
#endif
#ifdef VEC2D
inline Vec vec(){return Vec(0,0);};
#endif
inline Vec2 vec(double x, double y){return Vec2(x,y);};
inline Vec3 vec(double x, double y, double z){return Vec3(x,y,z);};

inline Vec3 cross(Vec3 v1, Vec3 v2){return v1.cross(v2);};
inline flt cross(Vec2 v1, Vec2 v2){return v1(0)*v2(1) - v2(0)*v1(1);};
inline Vec2 cross(Vec2 v, flt n){return Vec2(v(1)*n, -v(0)*n);};

inline Vec2 perp(Vec2 v){return Vec2(-v(1),v(0));};
inline Vec perpto(Vec r, Vec to){
    return r - to * ((r.dot(to)) / to.squaredNorm());
}

inline Vec2 rotate(Vec2 v, uint i){
    if(i % 4 == 0) return v;
    else if(i % 4 == 3) return Vec2(v(1), -v(0));
    else if(i % 4 == 2) return Vec2(-v(0), -v(1));
    else return Vec2(-v(1), v(0));
}
inline Vec2 flip(Vec2 v){return Vec2(v(1), v(0));};
inline Vec2 rotate_flip(Vec2 v, uint i){
    if((i / 4) % 2 == 1) return rotate(flip(v), i%4);
    return rotate(v, i%4);
};

inline Vec2 rotate_flip_inv(Vec2 v, uint i){
    Vec2 inv = rotate(v, 4-(i%4));
    if((i / 4) % 2 == 0) return inv;
    return flip(inv);
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
        Eigen::Matrix<flt, 1, 2> generate();
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


#ifdef VEC3D
Matrix best_rotation_matrix(Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &from, Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &to);
#endif
#endif
