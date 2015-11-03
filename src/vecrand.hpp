#ifndef VECRAND_H
#define VECRAND_H
#ifdef VEC2D
#define NDIM 2
#else
#ifndef VEC3D
#define VEC3D
#endif
#endif
#ifdef VEC3D
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
    /**
    The basic floating point type used in the simulations. The entire package can be compiled with
    LONGFLOAT defined to use long doubles.
    */
    typedef double flt;
#endif

typedef std::complex<flt> cmplx; // need the std:: for SWIG complex.i, not sure why

    /**
    The basic physics vector. See Vector3 for full methods.
    */
typedef Eigen::Matrix<flt, NDIM, 1> Vec;
typedef Eigen::Matrix<flt, 2, 1> Vec2;
typedef Eigen::Matrix<flt, 3, 1> Vec3;
typedef Eigen::Matrix<flt, NDIM, NDIM> Matrix;
typedef Eigen::Matrix<flt, 2, 2> Matrix2;
typedef Eigen::Matrix<flt, 3, 3> Matrix3;
typedef Eigen::Matrix<flt, NDIM, 2> VecPair;

using namespace std;

/**
A constant equal to \f$\frac{1}{d}\f$, where \f$d\f$ is the number of dimensions.
*/
const flt OVERNDIM = ((flt) 1.0)/NDIM;

/**
A one-liner for creating a Vec object, occasionally useful from within Python.
*/
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
};
inline Vec2 rotate_inv(Vec2 v, uint i){
    return rotate(v, 4 - i);
};
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

inline Vec3 rotate(Vec3 v, uint i){
    i %= 24;
    if(i == 0) return v;
    else if(i ==  1) return Vec3( v(0),-v(1),-v(2));
    else if(i ==  2) return Vec3( v(1), v(0),-v(2));
    else if(i ==  3) return Vec3( v(2),-v(1), v(0));
    else if(i ==  4) return Vec3(-v(2),-v(1),-v(0));
    else if(i ==  5) return Vec3(-v(1),-v(0),-v(2));
    else if(i ==  6) return Vec3(-v(0), v(1),-v(2));
    else if(i ==  7) return Vec3(-v(0), v(2), v(1));
    else if(i ==  8) return Vec3(-v(0),-v(2),-v(1));
    else if(i ==  9) return Vec3(-v(0),-v(1), v(2));
    else if(i == 10) return Vec3( v(0), v(2),-v(1));
    else if(i == 11) return Vec3( v(1), v(2), v(0));
    else if(i == 12) return Vec3( v(1),-v(2),-v(0));
    else if(i == 13) return Vec3( v(1),-v(0), v(2));
    else if(i == 14) return Vec3( v(2), v(1),-v(0));
    else if(i == 15) return Vec3( v(2),-v(0),-v(1));
    else if(i == 16) return Vec3(-v(2),-v(0), v(1));
    else if(i == 17) return Vec3(-v(1), v(2),-v(0));
    else if(i == 18) return Vec3(-v(1),-v(2), v(0));
    else if(i == 19) return Vec3(-v(2), v(1), v(0));
    else if(i == 20) return Vec3(-v(1), v(0), v(2));
    else if(i == 21) return Vec3(-v(2), v(0),-v(1));
    else if(i == 22) return Vec3( v(2), v(0), v(1));
    else if(i == 23) return Vec3( v(0),-v(2), v(1));
    throw std::runtime_error("This should be impossible to reach.");
};

inline Vec3 rotate_inv(Vec3 v, uint i){
    i %= 24;
    if(i < 10) return rotate(v, i);
    else return rotate(v, 33 - i);
};
inline Vec3 flip(Vec3 v){return -v;};
inline Vec3 rotate_flip(Vec3 v, uint i){
    if((i / 24) % 2 == 1) return rotate(flip(v), i%24);
    return rotate(v, i%24);
};

inline Vec3 rotate_flip_inv(Vec3 v, uint i){
    Vec3 inv = rotate_inv(v, i%24);
    if((i / 24) % 2 == 0) return inv;
    return flip(inv);
};


inline uint vecsize(){return sizeof(Vec);}

typedef boost::mt19937 engine;
typedef boost::normal_distribution<flt> normdistribution;
typedef boost::uniform_01<flt, flt> lindistribution;
typedef boost::variate_generator<engine&, normdistribution > normgenerator;
typedef boost::variate_generator<engine&, lindistribution > lingenerator;

/** Generate a random number between 0 and 1, using the "global" random number generator.
*/
flt rand01();
/** Generate a random vector from a Gaussian distribution, i.e.\
\f$P(x)=\frac{1}{\sqrt{2 \pi}} e^{\frac{-x^2}{2\pi}}\f$, and similarly for \f$y\f$ and \f$z\f$.

In terms of spherical coordinates, directionality is uniform on a sphere, and the radial
distribution is a Chi Distribution with \f$\sigma=1\f$.
*/
Vec rand_vec();
/** Generate a random vector inside a box with sides of length 1.
*/
Vec rand_vec_boxed();
#ifdef VEC3D
/** Generate a random vector inside a sphere.
*/
Vec rand_vec_sphere(flt radius=1);
#endif

/** Seed the global random number generator with a given integer.
*/
unsigned int seed(unsigned int n);
/** Seed the global random number generator with the current time.

\returns the seed used.
*/
unsigned int seed();

class GaussVec {
    protected:
        normdistribution distro;
        normgenerator gauss;
    public:
        GaussVec(flt sigma);
        void set(flt sigma){distro = normdistribution(0,sigma);};
#ifdef VEC2D
        Vec generate(){return Vec(gauss(),gauss());};
#else
        Vec generate(){return Vec(gauss(),gauss(),gauss());};
#endif
};

/**
A class for generating two random numbers from a Gaussian distribution, with a given correlation.

This is used by some of the integrators.
*/
class BivariateGauss {
    protected:
        normdistribution distro;
        normgenerator gauss;
        flt x11;
        flt x21;
        flt x22;

    public:
        /**
        See `set()` for parameters.
        */
        BivariateGauss(const flt s1=1, const flt s2=1, const flt corr=0);
        /**
        Set the standard deviations and correlations.

        \param s1 Standard deviation of the first vector. s1 > 0.
        \param s2 Standard deviation of the second vector. s2 > 0.
        \param corr Normalized correlation between the two vectors, 0 <= corr <= 1.
        */
        void set(const flt s1, const flt s2, const flt corr);
        /** Randomly generate two correlated numbers. */
        Eigen::Matrix<flt, 1, 2> generate();
#ifdef VEC2D
        /**
        Generate a single Vec
        */
        Vec gen_vec(){return Vec(gauss(), gauss());};
#else
        /**
        Generate a single Vec
        */
        Vec gen_vec(){return Vec(gauss(), gauss(), gauss());};
#endif
        /** Randomly generate two correlated Vec objects. */
        VecPair gen_vecs();
};

long double to_LD(double e); /**< Go to and from Long Doubles. Useful from Python.*/
double from_LD(long double e); /**< Go to and from Long Doubles. Useful from Python.*/
vector<long double> LDVector(vector<double> dists); /**< Go to and from Long Doubles. Useful from Python.*/

/// Is there a NaN in this matrix?
/// Eigen has an allFinite and hasNan function as of version 3.2, but
/// Ubuntu libeigen3-dev is not up that far.
template<typename Derived>
inline bool hasNaN(const Eigen::DenseBase<Derived> &m) {
    // return m.hasNan();
    return !((m.derived().array()==m.derived().array()).all());
}


/// Are all values finite, i.e. not NaN and not +/- Inf?
/// Eigen has an allFinite and hasNan function as of version 3.2, but
/// Ubuntu libeigen3-dev is not up that far.
template<typename Derived>
inline bool allFinite(const Eigen::DenseBase<Derived> &m){
    // return m.allFinite();
    return !hasNaN(((m.derived()-m.derived())));
}

template<typename Derived>
void finite_or_throw(const Eigen::DenseBase<Derived> &m){
    if(!allFinite(m)) {
        std::cerr << "Matrix !Finite ERROR" << std::endl;
        std::cerr << m << std::endl;
        throw std::invalid_argument("Matrix was not finite, cannot continue.");
    }
}


#ifdef VEC3D
Matrix best_rotation_matrix(Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &from, Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &to);
#endif
#endif
