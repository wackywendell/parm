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

#include "vec.hpp"

using namespace std;

#ifdef LONGFLOAT
    typedef long double flt;

    // on MAC OS, g++ is short for clang, and they already define the following functions pretty much
    // the same way I do
    #ifndef __clang__
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

#ifdef VEC2D
    typedef Vector2<flt> Vec;
#endif
#ifdef VEC3D
    /**
    The basic physics vector. See Vector3 for full methods.
    */
    typedef Vector3<flt> Vec;
#endif
typedef NumVector<flt, 2> Pair;
typedef NVector<Vec, 2> VecPair;
using namespace std;

/**
A constant equal to \f$\frac{1}{d}\f$, where \f$d\f$ is the number of dimensions.
*/
const flt OVERNDIM = ((flt) 1.0)/NDIM;

/**
A one-liner for creating a Vec object, occasionally useful from within Python.
*/
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

/** Generate a random number between 0 and 1, using the "global" random number generator.
*/
flt rand01();
/** Generate a random vector from a Gaussian distribution, i.e.\
\f$P(x)=\frac{1}{\sqrt{2 \pi}} e^{\frac{-x^2}{2\pi}}\f$, and similarly for \f$y\f$ and \f$z\f$.

In terms of spherical coordinates, directionality is uniform on a sphere, and the radial
distribution is a Chi Distribution with \f$\sigma=1\f$.
*/
Vec randVec();
/** Generate a random vector inside a box with sides of length 1.
*/
Vec randVecBoxed();
#ifdef VEC3D
/** Generate a random vector inside a sphere.
*/
Vec randVecSphere(flt radius=1);
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
        Pair generate();
#ifdef VEC2D
        /**
        Generate a single Vec
        */
        Vec genVec(){return Vec(gauss(), gauss());};
#else
        /**
        Generate a single Vec
        */
        Vec genVec(){return Vec(gauss(), gauss(), gauss());};
#endif
        /** Randomly generate two correlated Vec objects. */
        VecPair genVecs();
};

long double toLD(double e); /**< Go to and from Long Doubles. Useful from Python.*/
double fromLD(long double e); /**< Go to and from Long Doubles. Useful from Python.*/
vector<long double> LDVector(vector<double> dists); /**< Go to and from Long Doubles. Useful from Python.*/

#endif
