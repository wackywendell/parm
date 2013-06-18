#ifndef VECRAND_H
#define VECRAND_H
#ifdef VEC2D
#define NDIM 2
#else
#define VEC3D
#define NDIM 3
#endif
#include <iostream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include "vec.hpp"

#ifdef VEC2D
typedef Vector2<double> Vec;
#else
typedef Vector<double> Vec;
#endif
typedef Numvector<double, 2> Pair;
typedef Nvector<Vec, 2> VecPair;
using namespace std;

typedef boost::mt19937 engine;
typedef boost::normal_distribution<> normdistribution;
typedef boost::uniform_01<> lindistribution;
typedef boost::variate_generator<engine&, normdistribution > normgenerator;
typedef boost::variate_generator<engine&, lindistribution > lingenerator;

Vec randVec();
Vec randVecBoxed();

void seed(unsigned int n);
void seed();

class gaussVec {
    protected:
        engine e;
        normdistribution distro;
        normgenerator gauss;
    public:
        gaussVec(double sigma) : distro(0,sigma), gauss(e,distro){};
        void set(double sigma){distro = normdistribution(0,sigma);};
#ifdef VEC2D
        Vec generate(){return Vec(gauss(),gauss());};
#else
        Vec generate(){return Vec(gauss(),gauss(),gauss());};
#endif
        void seed(unsigned int n){e.seed(n);};
        void seed(){e.seed(static_cast<unsigned int>(time(0)));};
};

class bivariateGauss {
    protected:
        engine e;
        normdistribution distro;
        normgenerator gauss;
        double x11;
        double x21;
        double x22;
        
    public:
        bivariateGauss(const double s1=1, const double s2=1, const double corr=0)
                : distro(0,1), gauss(e,distro){set(s1, s2, corr);};
        void set(const double s1, const double s2, const double corr);
        Pair generate();
#ifdef VEC2D
        Vec genVec(){return Vec(gauss(), gauss());};
#else
        Vec genVec(){return Vec(gauss(), gauss(), gauss());};
#endif
        VecPair genVecs();
        void seed(unsigned int n){e.seed(n);};
        void seed(){e.seed(static_cast<unsigned int>(time(0)));};
};

#endif
