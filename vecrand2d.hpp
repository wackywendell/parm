#ifndef VECRAND_H
#define VECRAND_H
#define VEC2D
#define NDIM 2
#include <iostream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include "vec.hpp"

typedef Vector2<double> Vec;
typedef Numvector<double, 2> Pair;
typedef Nvector<Vec, 2> VecPair;
using namespace std;

typedef boost::mt19937 engine;
typedef boost::normal_distribution<> distribution;
typedef boost::variate_generator<engine&, distribution > generator;

Vec randVec();

void seed(unsigned int n);
void seed();

class gaussVec {
    protected:
        engine e;
        distribution distro;
        generator gauss;
    public:
        gaussVec(double sigma) : distro(0,sigma), gauss(e,distro){};
        void set(double sigma){distro = distribution(0,sigma);};
        Vec generate(){return Vec(gauss(),gauss());};
        void seed(unsigned int n){e.seed(n);};
        void seed(){e.seed(static_cast<unsigned int>(time(0)));};
};

class bivariateGauss {
    protected:
        engine e;
        distribution distro;
        generator gauss;
        double x11;
        double x21;
        double x22;
        
    public:
        bivariateGauss(const double s1=1, const double s2=1, const double corr=0)
                : distro(0,1), gauss(e,distro){set(s1, s2, corr);};
        void set(const double s1, const double s2, const double corr);
        Pair generate();
        Vec genVec(){return Vec(gauss(), gauss());};
        VecPair genVecs();
        void seed(unsigned int n){e.seed(n);};
        void seed(){e.seed(static_cast<unsigned int>(time(0)));};
};

//~ int main(){
    //~ randengine.seed(static_cast<unsigned int>(time(0)));
    //~ for(int i=0; i<10; i++) cout << rand3d() << endl;
//~ }

#endif
