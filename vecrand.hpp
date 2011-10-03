#ifndef VECRAND_H
#define VECRAND_H

#include <iostream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include "vec.hpp"

typedef Vector<double> Vec;
using namespace std;

typedef boost::mt19937 engine;
typedef boost::normal_distribution<> distribution;
typedef boost::variate_generator<engine&, distribution > generator;

Vec rand3d();

void seed(unsigned int n);
void seed();

class gaussVec {
    protected:
        engine e;
        distribution distro;
        generator gauss;
    public:
        gaussVec(double sigma) : distro(0,sigma/sqrt(3)), gauss(e,distro){};
        Vec generate(){return Vec(gauss(),gauss(),gauss());};
        void seed(unsigned int n){e.seed(n);};
        void seed(){e.seed(static_cast<unsigned int>(time(0)));};
};

//~ int main(){
    //~ randengine.seed(static_cast<unsigned int>(time(0)));
    //~ for(int i=0; i<10; i++) cout << rand3d() << endl;
//~ }

#endif
