#ifndef VECRAND_H
#define VECRAND_H

#include <iostream>
#include <ctime>
#include <stdlib.h>

#include "vec.hpp"

typedef Vector<double> Vec;
typedef Numvector<double, 2> Pair;
typedef Nvector<Vec, 2> VecPair;

Vec rand3d();

void seed(unsigned int n){srand48(n);};
void seed(){srand48(time(0));};

class bivariateGauss {
    protected:
        double x11;
        double x21;
        double x22;
        
    public:
        bivariateGauss(const double s1=1, const double s2=1, const double corr=0)
                {set(s1, s2, corr);};
        void set(const double s1, const double s2, const double corr);
        Pair generate();
        VecPair genVecs();
        void seed(uint n){srandom(n);};
        void seed(){srandom(time(0));};
};

//~ int main(){
    //~ randengine.seed(static_cast<unsigned int>(time(0)));
    //~ for(int i=0; i<10; i++) cout << rand3d() << endl;
//~ }

#endif
