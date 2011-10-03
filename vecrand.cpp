#include "vecrand.hpp"

engine randengine;
distribution normaldist(0,1);
generator gauss(randengine, normaldist);

Vec rand3d(){
    return Vec(gauss(), gauss(), gauss());
}

void seed(unsigned int n){
    randengine.seed(n);
}

void seed(){
    randengine.seed(static_cast<unsigned int>(time(0)));
}
