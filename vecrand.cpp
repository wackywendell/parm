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

void bivariateGauss::set(const double s1, const double s2, const double corr){
    // Taken from Allen and Tildesley, 348
    x11 = s1;
    x21 = s2 * corr;
    x22 = s2 * sqrt(1 - corr*corr);
}

Pair bivariateGauss::generate(){
    double x1 = gauss();
    double x2 = gauss();
    // Taken from Allen and Tildesley, 348
    Pair p;
    p[0] = x11*x1;
    p[1] = x21*x1 + x22*x2;
    return p;
}

VecPair bivariateGauss::genVecs(){
    Vec x1 = Vec(gauss(), gauss(), gauss());
    Vec x2 = Vec(gauss(), gauss(), gauss());
    // Taken from Allen and Tildesley, 348
    VecPair p;
    p[0] = x1*x11;
    p[1] = x1*x21 + x2*x22;
    return p;
}
