#include "vecrand.hpp"

Pair gauss2() {
    // Taken from Allen and Tildesley, 347
    Pair p;
    
    double x1 = drand48(), x2 = drand48();
    double r = sqrt(-2*log(x1));
    double theta = 2*M_PI*x2;
    p[0] = r*cos(theta);
    p[1] = r*sin(theta);
    return p;
    //~ for (int i=0; i<12; i++){
        //~ p[0] += drand48();
        //~ p[1] += drand48();
    //~ }
    //~ p[0] -= 6;
    //~ p[1] -= 6;
    //~ return p;
};

double gauss() {
    //~ double x=0;
    //~ for (int i=0; i<12; i++){
        //~ x += drand48();
    //~ }
    //~ return x-6;
    
    double x1 = drand48(), x2 = drand48();
    double r = sqrt(-2*log(x1));
    double theta = 2*M_PI*x2;
    return r*cos(theta);
}

Vec rand3d() {
    Pair p1 = gauss2();
    return Vec(p1[0], p1[1], gauss());
}

void bivariateGauss::set(const double s1, const double s2, const double corr){
    // Taken from Allen and Tildesley, 348
    assert(s1 >= 0);
    assert(s2 >= 0);
    assert(corr >= 0);
    assert(corr <= 1);
    x11 = s1;
    x21 = s2 * corr;
    x22 = s2 * sqrt(1 - corr*corr);
}

Pair bivariateGauss::generate(){
    // Taken from Allen and Tildesley, 348
    Pair p = gauss2();
    p[1] = x21*p[0] + x22*p[1];
    p[0] = x11*p[0];
    
    return p;
}

VecPair bivariateGauss::genVecs(){
    VecPair vp;
    //~ vp[0] = Vec(p0[0], p1[0], p2[0]);
    //~ vp[1] = Vec(p0[1], p1[1], p2[1]);
    //~ return vp;
    
    Pair p0=generate(), p1=generate(), p2=generate();
    
    vp[0] = Vec(p0[0], p1[0], p2[0]);
    vp[1] = Vec(p0[1], p1[1], p2[1]);
    return vp;
}
