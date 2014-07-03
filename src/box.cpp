#include "box.hpp"

Vec LeesEdwardsBox::diff(Vec r1, Vec r2){
    flt Ly = boxsize[1];
    flt dy = r1[1]-r2[1];
    int im = (int) round(dy / Ly);
    dy = dy - (im*Ly);
    
    flt Lx = boxsize[0];
    flt dx = r1[0] - r2[0];
    dx = dx - round((dx/Lx)-im*gamma)*Lx-im*gamma*Lx;
    
    #ifdef VEC2D
    return Vec(dx, dy);
    #endif
    #ifdef VEC3D
    flt dz = remainder(r1[2] - r2[2], boxsize[2]);
    return Vec(dx, dy, dz);
    #endif
};

array<int,NDIM> LeesEdwardsBox::box_round(Vec r1, Vec r2){
    Vec dr = r1 - r2;
    array<int,NDIM> boxes;
    boxes[1] = (int) round(dr[1] / boxsize[1]);
    boxes[0] = (int) round((dr[0]/boxsize[0])-boxes[1]*gamma);
    #ifdef VEC3D
    boxes[2] = (int) round(dr[2] / boxsize[2]);
    #endif
    return boxes;
};

Vec LeesEdwardsBox::diff(Vec r1, Vec r2, array<int,NDIM> boxes){
    Vec dr = r1 - r2;
    dr[0] -= (boxes[0] + boxes[1]*gamma)*boxsize[0];
    dr[1] -= boxes[1]*boxsize[1];
    #ifdef VEC3D
    dr[2] -= boxes[2]*boxsize[2];
    #endif
    return dr;
}

SCbox::SCbox(flt L, flt R) : L(L), R(R){};

flt SCbox::V(){
    #ifdef VEC2D
    return R*L + R*R*M_PI;
    #else
    flt endcap = R*R*M_PI;
    return endcap*L + endcap*R*4/3;
    #endif
};

Vec SCbox::dist(Vec r1){
    // Vector to nearest point on central line
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;
    return r1;
};

Vec SCbox::edgedist(Vec r1){
    // Vector to nearest edge
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;
    
    flt dmag = r1.mag(); // distance to center
    if(dmag == 0){
        #ifdef VEC2D
        return Vec(0, R);
        #else
        return Vec(0, R, 0);
        #endif
    }
    return r1 * ((R-dmag)/dmag);
}

bool SCbox::inside(Vec r1, flt buffer){
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;
    flt newR = R - buffer;
    return r1.sq() < newR*newR;
};

Vec SCbox::randLoc(flt min_dist_to_wall){
    if(min_dist_to_wall >= R) return Vec();
    Vec v;
    flt Rmin = R - min_dist_to_wall;
    flt Rminsq = pow(Rmin, 2.0);
    while(true){
        v = randVecBoxed();
        v[0] -= 0.5;
        v[0] *= (L + 2*Rmin);
        
        v[1] -= 0.5;
        v[1] *= 2*Rmin;
        
        #ifndef VEC2D
        v[2] -= 0.5;
        v[2] *= 2*Rmin;
        #endif
        
        flt distsq = pow(v[1],2);
        if(abs(v[0]) >= L/2.0) distsq += pow(abs(v[0]) - L/2.0, 2.0);
        
        
        #ifndef VEC2D
        distsq += pow(v[2],2);
        #endif
        
        if(distsq <= Rminsq) break;
    };
    return v;
};

Vec atomgroup::com() const{
    Vec v = Vec();
    for(unsigned int i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        v += (*this)[i].x * curmass;
    }
    return v / mass();
};

flt atomgroup::mass() const{
    flt m = 0;
    for(uint i=0; i<size(); i++){
        m += (*this)[i].m;
    }
    return m;
};

Vec atomgroup::comv() const {
    return momentum() / mass();
};

Vec atomgroup::momentum() const{
    Vec tot = Vec();
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        tot += (*this)[i].v * curmass;
    }
    return tot;
};

flt atomgroup::gyradius() const{
    Vec avgr = Vec();
    for(uint i = 0; i<size(); i++){
        avgr += (*this)[i].x;
    }
    avgr /= size(); // now avgr is the average location, akin to c.o.m.
    flt Rgsq = 0;
    for(uint i = 0; i<size(); i++){
        Rgsq += ((*this)[i].x - avgr).sq();
    }
    
    return sqrt(Rgsq/size());
};

#ifdef VEC3D
Vec atomgroup::angmomentum(const Vec &loc, Box &box) const{
    Vec tot = Vec();
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        Vec newloc = box.diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};
#elif defined VEC2D
flt atomgroup::angmomentum(const Vec &loc, Box &box) const{
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        newloc = box.diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * (*this)[i].m; // r x v m = r x p
    }
    return tot;
};
#endif

#ifdef VEC3D
flt atomgroup::moment(const Vec &loc, const Vec &axis, Box &box) const{
    if (axis.sq() == 0) return 0;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        newloc = box.diff((*this)[i].x, loc).perpto(axis);
        tot += newloc.dot(newloc) * (*this)[i].m;
    }
    return tot;
};

Matrix<flt> atomgroup::moment(const Vec &loc, Box &box) const{
    Matrix<flt> I;
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        Vec r = box.diff((*this)[i].x, loc);
        flt x = r.getx(), y = r.gety(), z = r.getz();
        I[0][0] += curmass * (y*y + z*z);
        I[1][1] += curmass * (x*x + z*z);
        I[2][2] += curmass * (x*x + y*y);
        I[0][1] -= curmass * (x*y);
        I[0][2] -= curmass * (x*z);
        I[1][2] -= curmass * (y*z);
    }
    I[1][0] = I[0][1];
    I[2][0] = I[0][2];
    I[2][1] = I[1][2];
    return I;
};

Vec atomgroup::omega(const Vec &loc, Box &box) const{
    Matrix<flt> Inv = moment(loc, box).SymmetricInverse();
    return Inv * (angmomentum(loc, box));
};

void atomgroup::addOmega(Vec w, Vec loc, Box &box){
    for(uint i=0; i<size(); i++){
        Vec r = box.diff((*this)[i].x, loc);
        (*this)[i].v -= r.cross(w);
    }
};
#elif defined VEC2D
flt atomgroup::moment(const Vec &loc, Box &box) const{
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        newloc = box.diff((*this)[i].x, loc);
        tot += newloc.dot(newloc) * (*this)[i].m;
    }
    return tot;
};

void atomgroup::addOmega(flt w, Vec loc, Box &box){
    for(uint i=0; i<size(); i++){
        Vec r = box.diff((*this)[i].x, loc);
        (*this)[i].v -= r.perp().norm()*w;
    }
};
#endif

flt atomgroup::kinetic(const Vec &originvelocity) const{
    flt totE = 0;
    Vec curv;
    for(uint i=0; i<size(); i++){
        curv = (*this)[i].v - originvelocity;
        totE += (*this)[i].m/2 * curv.dot(curv);
    }
    return totE;
};

void atomgroup::addv(Vec v){
    for(uint i=0; i<size(); i++){
        (*this)[i].v += v;
    }
};

void atomgroup::randomize_velocities(flt T){
    for(uint i=0; i<size(); i++){
        (*this)[i].v = randVec() * sqrt(T*2.0/(*this)[i].m);
    }
};

void atomgroup::resetForces(){
    for(uint i=0; i<size(); i++){
        (*this)[i].f = Vec();
    }
};

//~ void atomgroup::vverlet1(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].x += (*this)[i].v * dt + (*this)[i].a * (dt*dt/2);
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

//~ void atomgroup::vverlet2(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

void LeesEdwardsBox::shear(flt dgamma, atomgroup &atoms){
    // TODO: is this right for non-square boxes?
    // Is gamma Δx / Lx or Δx / Ly?
    for(uint i=0; i<atoms.size(); i++){
        //flt dy = remainder(atoms[i].x[1], boxsize[1]);
        flt dy = atoms[i].x[1];
        atoms[i].x[0] += dy * dgamma;
    }
    
    gamma += dgamma;
};
