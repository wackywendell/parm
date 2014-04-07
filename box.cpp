#include "box.hpp"

Vec LeesEdwardsBox::diff(Vec r1, Vec r2){
    flt Ly = boxsize[1];
    flt dy = r1[1]-r2[1];
    int im = (int) roundflt(dy / Ly);
    dy = dy - (im*Ly);
    
    flt Lx = boxsize[0];
    flt dx = r1[0] - r2[0];
    dx = dx - roundflt((dx/Lx)-im*gamma)*Lx-im*gamma*Lx;
    
    #ifdef VEC2D
    return Vec(dx, dy);
    #endif
    #ifdef VEC3D
    flt dz = remflt(r1[2] - r2[2], boxsize[2]);
    return Vec(dx, dy, dz);
    #endif
};

array<int,NDIM> LeesEdwardsBox::box_round(Vec r1, Vec r2){
    Vec dr = r1 - r2;
    array<int,NDIM> boxes;
    boxes[1] = (int) roundflt(dr[1] / boxsize[1]);
    boxes[0] = (int) roundflt((dr[0]/boxsize[0])-boxes[1]*gamma);
    #ifdef VEC3D
    boxes[2] = (int) roundflt(dr[2] / boxsize[2]);
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

Vec atomgroup::com() const{
    flt curmass = 0;
    Vec v = Vec();
    for(unsigned int i=0; i<size(); i++){
        curmass = this->getmass(i);
        v += (*this)[i].x * curmass;
    }
    return v / mass();
};

flt atomgroup::mass() const{
    flt m = 0;
    for(uint i=0; i<size(); i++){
        m += getmass(i);
    }
    return m;
};

Vec atomgroup::comv() const {
    return momentum() / mass();
};

Vec atomgroup::momentum() const{
    flt curmass;
    Vec tot = Vec();
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        tot += (*this)[i].v * curmass;
    }
    return tot;
};
#ifdef VEC3D
Vec atomgroup::angmomentum(const Vec &loc, Box &box) const{
    Vec tot = Vec();
    flt curmass;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box.diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};
#elif defined VEC2D
flt atomgroup::angmomentum(const Vec &loc, Box &box) const{
    flt tot = 0;
    flt curmass;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box.diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};
#endif

#ifdef VEC3D
flt atomgroup::moment(const Vec &loc, const Vec &axis, Box &box) const{
    if (axis.sq() == 0) return 0;
    flt curmass;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box.diff((*this)[i].x, loc).perpto(axis);
        tot += newloc.dot(newloc) * curmass;
    }
    return tot;
};

Matrix<flt> atomgroup::moment(const Vec &loc, Box &box) const{
    flt curmass;
    Matrix<flt> I;
    Vec r;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        r = box.diff((*this)[i].x, loc);
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
    flt curmass;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box.diff((*this)[i].x, loc);
        tot += newloc.dot(newloc) * curmass;
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
    flt curmass;
    flt totE = 0;
    Vec curv;
    for(uint i=0; i<size(); i++){
        curmass = getmass(i);
        curv = (*this)[i].v - originvelocity;
        totE += curmass/2 * curv.dot(curv);
    }
    return totE;
};

void atomgroup::addv(Vec v){
    for(uint i=0; i<size(); i++){
        (*this)[i].v += v;
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

void atomgroup::setAccel(){
    for(uint i=0; i<size(); i++){
        (*this)[i].a = (*this)[i].f / getmass(i);
    }
};

atomid atomvec::get_id(atom* a){
    uint n = (uint) (a - atoms);
    if (n >= sz or a < atoms) return atomid();
    return atomid(atoms + n, n);
};

//~ void atomgroup::vverlet2(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

metagroup::metagroup(vector<atomgroup*> groups){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i < g.size(); i++) atoms.push_back(& (g[i]));
    }
};

metagroup::metagroup(vector<sptr<atomgroup> > groups){
    vector<sptr<atomgroup> >::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i < g.size(); i++) atoms.push_back(& (g[i]));
    }
};

atomid metagroup::get_id(atom* a){
    for(uint n=0; n<atoms.size(); n++)
        if(atoms[n] == a) return atomid(a, n);
    return atomid();
};

void LeesEdwardsBox::shear(flt dgamma, atomgroup &atoms){
    // TODO: is this right for non-square boxes?
    // Is gamma Δx / Lx or Δx / Ly?
    for(uint i=0; i<atoms.size(); i++){
        //flt dy = remflt(atoms[i].x[1], boxsize[1]);
        flt dy = atoms[i].x[1];
        atoms[i].x[0] += dy * dgamma;
    }
    
    gamma += dgamma;
};
