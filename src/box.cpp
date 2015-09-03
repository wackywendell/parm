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

SCBox::SCBox(flt L, flt R) : L(L), R(R){};

flt SCBox::V(){
    #ifdef VEC2D
    return R*L + R*R*M_PI;
    #else
    flt endcap = R*R*M_PI;
    return endcap*L + endcap*R*4/3;
    #endif
};

Vec SCBox::dist(Vec r1){
    // Vector to nearest point on central line
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;
    return r1;
};

Vec SCBox::edge_dist(Vec r1){
    // Vector to nearest edge
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;

    flt dmag = r1.norm(); // distance to center
    if(dmag == 0){
        #ifdef VEC2D
        return Vec(0, R);
        #else
        return Vec(0, R, 0);
        #endif
    }
    return r1 * ((R-dmag)/dmag);
}

bool SCBox::inside(Vec r1, flt buffer){
    if(r1[0] < -L/2) r1[0] += L/2;
    else if(r1[0] > L/2) r1[0] -= L/2;
    else r1[0] = 0;
    flt newR = R - buffer;
    return r1.squaredNorm() < newR*newR;
};

Vec SCBox::rand_loc(flt min_dist_to_wall){
    if(min_dist_to_wall >= R) return Vec::Zero();
    Vec v = Vec::Zero();
    flt Rmin = R - min_dist_to_wall;
    flt Rminsq = pow(Rmin, 2.0);
    while(true){
        v = rand_vec_boxed();
        v[0] -= 0.5;
        v[0] *= (L + 2*Rmin);

        v[1] -= 0.5;
        v[1] *= 2*Rmin;

        #ifndef VEC2D
        v[2] -= 0.5;
        v[2] *= 2*Rmin;
        #endif

        flt distance_squared = pow(v[1],2);
        if(abs(v[0]) >= L/2.0) distance_squared += pow(abs(v[0]) - L/2.0, 2.0);


        #ifndef VEC2D
        distance_squared += pow(v[2],2);
        #endif

        if(distance_squared <= Rminsq) break;
    };
    return v;
};

Vec AtomGroup::com() const{
    Vec v = Vec::Zero();
    for(unsigned int i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        flt curmass = (*this)[i].m;
        v += (*this)[i].x * curmass;
    }
    return v / mass();
};

flt AtomGroup::mass() const{
    flt m = 0;
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        m += (*this)[i].m;
    }
    return m;
};

Vec AtomGroup::com_velocity() const {
    return momentum() / mass();
};

Vec AtomGroup::momentum() const{
    Vec tot = Vec::Zero();
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        flt curmass = a.m;
        tot += a.v * curmass;
    }
    return tot;
};

Vec AtomGroup::com_force() const{
    Vec tot = Vec::Zero();
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        tot += a.f;
    }
    return tot;
};

flt AtomGroup::gyradius() const{
    Vec avgr = Vec::Zero();
    for(uint i = 0; i<size(); i++){
        avgr += (*this)[i].x;
    }
    avgr /= size(); // now avgr is the average location, akin to c.o.m.
    flt Rgsq = 0;
    for(uint i = 0; i<size(); i++){
        Rgsq += ((*this)[i].x - avgr).squaredNorm();
    }

    return sqrt(Rgsq/size());
};

#ifdef VEC3D
Vec AtomGroup::angular_momentum(const Vec loc) const{
    Vec tot = Vec::Zero();
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        if(curmass <= 0 or isinf(curmass)) continue;
        Vec newloc = (*this)[i].x - loc;
        tot += cross(newloc, (*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};

Vec AtomGroup::torque(const Vec loc) const{
    Vec tot = Vec::Zero();
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        if(curmass <= 0 or isinf(curmass)) continue;
        Vec newloc = (*this)[i].x - loc;
        tot += cross(newloc, (*this)[i].f); // r x f
    }
    return tot;
};
#elif defined VEC2D
flt AtomGroup::angular_momentum(const Vec loc) const{
    flt tot = 0;
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        Vec newloc = a.x - loc;
        tot += cross(newloc, a.v) * a.m; // r x v m = r x p
    }
    return tot;
};

flt AtomGroup::torque(const Vec loc) const{
    flt tot = 0.;
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        if(curmass <= 0 or isinf(curmass)) continue;
        Vec newloc = (*this)[i].x - loc;
        tot += cross(newloc, (*this)[i].f); // r x f
    }
    return tot;
};
#endif

#ifdef VEC3D
flt AtomGroup::moment_about(const Vec loc, const Vec axis) const{
    if (axis.squaredNorm() == 0) return 0;
    flt tot = 0;
    Vec newloc = Vec::Zero();
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        newloc = perpto(a.x - loc, axis);
        tot += newloc.dot(newloc) * a.m;
    }
    return tot;
};

Matrix AtomGroup::moment(const Vec loc) const{
    Matrix I = Matrix::Zero();
    for(uint i=0; i<size(); i++){
        flt curmass = (*this)[i].m;
        if(curmass <= 0 or isinf(curmass)) continue;
        Vec r = (*this)[i].x - loc;
        flt x = r(0), y = r(1), z = r(2);
        I(0,0) += curmass * (y*y + z*z);
        I(1,1) += curmass * (x*x + z*z);
        I(2,2) += curmass * (x*x + y*y);
        I(0,1) -= curmass * (x*y);
        I(0,2) -= curmass * (x*z);
        I(1,2) -= curmass * (y*z);
    }
    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);
    return I;
};

Vec AtomGroup::omega(const Vec loc) const{
    Matrix Inv = moment(loc).inverse();
    Vec L = angular_momentum(loc);
    return (Inv * L).transpose();
};

void AtomGroup::add_omega(Vec w, Vec loc){
    for(uint i=0; i<size(); i++){
        Vec r = (*this)[i].x - loc;
        (*this)[i].v -= r.cross(w);
    }
};
#elif defined VEC2D
flt AtomGroup::moment(const Vec loc) const{
    flt tot = 0;
    Vec newloc = Vec::Zero();
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        newloc = a.x - loc;
        tot += newloc.dot(newloc) * a.m;
    }
    return tot;
};

void AtomGroup::add_omega(flt w, Vec loc){
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m <= 0 or isinf(a.m)) continue;
        Vec r = a.x - loc;
        a.v -= perp(r).normalized()*w;
    }
};
#endif

flt AtomGroup::kinetic_energy(const Vec originvelocity) const{
    flt totE = 0;
    Vec curv = Vec::Zero();
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m == 0 or isinf(a.m)) continue;
        curv = a.v - originvelocity;
        totE += a.m/2 * curv.dot(curv);
    }
    return totE;
};

void AtomGroup::add_velocity(Vec v){
    for(uint i=0; i<size(); i++){
        (*this)[i].v += v;
    }
};

void AtomGroup::randomize_velocities(flt T){
    for(uint i=0; i<size(); i++){
        Atom& a = (*this)[i];
        if(a.m == 0 or isinf(a.m)) continue;
        a.v = rand_vec() * sqrt(T/a.m);
    }
};

void AtomGroup::reset_forces(){
    for(uint i=0; i<size(); i++){
        (*this)[i].f = Vec::Zero();
    }
};

//~ void AtomGroup::vverlet1(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].x += (*this)[i].v * dt + (*this)[i].a * (dt*dt/2);
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

//~ void AtomGroup::vverlet2(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

void LeesEdwardsBox::shear(flt dgamma, AtomGroup &atoms){
    // TODO: is this right for non-square boxes?
    // Is gamma Δx / Lx or Δx / Ly?
    for(uint i=0; i<atoms.size(); i++){
        //flt dy = remainder(atoms[i].x[1], boxsize[1]);
        flt dy = atoms[i].x[1];
        atoms[i].x[0] += dy * dgamma;
    }

    gamma += dgamma;
};
