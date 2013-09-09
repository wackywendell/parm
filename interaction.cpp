#include "interaction.hpp"

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
Vec atomgroup::angmomentum(const Vec &loc, Box *box) const{
    Vec tot = Vec();
    flt curmass;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box->diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};
#elif defined VEC2D
flt atomgroup::angmomentum(const Vec &loc, Box *box) const{
    flt tot = 0;
    flt curmass;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box->diff((*this)[i].x, loc);
        tot += newloc.cross((*this)[i].v) * curmass; // r x v m = r x p
    }
    return tot;
};
#endif

#ifdef VEC3D
flt atomgroup::moment(const Vec &loc, const Vec &axis, Box *box) const{
    if (axis.sq() == 0) return 0;
    flt curmass;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box->diff((*this)[i].x, loc).perpto(axis);
        tot += newloc.dot(newloc) * curmass;
    }
    return tot;
};

Matrix<flt> atomgroup::moment(const Vec &loc, Box *box) const{
    flt curmass;
    Matrix<flt> I;
    Vec r;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        r = box->diff((*this)[i].x, loc);
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

Vec atomgroup::omega(const Vec &loc, Box *box) const{
    Matrix<flt> Inv = moment(loc, box).SymmetricInverse();
    return Inv * (angmomentum(loc, box));
};

void atomgroup::addOmega(Vec w, Vec loc, Box *box){
    for(uint i=0; i<size(); i++){
        Vec r = box->diff((*this)[i].x, loc);
        (*this)[i].v -= r.cross(w);
    }
};
#elif defined VEC2D
flt atomgroup::moment(const Vec &loc, Box *box) const{
    flt curmass;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<size(); i++){
        curmass = this->getmass(i);
        newloc = box->diff((*this)[i].x, loc);
        tot += newloc.dot(newloc) * curmass;
    }
    return tot;
};

void atomgroup::addOmega(flt w, Vec loc, Box *box){
    for(uint i=0; i<size(); i++){
        Vec r = box->diff((*this)[i].x, loc);
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
    uint n = a - atoms;
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

atomid metagroup::get_id(atom* a){
    for(uint n=0; n<atoms.size(); n++)
        if(atoms[n] == a) return atomid(a, n);
    return atomid();
};


//~ flt LJforce::energy(flt r){
    //~ flt l = r / sigma;
    //~ return 4*epsilon*(pow(l,-12) - pow(l,-6));
//~ }
//~ 
//~ flt LJforce::forces(flt r){
    //~ flt l = r / sigma;
    //~ return 4*epsilon*(12*pow(l,-13)-6*pow(l,-7))/sigma;
//~ }

//~ LJcutoff::LJcutoff(const flt ep, const flt sig, const flt cut) : 
        //~ LJforce(ep, sig){
    //~ setcut(cut);
//~ }
//~ 
//~ void LJcutoff::setcut(const flt cut){
    //~ cutoff = cut;
    //~ flt l = cutoff / sigma;
    //~ cutoffenergy = 4*epsilon*(pow(l,-12) - pow(l,-6));
//~ }
//~ 
//~ flt LJcutoff::energy(const flt r){
    //~ if(r > cutoff) return 0;
    //~ flt l = r / sigma;
    //~ return 4*epsilon*(pow(l,-12) - pow(l,-6)) - cutoffenergy;
//~ }
//~ 
//~ flt LJcutoff::forces(const flt r){
    //~ if(r > cutoff) return 0;
    //~ flt l = r / sigma;
    //~ return 4*epsilon*(12*pow(l,-13)-6*pow(l,-7))/sigma;
//~ }
//~ 
//~ LJcutrepulsive::LJcutrepulsive(const flt ep, const flt sig, const flt cut) : 
        //~ LJforce(ep, sig){
    //~ setcut(cut);
//~ }
//~ 
//~ void LJcutrepulsive::setcut(const flt cut){
    //~ cutoff = cut;
    //~ flt l = cutoff / sigma;
    //~ cutoffenergy = 4*(pow(l,-12));
//~ }
//~ 
//~ flt LJcutrepulsive::energy(const flt r, const flt eps, const flt cutoff){
    //~ if(cutoff>0 and r > cutoff) return 0;
    //~ return 4*eps*(pow(r,-12)) - energy(cutoff,eps,0);
//~ }
//~ 
//~ flt LJcutrepulsive::forces(const flt r, const flt sig, const flt eps, const flt cutoff){
    //~ if(r > cutoff) return 0;
    //~ return 4*eps*(12*pow(r,-13))/sig;
//~ }

flt spring::energy(const Vec r){
    flt m = r.mag();
    flt l = m - x0;
    return .5 * springk * l*l;
}

Vec spring::forces(const Vec r){
    flt m = r.mag();
    flt fmag = (x0 - m) * springk;
    return r * (fmag / m);
}


electricScreened::electricScreened(const flt screenLength, const flt q1, 
            const flt q2, const flt cutoff) : screen(screenLength), 
            q1(q1), q2(q2), cutoff(cutoff){};

flt electricScreened::energy(const flt r, const flt qaqb, const flt screen, const flt cutoff){
    if(cutoff <= 0) return exp(-r/screen) * qaqb / r;
    if(r > cutoff) return 0;
    return exp(-r/screen) * qaqb / r - energy(r, qaqb, screen, 0);
};

Vec electricScreened::forces(const Vec r, const flt qaqb, const flt screen, const flt cutoff){
    flt d = r.mag();
    if(cutoff > 0 and d > cutoff) return Vec();
    flt fmag = exp(-d/screen) * (qaqb/d) * (1/d + 1/screen);
    return r * (fmag/d);
};

flt bondangle::energy(const Vec& r1, const Vec& r2){
    flt costheta = r1.dot(r2) / r1.mag() / r2.mag();
    if(!usecos) return springk*pow(acos(costheta) - theta0,2)/2;
    else return springk*pow(costheta - cos(theta0),2)/2;
}

Nvector<Vec, 3> bondangle::forces(const Vec& r1, const Vec& r2){
    flt r1mag = r1.mag();
    flt r2mag = r2.mag();
    
    flt costheta = r1.dot(r2) / r1mag / r2mag;
    flt theta = acos(costheta);
    //theta is now the angle between x1 and x2
    
    flt fmag;
    if(usecos) fmag = -springk*(cos(theta0) - costheta)*sin(theta);
    else fmag = springk*(theta0 - theta); // torque magnitude
    // We have V = \frac{1}{2}k(\theta-\theta_{0})^{2}
    // Then -f = grad V = \frac{k}{r}(\theta-\theta_{0})\hat{\theta}
    // first we get the direction:
    Nvector<Vec, 3> force;
    force[0] = r2.perpto(r1);
    force[0].normalize();
    force[2] = r1.perpto(r2);
    force[2].normalize();
    
    // now we get magnitude: 
    force[0] *= fmag/r1mag;
    force[2] *= fmag/r2mag;
    
    //~ cout << force[2] << x2 << "force(2).x2: " << force[2].dot(x2) << endl;
    force[1] = -(force[0] + force[2]);
    // The direction of the force on the first atom (f0) is 
    // perpendicular to x1, and same for f2.
    // **TODO** its possible that x1 = +/-x2, and then x1.perp(x2) = 0
    // and then we get a divide by zero error.
    
    return force;
}

#ifdef VEC3D
dihedral::dihedral(const vector<flt> cvals, const vector<flt> svals, bool usepow) : 
                    coscoeffs(cvals), sincoeffs(svals), usepow(usepow){
}

Nvector<Vec,4> dihedral::forces(const Vec &r1, const Vec &r2, 
                   const Vec &r3) const {
    // Taken from Rapaport "Art of Molecular Dynamics Simulation" p.279
    // The expressions and notation are very close to that of the book.
    
    // Note that Rappaport defines it as such:
    /*
     The dihedral angle is defined as the angle between the
planes formed by atoms 1,2,3 and 2,3,4 measured in the plane normal to the 2–3
bond; it is zero when all four atoms are coplanar and atoms 1 and 4 are on opposite
sides of the bond. */
    
    // so we need a negative sign.
    
    // ri corresponds to Rapaport's b_i
    
    flt c[3][3];
    c[0][0] = r1.dot(r1);
    c[0][1] = c[1][0] = r1.dot(r2);
    c[0][2] = c[2][0] = r1.dot(r3);
    c[1][1] = r2.dot(r2);
    c[1][2] = c[2][1] = r2.dot(r3);
    c[2][2] = r3.dot(r3);
    
    flt p = c[0][2] * c[1][1] - c[0][1] * c[1][2];
    flt qa = c[0][0]*c[1][1] - c[0][1] * c[0][1];
    flt qb = c[1][1]*c[2][2] - c[1][2] * c[1][2];
    flt q = qa * qb;
    flt sqq = sqrt(q);
    
    flt t1 = p;
    flt t2 = c[0][0] * c[1][2] - c[0][1] * c[0][2];
    flt t3 = c[0][1] * c[0][1] - c[0][0] * c[1][1];
    flt t4 = c[1][1] * c[2][2] - c[1][2] * c[1][2];
    flt t5 = c[0][2] * c[1][2] - c[0][1] * c[2][2];
    flt t6 = -p;
    
    Nvector<Vec, 4> derivs;
    
    /*
     Rapaport: The dihedral angle is defined as the angle between the
     planes formed by atoms 1,2,3 and 2,3,4 measured in the plane 
     normal to the 2–3 bond; it is zero when all four atoms are 
     coplanar and atoms 1 and 4 are on opposite sides of the bond. 
     
     Note that this is the *opposite* of the chemical definition.
     
     flt const0 = c[1][1]/(sqq * qa);
     flt const3 = c[1][1]/(sqq * qb);
     
     We add a negative in at the beginning of those two to give us the 
     chemical definition.
     */
    
    flt const0 = -c[1][1]/(sqq * qa);
    flt const3 = -c[1][1]/(sqq * qb);
    
    
    derivs[0] = (r1 * t1 + r2 * t2 + r3 * t3) * const0;
    derivs[3] = (r1 * t4 + r2 * t5 + r3 * t6) * const3;
    
    derivs[1] = derivs[0] * (-1 - c[0][1]/c[1][1]) +
                            derivs[3] * (c[1][2]/c[1][1]);
    derivs[2] = derivs[0] * (c[0][1]/c[1][1]) -
                            derivs[3] * (1 + c[1][2]/c[1][1]);
    
     /*
     Rapaport says costheta = p/sqrt(q); we add a negative for the cosine.
     */
    
    flt dcostheta;
    if(sincoeffs.empty() and !usepow){        
        flt costheta = -p/sqq;
        // costheta =-1 corresponds to atoms 1 and 4 on opposite sides of the bond (zigzag)
        // costheta = 1 corresponds to a C shape

        dcostheta = dudcosthetaCOS(costheta); // F = -dU/d(costheta)
        
    
        //~ if(abs(costheta) < .7)
            //~ cout << "forces cos: " << costheta << " getcos: " << getcos(r1,r2,r3)
                 //~ << " dcostheta: " << dcostheta << '\n';
    } else {
        dcostheta = dudcostheta(getang(r1, r2, r3));
    }
        
    derivs *= -dcostheta;  // F = -dU/d(costheta)
    //~ assert(derivs[0].sq() < 1e8);
    //~ assert(derivs[1].sq() < 1e8);
    //~ assert(derivs[2].sq() < 1e8);
    //~ assert(derivs[3].sq() < 1e8);
        
    
    //~ flt mag = sqrt(derivs[0].sq() +derivs[1].sq() + derivs[2].sq() +
                    //~ derivs[3].sq());
    //~ 
    //~ std::cout << "costheta:" << costheta << " dcos:" << dcostheta
              //~ << " derivs:" << derivs  << " : " << mag << std::endl;
    return derivs;
    
    // pea79, dun92
    // Pear, M. R. and Weiner, J. H., Brownian dynamics study of a polymer chain of linked rigid bodies, J. Chem. Phys. 71 (1979) 212.
    // Dunn, J. H., Lambrakos, S. G., Moore, P. G., and Nagumo, M., An algorithm for calculating intramolecular angle-dependent forces on vector computers, J. Comp. Phys. 100 (1992) 17.


}

flt dihedral::dudcosthetaCOS(const flt costheta) const{
    assert(sincoeffs.empty());
    assert(!usepow);
    flt tot = 0;
    unsigned int cosmx = coscoeffs.size();
    for(unsigned int i=1; i < cosmx; i++){
        tot += coscoeffs[i] * i * pow(costheta, flt(i-1));
    }
    //~ cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    return tot;
}

flt dihedral::dudcostheta(const flt theta) const{
    flt tot = 0;
    unsigned int cosmx = coscoeffs.size();
    unsigned int sinmx = sincoeffs.size();
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;
    if(usepow) {
        flt costheta = cos(theta), sintheta = sin(theta);
        flt cottheta = -costheta / sintheta;
        for(unsigned int i=1; i < mx; i++){
            if (i < cosmx) tot += coscoeffs[i] * i * pow(costheta, flt(i-1));
            if (i < sinmx) tot += sincoeffs[i] * i * cottheta
                                    * pow(sintheta, flt(i-1));
        }
    } else {
        flt csctheta = 1 / sin(theta);
        for(unsigned int i=1; i < mx; i++){
            flt cositheta = cos(i*theta), sinitheta = sin(i*theta);
            if (i < cosmx) tot += coscoeffs[i] * csctheta * i * sinitheta;
            if (i < sinmx) tot -= sincoeffs[i] * csctheta * i * cositheta;
        }
    }
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    return tot;
}

flt dihedral::getcos(const Vec &r1, const Vec &r2, 
                   const Vec &r3){
   // The two normals to the planes
    Vec n1 = r1.cross(r2);
    Vec n2 = r2.cross(r3);
    //~ cout << r1 << ',' << r2  << ',' << r3 << "\n";
    flt n1mag = n1.mag();
    flt n2mag = n2.mag();
    
    if (n1mag == 0 or n2mag == 0) return -100; 
    // if one plane is ill-defined, then we have no torsion angle

    return (n1.dot(n2) / n1mag / n2mag);
};

flt dihedral::getang(const Vec &r1, const Vec &r2, 
                   const Vec &r3){
    
    return atan2(r1.dot(r2.cross(r3))*r2.mag(), (r1.cross(r2).dot(r2.cross(r3))));
};

flt dihedral::energy(const flt ang) const{
    
    flt costheta = (usepow ? cos(ang) : NAN);
    flt sintheta = (usepow ? sin(ang) : NAN);
    
    unsigned int cosmx = coscoeffs.size();
    unsigned int sinmx = sincoeffs.size();
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;
    
    flt tot = 0;
    for(unsigned int i=0; i < mx; i++){
        if(usepow) {
            if(i < cosmx) tot += coscoeffs[i] * pow(costheta, flt(i));
            if(i < sinmx) tot += sincoeffs[i] * pow(sintheta, flt(i));
        } else {
            if(i < cosmx) tot += coscoeffs[i] * cos(i * ang);
            if(i < sinmx) tot += sincoeffs[i] * sin(i * ang);
        }
    }
    
    return tot;
}
#endif
//~ flt interactgroup::energy(Box *box){
    //~ flt E=0;
    //~ vector<interaction*>::iterator it;
    //~ for(it = inters.begin(); it < inters.end(); it++){
        //~ E += (*it)->energy(box);
    //~ }
    //~ return E;
//~ };

//~ void interactgroup::setForces(Box *box){
    //~ vector<interaction*>::iterator it;
    //~ for(it = inters.begin(); it < inters.end(); it++){
        //~ (*it)->setForces(box);
    //~ }
//~ };

bondpairs::bondpairs(vector<bondgrouping> pairs) : pairs(pairs){};

flt bondpairs::energy(Box *box){
    flt E=0;
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = diff(it->a1->x, it->a2->x);
        E += spring(it->k, it->x0).energy(r);
    }
    return E;
}

void bondpairs::setForces(Box *box){
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        Vec r = diff(atom1.x, atom2.x);
        Vec f = spring(it->k, it->x0).forces(r);
        //~ assert(f.sq() < 10000000);
        atom1.f += f;
        atom2.f -= f;
    }
}

flt bondpairs::pressure(Box *box){
    return 0;
    vector<bondgrouping>::iterator it;
    flt P=0;
    for(it = pairs.begin(); it < pairs.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        Vec r = diff(atom1.x, atom2.x);
        Vec f = spring(it->k, it->x0).forces(r);
        P += f.dot(r);
    }
    return P;
}

flt bondpairs::mean_dists() const{
    flt dist=0;
    uint N=0;
    vector<bondgrouping>::const_iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = diff(it->a1->x, it->a2->x);
        dist += abs(r.mag() - it->x0);
        N++;
    }
    return dist/N;
}

flt bondpairs::std_dists() const{
    flt stds=0;
    uint N=0;
    vector<bondgrouping>::const_iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = diff(it->a1->x, it->a2->x);
        flt curdist = r.mag() - it->x0;
        stds += curdist*curdist;
        N++;
    }
    return sqrt(stds/N);
}

angletriples::angletriples(vector<anglegrouping> triples) : triples(triples){};

flt angletriples::energy(Box *box){
    flt E=0;
    vector<anglegrouping>::iterator it;
    for(it = triples.begin(); it < triples.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom2.x, atom3.x);
        E += bondangle(it->k, it->x0).energy(r1,r2);
    }
    return E;
}

void angletriples::setForces(Box *box){
    vector<anglegrouping>::iterator it;
    for(it = triples.begin(); it < triples.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom2.x, atom3.x);
        Nvector<Vec,3> f = bondangle(it->k, it->x0).forces(r1, r2);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        //~ assert(f[0].sq() < 1000000);
        //~ assert(f[1].sq() < 1000000);
        //~ assert(f[2].sq() < 1000000);
    }
};

flt angletriples::mean_dists() const{
    flt dist=0;
    uint N=0;
    vector<anglegrouping>::const_iterator it;
    //~ cout << "angle diffs: ";
    cout.precision(3);
    for(it = triples.begin(); it < triples.end(); it++){
        Vec r1 = diff(it->a1->x, it->a2->x);
        Vec r2 = diff(it->a3->x, it->a2->x);
        flt theta = acos(r1.dot(r2) / r1.mag() / r2.mag());
        //~ flt curdist = abs(theta - it->x0);
        dist += abs(theta - it->x0);
        //~ if(curdist > .2){
            //~ flt E = bondangle(it->k, it->x0).energy(r1,r2);
            //~ Vec f0 = bondangle(it->k, it->x0).forces(r1,r2)[0];
            //~ cout << "(" << theta << "," << it->x0 << "," << abs(theta - it->x0)
                 //~ << ";" << it->k << "," << E
                 //~ << ",f:" << f0.mag() <<  ")" << ", ";
         //~ }
        N++;
    }
    //~ cout << "Total: " << N << '\n';
    return dist/N;
};

flt angletriples::std_dists() const{
    flt stds=0;
    uint N=0;
    vector<anglegrouping>::const_iterator it;
    for(it = triples.begin(); it < triples.end(); it++){
        Vec r1 = diff(it->a1->x, it->a2->x);
        Vec r2 = diff(it->a3->x, it->a2->x);
        flt theta = acos(r1.dot(r2) / r1.mag() / r2.mag());
        flt curdist = theta - (it->x0);
        //~ cout << "angles: " << theta << " x0: " << it->x0
             //~ << " cur: " << curdist << "\n";
        stds += curdist*curdist;
        N++;
    }
    return sqrt(stds/N);
};

#ifdef VEC3D
dihedrals::dihedrals(vector<dihedralgrouping> ds) : groups(ds){};

flt dihedrals::energy(Box *box){
    flt E=0;
    vector<dihedralgrouping>::iterator it;
    for(it = groups.begin(); it < groups.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        atom & atom4 = *it->a4;
        Vec r1 = dihedralgrouping::diff(atom2.x, atom1.x);
        Vec r2 = dihedralgrouping::diff(atom3.x, atom2.x);
        Vec r3 = dihedralgrouping::diff(atom4.x, atom3.x);
        E += it->dih.energy(r1,r2,r3);
    }
    return E;
}

void dihedrals::setForces(Box *box){
    vector<dihedralgrouping>::iterator it;
    for(it = groups.begin(); it < groups.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        atom & atom4 = *it->a4;
        Vec r1 = dihedralgrouping::diff(atom2.x, atom1.x);
        Vec r2 = dihedralgrouping::diff(atom3.x, atom2.x);
        Vec r3 = dihedralgrouping::diff(atom4.x, atom3.x);
        Nvector<Vec,4> f = it->dih.forces(r1, r2, r3);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        atom4.f += f[3];
        //~ flt maxf = 1000000;
        //~ if(f[0].sq() > maxf or f[1].sq() > maxf or f[2].sq() > maxf 
            //~ or f[3].sq() > maxf){
                //~ cout << "dihedral overload: " << r1 << r2 << r3 << " :: " <<
                //~ f[0] << f[1] << f[2] << f[3] << "\n";
                //~ cout << "dihedral overload energy: " << dihedral(it->nums).energy(r1,r2,r3) << "\n";
            //~ }
    }
};

flt dihedrals::mean_dists() const{
    flt dist=0;
    uint N=0;
    vector<dihedralgrouping>::const_iterator it;
    for(it = groups.begin(); it < groups.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        atom & atom4 = *it->a4;
        Vec r1 = dihedralgrouping::diff(atom1.x, atom2.x);
        Vec r2 = dihedralgrouping::diff(atom2.x, atom3.x);
        Vec r3 = dihedralgrouping::diff(atom3.x, atom4.x);
        flt cosine = dihedral::getcos(r1, r2, r3);
        dist += cosine;
        N++;
    }
    return dist/N;
};
//~ 
//~ flt dihedrals::std_dists() const{
    //~ flt stds=0;
    //~ uint N=0;
    //~ vector<anglegrouping>::const_iterator it;
    //~ for(it = triples.begin(); it < triples.end(); it++){
        //~ Vec r1 = diff(it->a1->x, it->a2->x);
        //~ Vec r2 = diff(it->a3->x, it->a2->x);
        //~ flt theta = acos(r1.dot(r2) / r1.mag() / r2.mag());
        //~ flt curdist = theta - (it->x0);
        // cout << "angles: " << theta << " x0: " << it->x0
             // << " cur: " << curdist << "\n";
        //~ stds += curdist*curdist;
        //~ N++;
    //~ }
    //~ return sqrt(stds/N);
//~ };
#endif

void pairlist::clear(){
    map<const atomid, set<atomid> >::iterator mapiter;
    for(mapiter = pairs.begin(); mapiter != pairs.end(); mapiter++){
        mapiter->second.clear();
    }
};

neighborlist::neighborlist(Box *box, const flt innerradius, const flt outerradius) :
                box(box), critdist(innerradius), skinradius(outerradius),
                atoms(){};

neighborlist::neighborlist(Box *box, atomgroup &group, const flt innerradius,
            const flt outerradius, pairlist ignore) :
                box(box), critdist(innerradius), skinradius(outerradius),
                atoms(), ignorepairs(ignore), ignorechanged(false){
    lastlocs.resize(group.size());
    for(uint i=0; i<group.size(); i++){
        atomid a = group.get_id(i);
        ignorepairs.ensure(a);
        atoms.add(a.pointer());
        lastlocs[i] = a.x();
    }
    update_list(true);
    updatenum = 1;
};

atomid neighborlist::get_id(atom* a){
    for(uint i=0; i < atoms.size(); i++)
        if (atoms.get(i) == a) return atoms.get_id(i);
    return atomid();
};

void neighborlist::ignore(atom* a, atom* b){
    atomid aid = get_id(a), bid = get_id(b);
    assert(a != NULL);
    assert(b != NULL);
    ignore(aid, bid);
};

bool neighborlist::update_list(bool force){
    flt curdist = 0, bigdist = 0, biggestdist = 0;
    // biggestdist is the distance the furthest-moving atom has gone
    // bigdist is the next furthest
    
    if(not force and not ignorechanged){ // check if we need to update
        for(uint i=0; i < atoms.size(); i++){
            atom &atm = atoms[i];
            curdist = (atm.x - lastlocs[i]).mag();
            if(curdist > biggestdist){
                bigdist = biggestdist;
                biggestdist = curdist;
            }
            else if(curdist > bigdist){
                bigdist = curdist;
            }
            else continue; // this one hasn't moved enough to worry about
            
            // if we haven't continued, that means bigdist and/or biggestdist
            // have been changed, so we check
            if(bigdist + biggestdist >= skinradius - critdist){
                force = true;
                break; // we do need to update, stop checking
            }
        }
        // if we haven't found anything, then we're done; no need to update.
        if(not force){
            //~ if(skinradius - critdist > 5)
                //~ cout << "Not updating:" << bigdist + biggestdist << '-' 
                //~ << skinradius - critdist << '\n';
            return false;
        }
        
        //~ cout << "Updating:" << bigdist + biggestdist << '-' 
            //~ << skinradius - critdist << '\n';
    }
    
    // time to update
    updatenum++;
    ignorechanged = false;
    curpairs.clear();
    for(uint i=0; i<atoms.size(); i++){
        atomid a1=atoms.get_id(i);
        lastlocs[i] = a1.x();
        for(uint j=0; j<i; j++){
            atomid a2=atoms.get_id(j);
            if (ignorepairs.has_pair(a1, a2)) continue;
            if(box->diff(a1.x(), a2.x()).mag() < skinradius)
                curpairs.push_back(idpair(a1, a2));
        }
    }
    //~ cout << "neighborlist::update_list:: done.\n";
    // print stuff about the current update
    //~ set<atomid> curset = (ignorepairs.get_pairs(atoms.back()));
    //~ cout << "neighborlist | atoms: " << atoms.size() <<  "pairs: " << curpairs.size() << " ignored -1: "
         //~ << curset.size() << "\n";
    //~ cout << "ignored -1:";
    //~ for(set<atomid>::iterator it=curset.begin(); it!=curset.end(); it++)
        //~ cout << " " << it->n();
    //~ cout << '\n' << "paired:";
    //~ for(vector<idpair>::iterator it=begin(); it!=end(); it++)
        //~ if(it->first() == atoms.back()) cout << " " << it->last().n();
    //~ cout << '\n';
    
    return true;
}

ContactTracker::ContactTracker(Box *box, atomgroup *atoms, vector<flt> dists) :
    atoms(atoms), dists(dists), contacts(), breaks(0), formations(0),
        incontact(0){
    //~ cout << "Making contact tracker." << endl;
    uint N = atoms->size();
    dists.resize(N);
    contacts.resize(N);
    for(uint i=0; i<N; i++){
        contacts[i].resize(i, false);
    }
    update(box);
    breaks = 0;
    formations = 0;
};

void ContactTracker::update(Box *box){
    //~ cout << "Contact tracker update." << endl;
    uint N = atoms->size();
    incontact = 0;
    for(uint i=0; i<N; i++){
        Vec ri = atoms->get(i)->x;
        for(uint j=0; j<i; j++){
            Vec rj = atoms->get(j)->x;
            Vec dr = box->diff(ri,rj);
            bool curcontact = (dr.mag() <= ((dists[i] + dists[j])/2));
            if(curcontact && (!contacts[i][j])) formations++;
            else if(!curcontact && contacts[i][j]) breaks++;
            if(curcontact) incontact++;
            
            contacts[i][j] = curcontact;
        }
    }
    //~ cout << "Contact tracker update done." << endl;
};

LJsimple::LJsimple(flt cutoff, vector<LJatom> atms) : atoms(atms){};

flt LJsimple::energy(Box *box){
    flt E = 0;
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec dist = box->diff(pair.atom1.x(), pair.atom2.x());
        E += LJrepulsive::energy(dist, pair.sigma, pair.epsilon);
    }
    return E;
};

void LJsimple::setForces(Box *box){
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec r = box->diff(pair.atom1.x(), pair.atom2.x());
        Vec f = LJrepulsive::forces(r, pair.sigma, pair.epsilon);
        pair.atom1.f() += f;
        pair.atom2.f() -= f;
    }
};

flt LJsimple::pressure(Box *box){
    flt P=0;
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec r = box->diff(pair.atom1.x(), pair.atom2.x());
        Vec f = LJrepulsive::forces(r, pair.sigma, pair.epsilon);
        P += r.dot(f);
    }
    return P;
};

atomid LJsimple::get_id(atom* a){
    for(vector<LJatom>::iterator it=atoms.begin(); it!=atoms.end(); it++)
        if((*it) == a) return *it;
    return atomid();
};


Charges::Charges(flt screen, flt k, vector<Charged> atms) : atoms(atms),
                screen(screen), k(k){};

atomid Charges::get_id(atom* a){
    for(vector<Charged>::iterator it=atoms.begin(); it!=atoms.end(); it++)
        if((*it) == a) return *it;
    return atomid();
};

flt Charges::energy(Box *box){
    flt E = 0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec dist = box->diff(it->x(), it2->x());
        E += k*electricScreened::energy(dist.mag(), (it->q) * (it2->q), screen);
    }
    return E;
};

void Charges::setForces(Box *box){
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec r = box->diff(it->x(), it2->x());
        Vec f = electricScreened::forces(r, (it->q) * (it2->q), screen)*k;
        it->f() += f;
        it2->f() -= f;
    }
};

flt Charges::pressure(Box *box){
    flt P=0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec r = box->diff(it->x(), it2->x());
        Vec f = electricScreened::forces(r, (it->q) * (it2->q), screen)*k;
        P += r.dot(f);
    }
    return P;
};

bool jamminglist::operator<(const jamminglist& other ){
    //return distsq < other.distsq;
    if(other.distsq  - distsq > 1e-8) return true;
    if(distsq  - other.distsq > 1e-8) return false;
    //~ cout << "\nWithin 1e-8\n";
    
    uint sz = size();
    uint osz = other.size();
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for(uint i=0; i<sz; i++){
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing sizes...\n";
    //~ cout << "Equal!\n";
    return false; // consider them equal
};

#ifdef VEC2D
/* There are two ways of looking at the different arrangements.
 * In both cases, we leave A the same as it was, and rotate / flip / translate B.
 * Also in both cases, we wrap A, then subtract off its COM (in an infinite box).
 * 
 * Method 1:
   * "unwrap" the box into 9 boxes (2D), and then choose a box for each
   * of the particles.
   * (N-1)⁹, right? (×8×N!, with rotations / flips / permutations)
 * Method 2:
   * Pick a particle (in B), move the box so that's on the left.
   * Pick another particle (maybe the same), move the box so its on the bottom.
   * Calculate COM of B without PBC, subtract that off.
   * Compare A and B, but using the PBC to calculate distances.
   * N² (×8×N!)
 * Method 3:
   * Only try different rotations, but measure distances as Σ (\vec r_{ij} - \vec s_{ij})²
 
 * We'll go with method 3.
*/

bool jamminglistrot::operator<(const jamminglistrot& other ){
    //return distsq < other.distsq;
    if(other.distsq  - distsq > 1e-8) return true;
    if(distsq  - other.distsq > 1e-8) return false;
    
    uint sz = size();
    uint osz = other.size();
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << "\nWithin 1e-8. ";
    
    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for(uint i=0; i<sz; i++){
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing rotations...";
    if(rotation < other.rotation) return true;
    if(rotation > other.rotation) return false;
    
    //~ cout << " Comparing sizes...";
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << "Equal!\n";
    return false; // consider them equal
};

jammingtree2::jammingtree2(Box *box, vector<Vec>& A0, vector<Vec>& B0)
            : box(box), jlists(), A(A0), Bs(8, B0){
    for(uint rot=0; rot < 8; rot++){
        for(uint i=0; i<B0.size(); i++){
                        Bs[rot][i] = B0[i].rotate_flip(rot); }
        if(A0.size() <= B0.size()) jlists.push_back(jamminglistrot(rot));
        //~ cout << "Created, now size " << jlists.size() << endl;
    }
    
};

flt jammingtree2::distance(jamminglistrot& jlist){
    flt dist = 0;
    uint rot = jlist.rotation;
    for(uint i=1; i<jlist.size(); i++){
        uint si = jlist.assigned[i];
        for(uint j=0; j<i; j++){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            dist += box->diff(rij, sij).sq();
        }
    }
    return dist / jlist.assigned.size();
};

list<jamminglistrot> jammingtree2::expand(jamminglistrot curjlist){
    vector<uint>& curlist = curjlist.assigned;
    list<jamminglistrot> newlists = list<jamminglistrot>();
    if(curlist.size() >= A.size()){
        return newlists;
    }
    
    uint N = Bs[curjlist.rotation].size();
    for(uint i=0; i < N; i++){
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;
        
        jamminglistrot newjlist = jamminglistrot(curjlist, i, 0);
        newjlist.distsq = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool jammingtree2::expand(){
    jamminglistrot curjlist = jlists.front();
    list<jamminglistrot> newlists = expand(curjlist);
    
    if(newlists.size() <= 0){
        //~ cout << "No lists made\n";
        return false;
    }
    //~ cout << "Have " << newlists.size() << "\n";
    newlists.sort();
    //~ cout << "Sorted.\n";
    jlists.pop_front();
    //~ cout << "Popped.\n";
    jlists.merge(newlists);
    //~ cout << "Merged to size " << jlists.size() << "best dist now " << jlists.front().distsq << "\n";
    return true;
};


jammingtreeBD::jammingtreeBD(Box *box, vector<Vec>& A, vector<Vec>& B, 
                    uint cutoffA, uint cutoffB) :
            jammingtree2(box, A, B), cutoff1(cutoffA), cutoff2(cutoffB){
    if(cutoffA > cutoffB){jlists.clear();}
    if(A.size() - cutoffA > B.size() - cutoffB){jlists.clear();}
};

list<jamminglistrot> jammingtreeBD::expand(jamminglistrot curjlist){
    vector<uint>& curlist = curjlist.assigned;
    list<jamminglistrot> newlists = list<jamminglistrot>();
    if(curlist.size() >= A.size()){
        return newlists;
    }
    
    uint N = Bs[curjlist.rotation].size();
    uint start = 0, end=cutoff2;
    if(curlist.size() >= cutoff1){
        start=cutoff2;
        end = N;
    }
    for(uint i=start; i < end; i++){
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;
        
        jamminglistrot newjlist = jamminglistrot(curjlist, i, 0);
        newjlist.distsq = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool jammingtreeBD::expand(){
    jamminglistrot curjlist = jlists.front();
    list<jamminglistrot> newlists = expand(curjlist);
    
    if(newlists.size() <= 0) return false;
    newlists.sort();
    jlists.pop_front();
    jlists.merge(newlists);
    return true;
};


vector<Vec> jammingtree2::locationsB(jamminglistrot jlist){
    uint rot = jlist.rotation;
    vector<Vec> locs = vector<Vec>(jlist.size());
    
    uint N = jlist.size();
    for(uint i=0; i<N; i++){
        uint si = jlist.assigned[i];
        locs[i] = A[i];
        for(uint j=0; j<N; j++){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            locs[i] -= box->diff(rij, sij)/N;
        }
    }
    return locs;
};


vector<Vec> jammingtree2::locationsA(jamminglistrot jlist){
    uint rot = jlist.rotation;
    vector<Vec> locs = vector<Vec>(Bs[rot].size(), Vec(NAN,NAN));
    
    uint N = jlist.size();
    for(uint i=0; i<N; i++){
        uint si = jlist.assigned[i];
        locs[si] = Bs[rot][si];
        for(uint j=0; j<N; j++){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            locs[si] += box->diff(rij, sij)/N;
        }
        
        // this is an inverse rotateflip
        locs[si] = locs[si].rotate_flip_inv(rot);
    }
    return locs;
};

Vec jammingtree2::straight_diff(Box *bx, vector<Vec>& As, vector<Vec>& Bs){
    uint N = As.size();
    if(Bs.size() != N) return Vec(NAN,NAN);
    
    Vec loc = Vec();
    for(uint i=0; i<N; i++){
        for(uint j=0; j<N; j++){
            Vec rij = bx->diff(As[i], As[j]);
            Vec sij = bx->diff(Bs[i], Bs[j]);
            loc += bx->diff(rij, sij);
        }
    }
    return loc / N;
};

flt jammingtree2::straight_distsq(Box *bx, vector<Vec>& As, vector<Vec>& Bs){
    uint N = As.size();
    if(Bs.size() != N) return NAN;
    
    flt dist = 0;
    for(uint i=0; i<N; i++){
        for(uint j=0; j<N; j++){
            Vec rij = bx->diff(As[i], As[j]);
            Vec sij = bx->diff(Bs[i], Bs[j]);
            dist += bx->diff(rij, sij).sq();
        }
    }
    return dist / N;
};

#endif
