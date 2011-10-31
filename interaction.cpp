#include "interaction.hpp"

Vec atomgroup::com() const{
    flt curmass = 0;
    Vec v = Vec();
    for(unsigned int i=0; i<N(); i++){
        curmass = this->getmass(i);
        v += (*this)[i].x * curmass;
    }
    return v / mass();
};

flt atomgroup::mass() const{
    static flt m = -1;
    if(m > 0) return m;
    for(uint i=0; i<N(); i++){
        m += getmass(i);
    }
    return m;
};

Vec atomgroup::comvel() const {
    return momentum() / mass();
};

Vec atomgroup::momentum() const{
    flt curmass;
    Vec tot = Vec();
    for(uint i=0; i<N(); i++){
        curmass = this->getmass(i);
        tot += (*this)[i].v * curmass;
    }
    return tot;
};

Vec atomgroup::angmomentum(const Vec &loc) const{
    flt curmass;
    Vec tot = Vec();
    Vec newloc;
    for(uint i=0; i<N(); i++){
        curmass = this->getmass(i);
        newloc = (*this)[i].x - loc;
        tot += newloc.cross((*this)[i].v) * curmass;
    }
    return tot;
};

flt atomgroup::mominertia(const Vec &loc, const Vec &axis) const{
    flt curmass;
    flt tot = 0;
    Vec newloc;
    for(uint i=0; i<N(); i++){
        curmass = this->getmass(i);
        newloc = diff((*this)[i].x, loc).perp(axis);
        tot += newloc.dot(newloc) * curmass;
    }
    return tot;
};

flt atomgroup::kinetic(const Vec &originvelocity) const{
    flt curmass;
    flt totE = 0;
    Vec curv;
    for(uint i=0; i<N(); i++){
        curmass = getmass(i);
        curv = (*this)[i].v - originvelocity;
        totE += curmass/2 * curv.dot(curv);
    }
    return totE;
};

void atomgroup::resetForces(){
    for(uint i=0; i<N(); i++){
        (*this)[i].f = Vec(0,0,0);
    }
};

//~ void atomgroup::vverlet1(const flt dt){
    //~ for(uint i=0; i<N(); i++){
        //~ (*this)[i].x += (*this)[i].v * dt + (*this)[i].a * (dt*dt/2);
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

void atomgroup::setAccel(){
    for(uint i=0; i<N(); i++){
        (*this)[i].a = (*this)[i].f / getmass(i);
    }
};

//~ void atomgroup::vverlet2(const flt dt){
    //~ for(uint i=0; i<N(); i++){
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

LJforce::LJforce(const flt ep, const flt sig) : interactpair(),
        epsilon(ep), sigma(sig) {
}

flt LJforce::energy(flt r){
    flt l = r / sigma;
    return 4*epsilon*(pow(l,-12) - pow(l,-6));
}

flt LJforce::forces(flt r){
    flt l = r / sigma;
    return 4*epsilon*(12*pow(l,-13)-6*pow(l,-7))/sigma;
}

LJcutoff::LJcutoff(const flt ep, const flt sig, const flt cut) : 
        LJforce(ep, sig){
    setcut(cut);
}

void LJcutoff::setcut(const flt cut){
    cutoff = cut;
    flt l = cutoff / sigma;
    cutoffenergy = 4*epsilon*(pow(l,-12) - pow(l,-6));
}

flt LJcutoff::energy(const flt r){
    if(r > cutoff) return 0;
    flt l = r / sigma;
    return 4*epsilon*(pow(l,-12) - pow(l,-6)) - cutoffenergy;
}

flt LJcutoff::forces(const flt r){
    if(r > cutoff) return 0;
    flt l = r / sigma;
    return 4*epsilon*(12*pow(l,-13)-6*pow(l,-7))/sigma;
}

flt spring::energy(const Vec& r){
    flt m = r.mag();
    flt l = m - x0;
    return .5 * springk * l*l;
}

Vec spring::forces(const Vec& r){
    flt m = r.mag();
    flt fmag = (x0 - m) * springk;
    return r * (fmag / m);
}

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
    force[0] = r2.perp(r1);
    force[0].normalize();
    force[2] = r1.perp(r2);
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

dihedral::dihedral(const vector<flt> vals) : torsions(vals){
}

Nvector<Vec,4> dihedral::forces(const Vec &r1, const Vec &r2, 
                   const Vec &r3) const{
    // Taken from Rapaport "Art of Molecular Dynamics Simulation" p.279
    // The expressions and notation are very close to that of the book.
    
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
    
    flt t1 = p;
    flt t2 = c[0][0] * c[1][2] - c[0][1] * c[0][2];
    flt t3 = c[0][1] * c[0][1] - c[0][0] * c[1][1];
    flt t4 = c[1][1] * c[2][2] - c[1][2] * c[1][2];
    flt t5 = c[0][2] * c[1][2] - c[0][1] * c[2][2];
    flt t6 = -p;
    
    Nvector<Vec, 4> derivs;
    
    flt const0 = c[1][1]/(sqrt(q) * qa);
    derivs[0] = (r1 * t1 + r2 * t2 + r3 * t3) * const0;
    flt const3 = c[1][1]/(sqrt(q) * qb);
    derivs[3] = (r1 * t4 + r2 * t5 + r3 * t6) * const3;
    
    derivs[1] = derivs[0] * (-1 - c[0][1]/c[1][1]) +
                            derivs[3] * (c[1][2]/c[1][1]);
    derivs[2] = derivs[0] * (c[0][1]/c[1][1]) -
                            derivs[3] * (1 + c[1][2]/c[1][1]);
    
    //~ flt costheta = p/sqrt(q);
    flt dcostheta = dudcostheta(p/sqrt(q)); // actually -dU/d(costheta)
    
    derivs *= dcostheta;
    
    //~ flt mag = sqrt(derivs[0].sq() +derivs[1].sq() + derivs[2].sq() +
                    //~ derivs[3].sq());
    //~ 
    //~ std::cout << "costheta:" << costheta << " dcos:" << dcostheta
              //~ << " derivs:" << derivs  << " : " << mag << std::endl;
    return derivs;
}

flt dihedral::dudcostheta(const flt costheta) const{
    flt tot = 0;
    for(unsigned int i=1; i<torsions.size(); i++){
        tot -= torsions[i] * i * pow(costheta, flt(i-1));
    }
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    return tot;
}

flt dihedral::energy(const Vec &r1, const Vec &r2, 
                   const Vec &r3) const{
    
    // The two normals to the planes
    Vec n1 = r1.cross(r2);
    Vec n2 = r2.cross(r3);
    
    flt n1mag = n1.mag();
    flt n2mag = n2.mag();
    
    if (n1mag == 0 or n2mag == 0) return 1; 
    // if one plane is ill-defined, then we have no torsion angle

    flt costheta = -n1.dot(n2) / n1mag / n2mag;
    
    flt tot = 0;
    for(unsigned int i=0; i<torsions.size(); i++){
        //~ std::cout << "curE: " << torsions[i] * pow(costheta, i) << " E:" << tot << '\n';
        tot += torsions[i] * pow(costheta, flt(i));
    }
    //~ std::cout << "costheta (E)" << costheta << " E:" << tot << '\n';
    
    return tot;
}

interMolPair::interMolPair(vector<atomgroup*> groupvec, interactpair* ipair)
    : groups(groupvec), pair(ipair) {};

flt interMolPair::energy(){
    flt E = 0;
    vector<atomgroup*>::iterator g1;
    vector<atomgroup*>::iterator g2;
    for(g1 = groups.begin(); g1 < groups.end(); g1++){
        g2 = g1;
        for(g2++; g2 < groups.end(); g2++)
        for(uint i1 = 0; i1 < (*g1)->N(); i1++)
        for(uint i2 = 0; i2 < (*g2)->N(); i2++){
            Vec r = diff((**g1)[i1].x, (**g2)[i2].x);
            E += pair->energy(r);
        }
    }
    return E;
};

void interMolPair::setForces(){
    vector<atomgroup*>::iterator g1;
    vector<atomgroup*>::iterator g2;
    for(g1 = groups.begin(); g1 < groups.end(); g1++){
        g2 = g1;
        for(g2++; g2 < groups.end(); g2++)
        for(uint i1 = 0; i1 < (*g1)->N(); i1++)
        for(uint i2 = 0; i2 < (*g2)->N(); i2++){
            atom & a1 = (**g1)[i1];
            atom & a2 = (**g2)[i2];
            Vec r = diff((**g1)[i1].x, (**g2)[i2].x);
            Vec force = pair->forces(r);
            a1.f += force;
            a2.f -= force;
        }
    }
};


intraMolNNPair::intraMolNNPair(vector<atomgroup*> groupvec, interactpair* ipair)
    : groups(groupvec), pair(ipair) {};

flt intraMolNNPair::energy(){
    flt E = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-1; i++){
            Vec r = diff(g[i].x, g[i+1].x);
            E += pair->energy(r);
        }
    }
    return E;
};

void intraMolNNPair::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-1; i++){
            atom & a1 = g[i];
            atom & a2 = g[i+1];
            Vec r = diff(g[i].x, g[i+1].x);
            Vec force = pair->forces(r);
            a1.f += force;
            a2.f -= force;
        }
    }
};


intraMolPairs::intraMolPairs(vector<atomgroup*> groupvec, 
            interactpair* ipair, cuint s)
    : groups(groupvec), pair(ipair), skip(s) {};

flt intraMolPairs::energy(){
    flt E = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-skip-1; i++)
        for(uint j = i+skip+1; j < g.N(); j++){
            Vec r = diff(g[i].x, g[j].x);
            E += pair->energy(r);
        }
    }
    return E;
};

void intraMolPairs::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-skip-1; i++)
        for(uint j = i+skip+1; j < g.N(); j++){
            atom & a1 = g[i];
            atom & a2 = g[j];
            Vec r = diff(g[i].x, g[j].x);
            Vec force = pair->forces(r);
            a1.f += force;
            a2.f -= force;
        }
    }
};


intraMolNNTriple::intraMolNNTriple(vector<atomgroup*> groupvec, interacttriple* itrip)
    : groups(groupvec), trip(itrip) {};

flt intraMolNNTriple::energy(){
    flt E = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-2; i++){
            E += trip->energy(diff(g[i+1].x, g[i].x), diff(g[i+1].x, g[i+2].x));
        }
    }
    return E;
};

void intraMolNNTriple::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-2; i++){
            atom & a0 = g[i];
            atom & a1 = g[i+1];
            atom & a2 = g[i+2];
            Vec r = diff(a0.x, a1.x);
            Nvector<Vec, 3> force = trip->forces(diff(a1.x, a0.x), diff(a1.x, a2.x));
            a0.f += force[0];
            a1.f += force[1];
            a2.f += force[2];
        }
    }
};

intraMolNNQuad::intraMolNNQuad(vector<atomgroup*> groupvec, interactquad* iquad)
    : groups(groupvec), quad(iquad) {};

flt intraMolNNQuad::energy(){
    flt E = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-3; i++){
            E += quad->energy(diff(g[i+1].x, g[i].x), diff(g[i+2].x, g[i+1].x),
                                                 diff(g[i+3].x, g[i+2].x));
        }
    }
    return E;
};

void intraMolNNQuad::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.N()-3; i++){
            //~ Nvector<Vec, 4> force = trip->forces(g.diff(i+1, i), 
                        //~ g.diff(i+2, i+1), g.diff(i+3, i+2));
            Nvector<Vec, 4> force = quad->forces(diff(g[i+1].x, g[i].x), 
                        diff(g[i+2].x, g[i+1].x), diff(g[i+3].x, g[i+2].x));
            for(uint j=0; j<4; j++) g[i+j].f += force[j];
        }
    }
};

flt singletpairs::energy(){
    flt E=0;
    vector<atompair>::iterator it;
    for(it = atoms.begin(); it < atoms.end(); it++){
        E += inter->energy(diff(it->first().x, it->last().x));
    }
    return E;
};

void singletpairs::setForces(){
    vector<atompair>::iterator it;
    for(it = atoms.begin(); it < atoms.end(); it++){
            atom & a1 = it->first();
            atom & a2 = it->last();
            Vec force = inter->forces(diff(a1.x, a2.x));
            a1.f += force;
            a2.f -= force;
    }
};

flt interactgroup::energy(){
    flt E=0;
    vector<interaction*>::iterator it;
    for(it = inters.begin(); it < inters.end(); it++){
        E += (*it)->energy();
    }
    return E;
};

void interactgroup::setForces(){
    vector<interaction*>::iterator it;
    for(it = inters.begin(); it < inters.end(); it++){
        (*it)->setForces();
    }
};

bondpairs::bondpairs(vector<bondgrouping> pairs) : pairs(pairs){};

flt bondpairs::energy(){
    flt E=0;
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = diff(it->a1->x, it->a2->x);
        E += spring(it->k, it->x0).energy(r);
    }
    return E;
}

void bondpairs::setForces(){
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        Vec r = diff(atom1.x, atom2.x);
        Vec f = spring(it->k, it->x0).forces(r);
        atom1.f += f;
        atom2.f -= f;
    }
}

angletriples::angletriples(vector<anglegrouping> triples) : triples(triples){};

flt angletriples::energy(){
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

void angletriples::setForces(){
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
    }
}


LJgroup::LJgroup(flt cutoff, vector<LJdata> atoms, set<atompair, atompaircomp> pairs) : 
        atoms(atoms), pairs(pairs){
    if(cutoff > 0) LJ = new LJcutoff(1,1,cutoff);
    else LJ = new LJforce(1,1);
};

flt LJgroup::energy(){
    flt E = 0;
    vector<LJdata>::iterator it1, it2;
    for(it1 = atoms.begin(); it1 < atoms.end(); it1++){
        for(it2 = it1; it2 < atoms.end(); it2++){
            atom *a1 = it1->atm, *a2 = it2->atm;
            if (a1 == a2) continue;
            if (has_pair(a1, a2)) continue;
            flt sigma = .5 * (it1->sigma + it2->sigma);
            flt epsilon = sqrt(it1->epsilon + it2->epsilon); //optimize
            E += epsilon * LJ->energy(diff(a1->x, a2->x).mag() / sigma); //optimize
        }
    }
    return E;
}

void LJgroup::setForces(){
    flt E = 0;
    vector<LJdata>::iterator it1, it2;
    for(it1 = atoms.begin(); it1 < atoms.end(); it1++){
        for(it2 = it1; it2 < atoms.end(); it2++){
            atom *a1 = it1->atm, *a2 = it2->atm;
            if (a1 == a2) continue;
            if (has_pair(a1, a2)) continue;
            flt sigma = .5 * (it1->sigma + it2->sigma);
            flt epsilon = sqrt(it1->epsilon + it2->epsilon); //optimize
            Vec r = diff(a1->x, a2->x);
            flt rmag = r.mag();
            flt fmag = LJ->forces(rmag / sigma) * epsilon / sigma;
            Vec f = r * (fmag / rmag);
            a1->f += f;
            a2->f -= f;
        }
    }
    return;
}