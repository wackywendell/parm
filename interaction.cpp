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

Vec atomgroup::angmomentum(const Vec &loc) const{
    flt curmass;
    Vec tot = Vec();
    Vec newloc;
    for(uint i=0; i<size(); i++){
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
    for(uint i=0; i<size(); i++){
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

void atomgroup::addrot(Vec omega, Vec loc){
    for(uint i=0; i<size(); i++){
        Vec r = (*this)[i].x - loc;
        (*this)[i].v += r.cross(omega);
    }
};

void atomgroup::resetForces(){
    for(uint i=0; i<size(); i++){
        (*this)[i].f = Vec(0,0,0);
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

//~ void atomgroup::vverlet2(const flt dt){
    //~ for(uint i=0; i<size(); i++){
        //~ (*this)[i].v += (*this)[i].a * (dt/2);
    //~ }
//~ };

metagroup::metagroup(vector<atomgroup*> groups){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i < g.size(); i++){
            atoms.push_back(& (g[i]));
            masses.push_back(g.getmass(i));
        }
    }
};

atomid metagroup::get_id(atom* a){
    uint n=0;
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


electricScreened::electricScreened(const flt screenLength, const flt q1, 
            const flt q2, const flt cutoff) : screen(screenLength), 
            q1(q1), q2(q2), cutoff(cutoff){};

flt electricScreened::energy(const flt r, const flt qaqb, const flt screen, const flt cutoff){
    if(cutoff <= 0) return exp(-r/screen) * qaqb / r;
    if(r > cutoff) return 0;
    return exp(-r/screen) * qaqb / r - energy(r, qaqb, screen, 0);
};

Vec electricScreened::forces(const Vec &r, const flt qaqb, const flt screen, const flt cutoff){
    flt d = r.mag();
    if(cutoff > 0 and d > cutoff) return Vec(0,0,0);
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
    
    flt costheta = p/sqrt(q);
    // costheta=1 corresponds to atoms 1 and 4 on opposite sides of the bond (zigzag)
    // costheta=-1 corresponds to a C shape
    
    
    flt dcostheta = -dudcostheta(costheta); // F = -dU/d(costheta)
    //~ if(abs(costheta) < .7)
        //~ cout << "forces cos: " << costheta << " getcos: " << getcos(r1,r2,r3)
             //~ << " dcostheta: " << dcostheta << '\n';
    
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
        tot += torsions[i] * i * pow(costheta, flt(i-1));
    }
    //~ cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
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

    return -(n1.dot(n2) / n1mag / n2mag);
};

flt dihedral::energy(const Vec &r1, const Vec &r2, 
                   const Vec &r3) const{
    
    // The two normals to the planes
    //~ Vec n1 = r1.cross(r2);
    //~ Vec n2 = r2.cross(r3);
    //~ 
    //~ flt n1mag = n1.mag();
    //~ flt n2mag = n2.mag();
    //~ 
    //~ if (n1mag == 0 or n2mag == 0) return -100; 
    // if one plane is ill-defined, then we have no torsion angle

    flt costheta = getcos(r1,r2,r3);
    
    flt tot = 0;
    for(unsigned int i=0; i<torsions.size(); i++){
        tot += torsions[i] * pow(costheta, flt(i));
        //~ if(abs(costheta) < .75)
        //~ std::cout << i << " t: " << torsions[i] << " curE: "
        //~ << torsions[i] * pow(costheta, i) << " E:" << tot << '\n';
    }
    //~ if(abs(costheta) < .75)
    //~ std::cout << "costheta: " << costheta << " E:" << tot << '\n';
    
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
        for(uint i1 = 0; i1 < (*g1)->size(); i1++)
        for(uint i2 = 0; i2 < (*g2)->size(); i2++){
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
        for(uint i1 = 0; i1 < (*g1)->size(); i1++)
        for(uint i2 = 0; i2 < (*g2)->size(); i2++){
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
        for(uint i = 0; i < g.size()-1; i++){
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
        for(uint i = 0; i < g.size()-1; i++){
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
        for(uint i = 0; i < g.size()-skip-1; i++)
        for(uint j = i+skip+1; j < g.size(); j++){
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
        for(uint i = 0; i < g.size()-skip-1; i++)
        for(uint j = i+skip+1; j < g.size(); j++){
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
        for(uint i = 0; i < g.size()-2; i++){
            E += trip->energy(diff(g[i+1].x, g[i].x), diff(g[i+1].x, g[i+2].x));
        }
    }
    return E;
};

void intraMolNNTriple::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git < groups.end(); git++){
        atomgroup & g = **git;
        for(uint i = 0; i < g.size()-2; i++){
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
        for(uint i = 0; i < g.size()-3; i++){
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
        for(uint i = 0; i < g.size()-3; i++){
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
};

flt angletriples::mean_dists() const{
    flt dist=0;
    uint N=0;
    vector<anglegrouping>::const_iterator it;
    for(it = triples.begin(); it < triples.end(); it++){
        Vec r1 = diff(it->a1->x, it->a2->x);
        Vec r2 = diff(it->a3->x, it->a2->x);
        flt theta = acos(r1.dot(r2) / r1.mag() / r2.mag());
        dist += abs(theta - it->x0);
        N++;
    }
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

dihedrals::dihedrals(vector<dihedralgrouping> ds) : groups(ds){};

flt dihedrals::energy(){
    flt E=0;
    vector<dihedralgrouping>::iterator it;
    for(it = groups.begin(); it < groups.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        atom & atom4 = *it->a4;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom3.x, atom2.x);
        Vec r3 = diff(atom4.x, atom3.x);
        E += dihedral(it->nums).energy(r1,r2,r3);
    }
    return E;
}

void dihedrals::setForces(){
    vector<dihedralgrouping>::iterator it;
    for(it = groups.begin(); it < groups.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        atom & atom4 = *it->a4;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom3.x, atom2.x);
        Vec r3 = diff(atom4.x, atom3.x);
        Nvector<Vec,4> f = dihedral(it->nums).forces(r1, r2, r3);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        atom4.f += f[3];
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
        Vec r1 = diff(atom1.x, atom2.x);
        Vec r2 = diff(atom2.x, atom3.x);
        Vec r3 = diff(atom3.x, atom4.x);
        flt cosine = dihedral(it->nums).getcos(r1, r2, r3);
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

void pairlist::clear(){
    map<const atomid, set<atomid> >::iterator mapiter;
    for(mapiter = pairs.begin(); mapiter != pairs.end(); mapiter++){
        mapiter->second.clear();
    }
};

neighborlist::neighborlist(atomgroup &group, const flt innerradius,
            const flt outerradius, pairlist ignore) :
                skinradius(outerradius), critdist(innerradius),
                ignorepairs(ignore), ignorechanged(false){
    atoms.resize(group.size());
    lastlocs.resize(group.size());
    for(uint i=0; i<group.size(); i++){
        atomid a = group.get_id(i);
        ignorepairs.ensure(a);
        atoms[a.n()] = a;
        lastlocs[a.n()] = a.x();
    }
    update_list(true);
    updatenum = 1;
};

atomid neighborlist::get_id(atom* a){
    for(vector<atomid>::iterator it=atoms.begin(); it!=atoms.end(); it++)
        if ((*it) == a) return *it;
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
    
    vector<atomid>::iterator atomit, atomj;
    
    if(not force and not ignorechanged){ // check if we need to update
        for(atomit=atoms.begin(); atomit!=atoms.end(); atomit++){
            curdist = (atomit->x() - lastlocs[atomit->n()]).mag();
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
    for(atomit=atoms.begin(); atomit!=atoms.end(); atomit++){
        lastlocs[atomit->n()] = atomit->x();
        for(atomj=atoms.begin(); atomj!=atomit; atomj++){
            if (ignorepairs.has_pair(*atomit, *atomj)) continue;
            if(diff(atomit->x(), atomj->x()).mag() < skinradius)
                curpairs.push_back(idpair(*atomit, *atomj));
        }
    }
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

LJsimple::LJsimple(flt cutoff, vector<LJatom> atms) : atoms(atms){};

flt LJsimple::energy(){
    flt E = 0;
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec dist = diff(pair.atom1.x(), pair.atom2.x());
        E += LJrepulsive::energy(dist, pair.sigma, pair.epsilon);
    }
    return E;
};

void LJsimple::setForces(){
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec r = diff(pair.atom1.x(), pair.atom2.x());
        Vec f = LJrepulsive::forces(r, pair.sigma, pair.epsilon);
        pair.atom1.f() += f;
        pair.atom2.f() -= f;
    }
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

flt Charges::energy(){
    flt E = 0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec dist = diff(it->x(), it2->x());
        E += k*electricScreened::energy(dist.mag(), (it->q) * (it2->q), screen);
    }
    return E;
};

void Charges::setForces(){
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec r = diff(it->x(), it2->x());
        Vec f = electricScreened::forces(r, (it->q) * (it2->q), screen)*k;
        it->f() += f;
        it2->f() -= f;
    }
};
