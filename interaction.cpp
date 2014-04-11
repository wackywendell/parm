#include "interaction.hpp"

flt spring::energy(const Vec r){
    flt m = r.mag();
    flt l = m - x0;
    return .5 * springk * l*l;
}

Vec spring::forces(const Vec r){
    flt m = r.mag();
    if(m < 1e-32) return Vec();
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
    if(costheta > 1){
        //~ cout << "resetting! " << costheta << endl;
        costheta = 1;
    }
    else if(costheta < -1){
        //~ cout << "resetting! " << costheta << endl;
        costheta = -1;
    }
    if(!usecos) return springk*powflt(acos(costheta) - theta0,2)/2;
    else return springk*powflt(costheta - cos(theta0),2)/2;
}

Nvector<Vec, 3> bondangle::forces(const Vec& r1, const Vec& r2){
    flt r1mag = r1.mag();
    flt r2mag = r2.mag();
    
    flt costheta = r1.dot(r2) / r1mag / r2mag;
    if(costheta > 1){
        //~ cout << "resetting! " << costheta << endl;
        costheta = 1;
    }
    else if(costheta < -1){
        //~ cout << "resetting! " << costheta << endl;
        costheta = -1;
    }
    flt theta = acos(costheta);
    //theta is now the angle between x1 and x2
    
    flt fmag;
    if(usecos) fmag = -springk*(cos(theta0) - costheta)*sin(theta);
    else fmag = springk*(theta0 - theta); // torque magnitude
    // We have V = \frac{1}{2}k(\theta-\theta_{0})^{2}
    // Then -f = grad V = \frac{k}{r}(\theta-\theta_{0})\hat{\theta}
    // first we get the direction:
    Nvector<Vec, 3> force;
    if (fmag == 0) {return force;}
    Vec f0 = r2.perpto(r1);
    Vec f2 = r1.perpto(r2);
    if ((f0.sq() <= 1e-30) or (f2.sq() <= 1e-30)) {return force;}
    force[0] = f0;
    force[0].normalize();
    force[2] = f2;
    force[2].normalize();
    
    // now we get magnitude: 
    force[0] *= fmag/r1mag;
    force[2] *= fmag/r2mag;
    
    //~ if(!(force[0].sq() < 1e6)){
        //~ cout << "theta: " << theta << "  theta0: " << theta0 << endl;
        //~ cout << "fmag: " << fmag << "  r1mag: " << r1mag << "  r2mag: " << r2mag << endl;
        //~ cout << "f0: " << f0 << "  force[0]: " << force[0] << endl;
        //~ cout << "f2: " << f2 << "  force[2]: " << force[2] << endl;
    //~ }
    
    //~ cout << force[2] << x2 << "force(2).x2: " << force[2].dot(x2) << endl;
    force[1] = -(force[0] + force[2]);
    // The direction of the force on the first atom (f0) is 
    // perpendicular to x1, and same for f2.
    // **FIXED** its possible that x1 = +/-x2, and then x1.perp(x2) = 0
    // and then we get a divide by zero error.
    
    //~ assert(force[0].sq() <= 1e7);
    //~ assert(force[1].sq() <= 1e7);
    //~ assert(force[2].sq() <= 1e7);
    
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
    flt sqq = sqrtflt(q);
    
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
     Rapaport says costheta = p/sqrtflt(q); we add a negative for the cosine.
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
        
    
    //~ flt mag = sqrtflt(derivs[0].sq() +derivs[1].sq() + derivs[2].sq() +
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
    unsigned int cosmx = (unsigned int)(coscoeffs.size());
    for(unsigned int i=1; i < cosmx; i++){
        tot += coscoeffs[i] * i * powflt(costheta, flt(i-1));
    }
    //~ cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    return tot;
}

flt dihedral::dudcostheta(const flt theta) const{
    flt tot = 0;
    unsigned int cosmx = (unsigned int)(coscoeffs.size());
    unsigned int sinmx = (unsigned int)(sincoeffs.size());
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;
    if(usepow) {
        flt costheta = cos(theta), sintheta = sin(theta);
        flt cottheta = -costheta / sintheta;
        for(unsigned int i=1; i < mx; i++){
            if (i < cosmx) tot += coscoeffs[i] * i * powflt(costheta, flt(i-1));
            if (i < sinmx) tot += sincoeffs[i] * i * cottheta
                                    * powflt(sintheta, flt(i-1));
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
    
    unsigned int cosmx = (unsigned int)(coscoeffs.size());
    unsigned int sinmx = (unsigned int)(sincoeffs.size());
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;
    
    flt tot = 0;
    for(unsigned int i=0; i < mx; i++){
        if(usepow) {
            if(i < cosmx) tot += coscoeffs[i] * powflt(costheta, flt(i));
            if(i < sinmx) tot += sincoeffs[i] * powflt(sintheta, flt(i));
        } else {
            if(i < cosmx) tot += coscoeffs[i] * cos(i * ang);
            if(i < sinmx) tot += sincoeffs[i] * sin(i * ang);
        }
    }
    
    return tot;
}
#endif
//~ flt interactgroup::energy(Box &box){
    //~ flt E=0;
    //~ vector<interaction*>::iterator it;
    //~ for(it = inters.begin(); it < inters.end(); it++){
        //~ E += (*it)->energy(box);
    //~ }
    //~ return E;
//~ };

//~ void interactgroup::setForces(Box &box){
    //~ vector<interaction*>::iterator it;
    //~ for(it = inters.begin(); it < inters.end(); it++){
        //~ (*it)->setForces(box);
    //~ }
//~ };

bondgrouping::bondgrouping(flt k, flt x0, atom* a1, atom* a2, 
        BondDiffType diff, OriginBox *box) :
            k(k), x0(x0), a1(a1), a2(a2), diff_type(diff){
    if(diff == FIXEDBOX){
        assert(box != NULL);
        fixed_box = box->box_round(a1->x, a2->x);
    }
};

Vec bondgrouping::diff(Box &box) const{
    switch(diff_type){
        case BOXED:
            return box.diff(a1->x, a2->x);
        case UNBOXED:
            return a1->x - a2->x;
        case FIXEDBOX:
            OriginBox &obox = (OriginBox &) box;
            return obox.diff(a1->x, a2->x, fixed_box);
    }
    return Vec()*NAN;
};

bondpairs::bondpairs(vector<bondgrouping> pairs, bool zeropressure) : 
        zeropressure(zeropressure), pairs(pairs){};

bondpairs::bondpairs(bool zeropressure) : 
        zeropressure(zeropressure){};

bool bondpairs::add_or_replace(bondgrouping b){
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        if(((b.a1 == it->a1) and (b.a2 == it->a2)) or
            ((b.a1 == it->a2) and (b.a2 == it->a1))){
                *it = b;
                return true;
            }
        }
    pairs.push_back(b);
    return false;
};

bool bondpairs::replace(flt k, flt x0, atom* a1, atom* a2){
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        if(((a1 == it->a1) and (a2 == it->a2)) or
            ((a1 == it->a2) and (a2 == it->a1))){
                it->k = k;
                it->x0 = x0;
                return true;
            }
        }
    return false;
};

flt bondpairs::energy(Box &box){
    flt E=0;
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = it->diff(box);
        E += spring(it->k, it->x0).energy(r);
    }
    return E;
}

void bondpairs::setForces(Box &box){
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        Vec r = it->diff(box);
        Vec f = spring(it->k, it->x0).forces(r);
        //~ assert(f.sq() < 10000000);
        atom1.f += f;
        atom2.f -= f;
    }
}

flt bondpairs::setForcesGetPressure(Box &box){
    if(zeropressure){
        setForces(box);
        return 0.0;
    }
    flt P=0;
    vector<bondgrouping>::iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        Vec r = it->diff(box);
        Vec f = spring(it->k, it->x0).forces(r);
        //~ assert(f.sq() < 10000000);
        atom1.f += f;
        atom2.f -= f;
        if(it->diff_type == UNBOXED) continue;
        P += f.dot(r);
    }
    return P;
}

flt bondpairs::pressure(Box &box){
    if(zeropressure) return 0;
    vector<bondgrouping>::iterator it;
    flt P=0;
    for(it = pairs.begin(); it < pairs.end(); it++){
        if(it->diff_type == UNBOXED) continue;
        Vec r = it->diff(box);
        Vec f = spring(it->k, it->x0).forces(r);
        P += f.dot(r);
    }
    return P;
}

flt bondpairs::mean_dists(Box &box) const{
    flt dist=0;
    uint N=0;
    vector<bondgrouping>::const_iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = it->diff(box);
        dist += abs(r.mag() - it->x0);
        N++;
    }
    return dist/N;
}

flt bondpairs::std_dists(Box &box) const{
    flt stds=0;
    uint N=0;
    vector<bondgrouping>::const_iterator it;
    for(it = pairs.begin(); it < pairs.end(); it++){
        Vec r = it->diff(box);
        flt curdist = r.mag() - it->x0;
        stds += curdist*curdist;
        N++;
    }
    return sqrtflt(stds/N);
}

angletriples::angletriples(vector<anglegrouping> triples) : triples(triples){};

flt angletriples::energy(Box &box){
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

void angletriples::setForces(Box &box){
    vector<anglegrouping>::iterator it;
    for(it = triples.begin(); it < triples.end(); it++){
        atom & atom1 = *it->a1;
        atom & atom2 = *it->a2;
        atom & atom3 = *it->a3;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom2.x, atom3.x);
        Nvector<Vec,3> f = bondangle(it->k, it->x0).forces(r1, r2);
        assert(f[0].sq() < 1e8);
        assert(f[1].sq() < 1e8);
        assert(f[2].sq() < 1e8);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        assert(f[0].sq() < 1e8);
        assert(f[1].sq() < 1e8);
        assert(f[2].sq() < 1e8);
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
    return sqrtflt(stds/N);
};

#ifdef VEC3D
dihedrals::dihedrals(vector<dihedralgrouping> ds) : groups(ds){};

flt dihedrals::energy(Box &box){
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

void dihedrals::setForces(Box &box){
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
    //~ return sqrtflt(stds/N);
//~ };
#endif

void pairlist::clear(){
    map<const atomid, set<atomid> >::iterator mapiter;
    for(mapiter = pairs.begin(); mapiter != pairs.end(); mapiter++){
        mapiter->second.clear();
    }
};

LJsimple::LJsimple(flt cutoff, vector<LJatom> atms) : atoms(atms){};

flt LJsimple::energy(Box &box){
    flt E = 0;
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec dist = box.diff(pair.atom1.x(), pair.atom2.x());
        E += LJrepulsive::energy(dist, pair.sigma, pair.epsilon);
    }
    return E;
};

void LJsimple::setForces(Box &box){
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec r = box.diff(pair.atom1.x(), pair.atom2.x());
        Vec f = LJrepulsive::forces(r, pair.sigma, pair.epsilon);
        pair.atom1.f() += f;
        pair.atom2.f() -= f;
    }
};

flt LJsimple::pressure(Box &box){
    flt P=0;
    vector<LJatom>::iterator it;
    vector<LJatom>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        LJpair pair = LJpair(*it, *it2);
        Vec r = box.diff(pair.atom1.x(), pair.atom2.x());
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

flt Charges::energy(Box &box){
    flt E = 0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec dist = box.diff(it->x(), it2->x());
        E += k*electricScreened::energy(dist.mag(), (it->q) * (it2->q), screen);
    }
    return E;
};

void Charges::setForces(Box &box){
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec r = box.diff(it->x(), it2->x());
        Vec f = electricScreened::forces(r, (it->q) * (it2->q), screen)*k;
        it->f() += f;
        it2->f() -= f;
    }
};

flt Charges::pressure(Box &box){
    flt P=0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for(it = atoms.begin(); it != atoms.end(); it++)
    for(it2 = atoms.begin(); it2 != it; it2++){
        if (ignorepairs.has_pair(*it, *it2)) continue;
        Vec r = box.diff(it->x(), it2->x());
        Vec f = electricScreened::forces(r, (it->q) * (it2->q), screen)*k;
        P += r.dot(f);
    }
    return P;
};

flt SoftWall::energy(Box &box){
    flt E=0;
    typename vector<WallAtom>::iterator it;
    for(it = group.begin(); it != group.end(); it++){
        atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm)*2;
        if(dist > it->sigma) continue;
        E += it->epsilon * powflt(1 - (dist/(it->sigma)), expt)/expt/2.0;
        // Note that normally you have ε(1-r/σ)^n for 2 particles.
        // We divide by 2 because now there is only one particle, pushing
        // on its mirror image; the force should be the same as if the 
        // mirror image was there, so the energy needs to be half
    }
    return E;
};

void SoftWall::setForces(Box &box){
    lastf = 0;
    typename vector<WallAtom>::iterator it;
    for(it = group.begin(); it != group.end(); it++){
        atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm)*2;
        if(dist > it->sigma) continue;
        flt f = it->epsilon * powflt(1 - (dist/(it->sigma)), expt - 1.0);
        lastf += -f;
        a.f += norm*f;
    }
};

flt SoftWall::setForcesGetPressure(Box &box){
    flt p=0;
    lastf = 0;
    typename vector<WallAtom>::iterator it;
    for(it = group.begin(); it != group.end(); it++){
        atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm)*2;
        if(dist > it->sigma) continue;
        flt f = it->epsilon * powflt(1 - (dist/(it->sigma)), expt - 1.0);
        p += dist*f/2;
        lastf += -f;
        a.f += norm*f;
    }
    return p;
};

flt SoftWall::pressure(Box &box){
    flt p = 0;
    typename vector<WallAtom>::iterator it;
    for(it = group.begin(); it != group.end(); it++){
        atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm)*2;
        if(dist > it->sigma) continue;
        flt f = it->epsilon * powflt(1 - (dist/(it->sigma)), expt - 1.0);
        p += dist*f/2;
    }
    return p;
};

SpheroCylinderDiff SCPair::NearestLoc(Box &box){
    // see Abreu, Charlles RA and Tavares, Frederico W. and Castier, Marcelo, "Influence of particle shape on the packing and on the segregation of spherocylinders via Monte Carlo simulations", Powder Technology 134, 1 (2003), pp. 167–180.
    // Uses that notation, just i -> 1, j -> 2, adds s1,s2
    SpheroCylinderDiff diff;
    
    atom &a1 = p1.first();
    atom &a1p = p1.last();
    atom &a2 = p2.first();
    atom &a2p = p2.last();
    Vec r1 = (a1.x + a1p.x)/2, r2 = (a2.x + a2p.x)/2;
    Vec s1 = (a1p.x - a1.x), s2 = (a2p.x - a2.x);
    //flt myl1 = s1.mag(), myl2 = s2.mag();
    Vec u1 = s1.norm(), u2=s2.norm();
    diff.r = box.diff(r2, r1);
    
    flt u1u2 = u1.dot(u2);
    //~ cout << "u1: " << u1 << "  u2: " << u2 << "  u1u2:" << u1u2 << "\n";
    
    flt u1u2sq = u1u2*u1u2;
    flt u1r12 = u1.dot(diff.r), u2r12 = u2.dot(diff.r);
    //~ cout << "r: " << diff.r << "  u1r12: " << u1r12 << "  u2r12: " << u2r12 << "\n";
    
    // Where the two lines would intersect
    flt lambda1p, lambda2p;
    
    if(abs(1-u1u2sq) < 1e-8){
        // They are too close to parallel, so we just say the "middle" 
        // of the two spherocylinders (r12*l1/(l1+l2)) projected onto their axes (u1)
        lambda1p = u1r12*l1/(l1+l2);
        // symmetry would be u2r21/2, but r21 = -r12
        lambda2p = -u2r12*l2/(l1+l2);
        //~ cout << "Too small, adjusted." << '\n';
    } else {
        // Where the two lines actually would intersect
        lambda1p = (u1r12 - (u1u2*u2r12))/(1-u1u2sq);
        lambda2p = ((u1u2*u1r12) - u2r12)/(1-u1u2sq);
    }

    //~ cout << "l1p: " << lambda1p << "  l2p: " << lambda2p << "\n";
    
    flt lambda1s=lambda1p, lambda2s=lambda2p;
    
    flt L1 = abs(lambda1p) - (l1/2);
    flt L2 = abs(lambda2p) - (l2/2);
    if(L1 > 0 or L2 > 0){
        if(L2 > L1){
            lambda2s = copysignflt(l2/2, lambda2p);
            lambda1s = u1r12 + (lambda2s*u1u2);
            //~ cout << "new l2s: " << lambda2s << "  l1s: " << lambda1s;
            lambda1s = confineRange(-l1/2, lambda1s, l1/2);
            //~ cout << " -> " << lambda1s << "\n";
        } else {
            lambda1s = copysignflt(l1/2, lambda1p);
            lambda2s = -u2r12 + (lambda1s*u1u2);
            //~ cout << "new l1s: " << lambda1s << "  l2s: " << lambda2s;
            lambda2s = confineRange(-l2/2, lambda2s, l2/2);
            //~ cout << " -> " << lambda2s << "\n";
        }
    }
    
    diff.lambda1 = lambda1s;
    diff.lambda2 = lambda2s;
    diff.delta = box.diff(r2 + (u2*lambda2s), r1 + (u1*lambda1s));
    
    return diff;
};

void SCPair::applyForce(Box &box, Vec f, SpheroCylinderDiff diff, flt IoverM1, flt IoverM2){
    atom &a1 = p1.first();
    atom &a1p = p1.last();
    atom &a2 = p2.first();
    atom &a2p = p2.last();
    Vec r1 = (a1.x + a1p.x)/2, r2 = (a2.x + a2p.x)/2;
    Vec s1 = (a1.x - a1p.x), s2 = (a2.x - a2p.x);
    flt M1 = a1.m + a1p.m;
    flt M2 = a2.m + a2p.m;
    
    a1.f -= f/2; // note that the force on a1 is half the total force, this carries through to atau1
    a1p.f -= f/2;
    a2.f += f/2;
    a2p.f += f/2;
    
    Vec t1 = s1*(diff.lambda1/l1);
    Vec atau1 = s1.cross(t1.cross(f)) / (-2*IoverM1*M1);
    //~ cout << "t1: " << t1 << "  atau1: " << atau1 << endl;
    // Formula says (t1×f)×s1 / 2I
    // -I because it should be (t1×f)×s1, but we wrote s1×(t1×f)
    // 4 and not 2 because f is the force on the whole thing, we only want half
    a1.f += atau1 * a1.m;
    a1p.f -= atau1 * a1p.m;
    
    Vec t2 = s2*(diff.lambda2/l2);
    Vec atau2 = s2.cross(t2.cross(f)) / (2*IoverM2*M2); // 2 because it should be -f
    a2.f += atau2 * a2.m;
    a2p.f -= atau2 * a2p.m;
    
    //~ cout << "t2: " << t2 << "  atau2: " << atau1 << endl;
};

flt SCSpringList::energy(Box &box){
    flt E = 0;
    for(uint i = 0; i < scs->pairs() - 1; i++){
        atompair pi = scs->pair(i);
        for(uint j = i+1; j < scs->pairs(); j++){
            atompair pj = scs->pair(j);
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, ls[i], ls[j]);
            SpheroCylinderDiff diff = scp.NearestLoc(box);
            //~ cout << "SCSpringList diff delta: " << diff.delta << '\n';
            //~ cout << "SCSpringList diff lambdas: " << diff.lambda1 << "    " << diff.lambda2 << '\n';
            E += scp.energy(box, diff);
            //~ cout << "SCSpringList energy: " << E << '\n';
        }
    }
    return E;
};

void SCSpringList::setForces(Box &box){
    for(uint i = 0; i < scs->pairs() - 1; i++){
        atompair pi = scs->pair(i);
        flt l1 = ls[i];
        for(uint j = i+1; j < scs->pairs(); j++){
            atompair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.NearestLoc(box);
            Vec f = scp.forces(box, diff);
            scp.applyForce(box, f, diff, l1*l1/4, l2*l2/4);
        }
    }
};
