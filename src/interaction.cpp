#include "interaction.hpp"

flt Spring::energy(const Vec r) {
    flt m = r.norm();
    flt l = m - x0;
    return .5 * springk * l * l;
}

Vec Spring::forces(const Vec r) {
    flt m = r.norm();
    if (m < 1e-32) return Vec::Zero();
    flt fmag = (x0 - m) * springk;
    return r * (fmag / m);
}

ElectricScreened::ElectricScreened(const flt screenLength, const flt q1,
                                   const flt q2, const flt cutoff)
    : screen(screenLength),
      q1(q1),
      q2(q2),
      cutoff(cutoff),
      cutoffE(exp(-cutoff / screenLength) * q1 * q2 / cutoff){};

flt ElectricScreened::energy(const flt r, const flt qaqb, const flt screen,
                             const flt cutoff) {
    if (cutoff <= 0) return exp(-r / screen) * qaqb / r;
    if (r > cutoff) return 0;
    return exp(-r / screen) * qaqb / r - energy(cutoff, qaqb, screen, 0);
};

Vec ElectricScreened::forces(const Vec r, const flt qaqb, const flt screen,
                             const flt cutoff) {
    flt d = r.norm();
    if (cutoff > 0 and d > cutoff) return Vec::Zero();
    flt fmag = exp(-d / screen) * (qaqb / d) * (1 / d + 1 / screen);
    return r * (fmag / d);
};

FixedForceRegionAtom::FixedForceRegionAtom(AtomID a, Vec dir, vector<flt> bound,
                                           vector<flt> Fs)
    : AtomID(a), direction(dir.normalized()), boundaries(bound), Fs(Fs) {
    if (boundaries.size() != Fs.size() - 1) {
        throw std::invalid_argument(
            "FixedForceRegionAtom: boundaries.size() != Fs.size() - 1");
    } else if (boundaries.size() <= 0) {
        throw std::invalid_argument(
            "FixedForceRegionAtom: boundaries.size() == 0, use FixedForceAtom");
    }
};

flt FixedForceRegionAtom::energy(Box &box) {
    flt dist = (*this)->x.dot(direction);
    if (dist < boundaries[0]) return (boundaries[0] - dist) * Fs[0];
    flt E = 0;
    for (uint i = 0; i < boundaries.size(); ++i) {
        if ((i + 1 >= boundaries.size()) or (dist < boundaries[i + 1])) {
            E -= Fs[i + 1] * (dist - boundaries[i]);
            break;
        }
        E -= Fs[i + 1] * (boundaries[i + 1] - boundaries[i]);
    }

    return E;
};

void FixedForceRegionAtom::set_force(Box &box) {
    flt dist = (*this)->x.dot(direction);
    flt F = Fs[Fs.size() - 1];
    for (uint i = 0; i < boundaries.size(); ++i) {
        if (dist < boundaries[i]) {
            F = Fs[i];
            break;
        }
    }
    (*this)->f += direction * F;
};

flt BondAngle::energy(const Vec &r1, const Vec &r2) {
    flt costheta = r1.dot(r2) / r1.norm() / r2.norm();
    if (costheta > 1)
        costheta = 1;
    else if (costheta < -1)
        costheta = -1;
    if (!usecos)
        return springk * pow(acos(costheta) - theta0, 2) / 2;
    else
        return springk * pow(costheta - cos(theta0), 2) / 2;
}

barray<Vec, 3> BondAngle::forces(const Vec &r1, const Vec &r2) {
    flt r1mag = r1.norm();
    flt r2mag = r2.norm();
    flt costheta = r1.dot(r2) / r1mag / r2mag;
    if (costheta > 1)
        costheta = 1;
    else if (costheta < -1)
        costheta = -1;
    flt theta = acos(costheta);
    // theta is now the angle between x1 and x2

    flt fmag;
    if (usecos)
        fmag = -springk * (cos(theta0) - costheta) * sin(theta);
    else
        fmag = springk * (theta0 - theta);  // torque magnitude
    // We have V = \frac{1}{2}k(\theta-\theta_{0})^{2}
    // Then -f = grad V = \frac{k}{r}(\theta-\theta_{0})\hat{\theta}
    // first we get the direction:
    barray<Vec, 3> force;
    if (fmag == 0) {
        return force;
    }
    Vec f0 = perpto(r2, r1);
    Vec f2 = perpto(r1, r2);
    if ((f0.squaredNorm() <= 1e-30) or (f2.squaredNorm() <= 1e-30)) {
        return force;
    }
    force[0] = f0;
    force[0].normalize();
    force[2] = f2;
    force[2].normalize();

    // now we get magnitude:
    force[0] *= fmag / r1mag;
    force[2] *= fmag / r2mag;

    //~ if(!(force[0].squaredNorm() < 1e6)){
    //~ cout << "theta: " << theta << "  theta0: " << theta0 << endl;
    //~ cout << "fmag: " << fmag << "  r1mag: " << r1mag << "  r2mag: " << r2mag
    //<< endl;
    //~ cout << "f0: " << f0 << "  force[0]: " << force[0] << endl;
    //~ cout << "f2: " << f2 << "  force[2]: " << force[2] << endl;
    //~ }

    //~ cout << force[2] << x2 << "force(2).x2: " << force[2].dot(x2) << endl;
    force[1] = -(force[0] + force[2]);
    // The direction of the force on the first Atom (f0) is
    // perpendicular to x1, and same for f2.
    // **FIXED** its possible that x1 = +/-x2, and then x1.perp(x2) = 0
    // and then we get a divide by zero error.

    //~ assert(force[0].squaredNorm() <= 1e7);
    //~ assert(force[1].squaredNorm() <= 1e7);
    //~ assert(force[2].squaredNorm() <= 1e7);

    return force;
}

#ifdef VEC3D
Dihedral::Dihedral(const vector<flt> cvals, const vector<flt> svals,
                   bool usepow)
    : cos_coefficients(cvals), sincoeffs(svals), usepow(usepow) {}

DihedralDerivs Dihedral::dr_dcostheta(const Vec &r1, const Vec &r2,
                                      const Vec &r3) {
    // Taken from Rapaport "Art of Molecular Dynamics Simulation" p.279
    // The expressions and notation are very close to that of the book.

    // Note that Rappaport defines it as such:
    /*
     The Dihedral angle is defined as the angle between the
planes formed by atoms 1,2,3 and 2,3,4 measured in the plane normal to the 2–3
bond; it is zero when all four atoms are coplanar and atoms 1 and 4 are on
opposite
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
    flt qa = c[0][0] * c[1][1] - c[0][1] * c[0][1];
    flt qb = c[1][1] * c[2][2] - c[1][2] * c[1][2];
    flt q = qa * qb;
    flt sqq = sqrt(q);

    flt t1 = p;
    flt t2 = c[0][0] * c[1][2] - c[0][1] * c[0][2];
    flt t3 = c[0][1] * c[0][1] - c[0][0] * c[1][1];
    flt t4 = c[1][1] * c[2][2] - c[1][2] * c[1][2];
    flt t5 = c[0][2] * c[1][2] - c[0][1] * c[2][2];
    flt t6 = -p;

    DihedralDerivs dd;

    /*
     Rapaport: The Dihedral angle is defined as the angle between the
     planes formed by atoms 1,2,3 and 2,3,4 measured in the plane
     normal to the 2–3 bond; it is zero when all four atoms are
     coplanar and atoms 1 and 4 are on opposite sides of the bond.

     Note that this is the *opposite* of the chemical definition.

     flt const0 = c[1][1]/(sqq * qa);
     flt const3 = c[1][1]/(sqq * qb);

     We add a negative in at the beginning of those two to give us the
     chemical definition.
     */

    flt const0 = -c[1][1] / (sqq * qa);
    flt const3 = -c[1][1] / (sqq * qb);

    dd.derivs[0] = (r1 * t1 + r2 * t2 + r3 * t3) * const0;
    dd.derivs[3] = (r1 * t4 + r2 * t5 + r3 * t6) * const3;

    dd.derivs[1] = dd.derivs[0] * (-1 - c[0][1] / c[1][1]) +
                   dd.derivs[3] * (c[1][2] / c[1][1]);
    dd.derivs[2] = dd.derivs[0] * (c[0][1] / c[1][1]) -
                   dd.derivs[3] * (1 + c[1][2] / c[1][1]);

    /*
    Rapaport says costheta = p/sqrt(q); we add a negative for the cosine.
    */
    dd.costheta = -p / sqq;
    // costheta =-1 corresponds to atoms 1 and 4 on opposite sides of the bond
    // (zigzag)
    // costheta = 1 corresponds to a C shape

    return dd;
}

barray<Vec, 4> Dihedral::forces(const Vec &r1, const Vec &r2,
                               const Vec &r3) const {
    DihedralDerivs dd = dr_dcostheta(r1, r2, r3);
    flt dcostheta;
    if (sincoeffs.empty() and !usepow) {
        dcostheta = dU_dcostheta_cos(dd.costheta);  // F = -dU/d(costheta)
    } else {
        dcostheta = dU_dcostheta(get_angle(r1, r2, r3));
    }

    for (uint i = 0; i < 4; i++)
        dd.derivs[i] *= -dcostheta;  // F = -dU/d(costheta)
                                     //~ assert(derivs[0].squaredNorm() < 1e8);
                                     //~ assert(derivs[1].squaredNorm() < 1e8);
                                     //~ assert(derivs[2].squaredNorm() < 1e8);
                                     //~ assert(derivs[3].squaredNorm() < 1e8);

    //~ flt mag = sqrt(derivs[0].squaredNorm() +derivs[1].squaredNorm() +
    // derivs[2].squaredNorm() +
    //~ derivs[3].squaredNorm());
    //~
    //~ std::cout << "costheta:" << costheta << " dcos:" << dcostheta
    //~ << " derivs:" << derivs  << " : " << mag << std::endl;
    return dd.derivs;

    // pea79, dun92
    // Pear, M. R. and Weiner, J. H., Brownian dynamics study of a polymer chain
    // of linked rigid bodies, J. Chem. Phys. 71 (1979) 212.
    // Dunn, J. H., Lambrakos, S. G., Moore, P. G., and Nagumo, M., An algorithm
    // for calculating intramolecular angle-dependent forces on vector
    // computers, J. Comp. Phys. 100 (1992) 17.
}

flt Dihedral::dU_dcostheta_cos(const flt costheta) const {
    assert(sincoeffs.empty());
    assert(!usepow);
    flt tot = 0;
    unsigned int cosmx = (unsigned int)(cos_coefficients.size());
    for (unsigned int i = 1; i < cosmx; ++i) {
        tot += cos_coefficients[i] * i * pow(costheta, flt(i - 1));
    }
    //~ cout << "dudcos tot: " << tot << ", cos: " << costheta << '\n';
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta <<
    //'\n';
    return tot;
}

flt Dihedral::dU_dcostheta(const flt theta) const {
    flt tot = 0;
    unsigned int cosmx = (unsigned int)(cos_coefficients.size());
    unsigned int sinmx = (unsigned int)(sincoeffs.size());
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;
    if (usepow) {
        flt costheta = cos(theta), sintheta = sin(theta);
        flt cottheta = -costheta / sintheta;
        for (unsigned int i = 1; i < mx; ++i) {
            if (i < cosmx)
                tot += cos_coefficients[i] * i * pow(costheta, flt(i - 1));
            if (i < sinmx)
                tot += sincoeffs[i] * i * cottheta * pow(sintheta, flt(i - 1));
        }
    } else {
        flt csctheta = 1 / sin(theta);
        for (unsigned int i = 1; i < mx; ++i) {
            flt cositheta = cos(i * theta), sinitheta = sin(i * theta);
            if (i < cosmx)
                tot += cos_coefficients[i] * csctheta * i * sinitheta;
            if (i < sinmx) tot -= sincoeffs[i] * csctheta * i * cositheta;
        }
    }
    //~ if(tot > 100) cout << "dudcos tot: " << tot << ", cos: " << costheta <<
    //'\n';
    return tot;
}

flt Dihedral::get_cos(const Vec &r1, const Vec &r2, const Vec &r3) {
    // The two normals to the planes
    Vec n1 = r1.cross(r2);
    Vec n2 = r2.cross(r3);
    //~ cout << r1 << ',' << r2  << ',' << r3 << "\n";
    flt n1mag = n1.norm();
    flt n2mag = n2.norm();

    if (n1mag == 0 or n2mag == 0) return -100;
    // if one plane is ill-defined, then we have no torsion angle

    return (n1.dot(n2) / n1mag / n2mag);
};

flt Dihedral::get_angle(const Vec &r1, const Vec &r2, const Vec &r3) {
    return atan2(r1.dot(r2.cross(r3)) * r2.norm(),
                 (r1.cross(r2).dot(r2.cross(r3))));
};

flt Dihedral::energy(const flt ang) const {
    flt costheta = (usepow ? cos(ang) : NAN);
    flt sintheta = (usepow ? sin(ang) : NAN);

    unsigned int cosmx = (unsigned int)(cos_coefficients.size());
    unsigned int sinmx = (unsigned int)(sincoeffs.size());
    unsigned int mx = cosmx > sinmx ? cosmx : sinmx;

    flt tot = 0;
    for (unsigned int i = 0; i < mx; ++i) {
        if (usepow) {
            if (i < cosmx) tot += cos_coefficients[i] * pow(costheta, flt(i));
            if (i < sinmx) tot += sincoeffs[i] * pow(sintheta, flt(i));
        } else {
            if (i < cosmx) tot += cos_coefficients[i] * cos(i * ang);
            if (i < sinmx) tot += sincoeffs[i] * sin(i * ang);
        }
    }

    return tot;
}
#endif

bool RandomForce::add(RandomForceAtom a, bool replace) {
    vector<RandomForceAtom>::iterator it;
    for (it = group.begin(); it < group.end(); ++it) {
        if ((*it) == (AtomRef)a) {
            if (replace) {
                *it = a;
                return true;
            } else {
                throw std::invalid_argument("Atoms already inserted.");
            }
        }
    }
    group.push_back(a);
    return false;
};

void RandomForce::set_forces(Box &box) {
    vector<RandomForceAtom>::iterator it;
    for (it = group.begin(); it < group.end(); ++it) {
        if (rand01() < 1.0 / it->freq) {
            RandomForceAtom &a = *it;
            Vec v = rand_vec();
            switch (a.force_type) {
                case FIXED:
                    a->f += v * (a.force_mag / v.norm());
                    continue;
                case UNIFORM:
                    a->f += v * (rand01() * a.force_mag / v.norm());
                    continue;
                case GAUSSIAN:
                    a->f += v * a.force_mag;
                    continue;
            }
        }
    };
};

BondGrouping::BondGrouping(flt k, flt x0, AtomID a1, AtomID a2,
                           BondDiffType diff, OriginBox *box)
    : k(k), x0(x0), a1(a1), a2(a2), diff_type(diff) {
    if (diff == FIXEDBOX) {
        assert(box != NULL);
        fixed_box = box->box_round(a1->x, a2->x);
    }
};

Vec BondGrouping::diff(Box &box) const {
    switch (diff_type) {
        case BOXED:
            return box.diff(a1->x, a2->x);
        case UNBOXED:
            return a1->x - a2->x;
        case FIXEDBOX:
            OriginBox &obox = (OriginBox &)box;
            return obox.diff(a1->x, a2->x, fixed_box);
    }
    return Vec::Zero() * NAN;
};

BondPairs::BondPairs(vector<BondGrouping> pairs, bool zeropressure)
    : zeropressure(zeropressure), pairs(pairs){};

BondPairs::BondPairs(bool zeropressure) : zeropressure(zeropressure){};

bool BondPairs::add(BondGrouping b, bool replace) {
    vector<BondGrouping>::iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        if (b.same_atoms(*it)) {
            if (replace) {
                *it = b;
                return true;
            } else {
                throw std::invalid_argument("Atoms already inserted.");
            }
        }
    }
    pairs.push_back(b);
    return false;
};

flt BondPairs::energy(Box &box) {
    flt E = 0;
    vector<BondGrouping>::iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        Vec r = it->diff(box);
        E += Spring(it->k, it->x0).energy(r);
    }
    return E;
}

void BondPairs::set_forces(Box &box) {
    vector<BondGrouping>::iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Vec r = it->diff(box);
        Vec f = Spring(it->k, it->x0).forces(r);
        //~ assert(f.squaredNorm() < 10000000);
        atom1.f += f;
        atom2.f -= f;
    }
}

flt BondPairs::set_forces_get_pressure(Box &box) {
    if (zeropressure) {
        set_forces(box);
        return 0.0;
    }
    flt P = 0;
    vector<BondGrouping>::iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Vec r = it->diff(box);
        Vec f = Spring(it->k, it->x0).forces(r);
        //~ assert(f.squaredNorm() < 10000000);
        atom1.f += f;
        atom2.f -= f;
        if (it->diff_type == UNBOXED) continue;
        P += f.dot(r);
    }
    return P;
}

flt BondPairs::pressure(Box &box) {
    if (zeropressure) return 0;
    vector<BondGrouping>::iterator it;
    flt P = 0;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        if (it->diff_type == UNBOXED) continue;
        Vec r = it->diff(box);
        Vec f = Spring(it->k, it->x0).forces(r);
        P += f.dot(r);
    }
    return P;
}

flt BondPairs::mean_dists(Box &box) const {
    flt dist = 0;
    uint N = 0;
    vector<BondGrouping>::const_iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        Vec r = it->diff(box);
        dist += abs(r.norm() - it->x0);
        N++;
    }
    return dist / N;
}

flt BondPairs::std_dists(Box &box) const {
    flt stds = 0;
    uint N = 0;
    vector<BondGrouping>::const_iterator it;
    for (it = pairs.begin(); it < pairs.end(); ++it) {
        Vec r = it->diff(box);
        flt curdist = r.norm() - it->x0;
        stds += curdist * curdist;
        N++;
    }
    return sqrt(stds / N);
}

AngleTriples::AngleTriples(vector<AngleGrouping> triples) : triples(triples){};

bool AngleTriples::add(AngleGrouping a, bool replace) {
    vector<AngleGrouping>::iterator it;
    for (it = triples.begin(); it < triples.end(); ++it) {
        if (a.same_atoms(*it)) {
            if (replace) {
                *it = a;
                return true;
            } else {
                throw std::invalid_argument("Atoms already inserted.");
            }
        }
    }
    triples.push_back(a);
    return false;
};

bool AngleTriples::add(flt k, AtomID a1, AtomID a2, AtomID a3, bool replace) {
    Vec r1 = diff(a2->x, a1->x);
    Vec r2 = diff(a2->x, a3->x);
    flt x0 = BondAngle::get_angle(r1, r2);
    return add(AngleGrouping(k, x0, a1, a2, a3), replace);
};

flt AngleTriples::energy(Box &box) {
    flt E = 0;
    vector<AngleGrouping>::iterator it;
    for (it = triples.begin(); it < triples.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Atom &atom3 = *it->a3;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom2.x, atom3.x);
        E += BondAngle(it->k, it->x0).energy(r1, r2);
    }
    return E;
}

void AngleTriples::set_forces(Box &box) {
    vector<AngleGrouping>::iterator it;
    for (it = triples.begin(); it < triples.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Atom &atom3 = *it->a3;
        Vec r1 = diff(atom2.x, atom1.x);
        Vec r2 = diff(atom2.x, atom3.x);
        barray<Vec, 3> f = BondAngle(it->k, it->x0).forces(r1, r2);
        assert(f[0].squaredNorm() < 1e8);
        assert(f[1].squaredNorm() < 1e8);
        assert(f[2].squaredNorm() < 1e8);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        assert(f[0].squaredNorm() < 1e8);
        assert(f[1].squaredNorm() < 1e8);
        assert(f[2].squaredNorm() < 1e8);
    }
};

flt AngleTriples::mean_dists() const {
    flt dist = 0;
    uint N = 0;
    vector<AngleGrouping>::const_iterator it;
    //~ cout << "angle diffs: ";
    cout.precision(3);
    for (it = triples.begin(); it < triples.end(); ++it) {
        Vec r1 = diff(it->a1->x, it->a2->x);
        Vec r2 = diff(it->a3->x, it->a2->x);
        flt theta = acos(r1.dot(r2) / r1.norm() / r2.norm());
        //~ flt curdist = abs(theta - it->x0);
        dist += abs(theta - it->x0);
        //~ if(curdist > .2){
        //~ flt E = BondAngle(it->k, it->x0).energy(r1,r2);
        //~ Vec f0 = BondAngle(it->k, it->x0).forces(r1,r2)[0];
        //~ cout << "(" << theta << "," << it->x0 << "," << abs(theta - it->x0)
        //~ << ";" << it->k << "," << E
        //~ << ",f:" << f0.norm() <<  ")" << ", ";
        //~ }
        N++;
    }
    //~ cout << "Total: " << N << '\n';
    return dist / N;
};

flt AngleTriples::std_dists() const {
    flt stds = 0;
    uint N = 0;
    vector<AngleGrouping>::const_iterator it;
    for (it = triples.begin(); it < triples.end(); ++it) {
        Vec r1 = diff(it->a1->x, it->a2->x);
        Vec r2 = diff(it->a3->x, it->a2->x);
        flt theta = acos(r1.dot(r2) / r1.norm() / r2.norm());
        flt curdist = theta - (it->x0);
        //~ cout << "angles: " << theta << " x0: " << it->x0
        //~ << " cur: " << curdist << "\n";
        stds += curdist * curdist;
        N++;
    }
    return sqrt(stds / N);
};

#ifdef VEC3D
Dihedrals::Dihedrals(vector<DihedralGrouping> ds) : groups(ds){};

flt Dihedrals::energy(Box &box) {
    flt E = 0;
    vector<DihedralGrouping>::iterator it;
    for (it = groups.begin(); it < groups.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Atom &atom3 = *it->a3;
        Atom &atom4 = *it->a4;
        Vec r1 = DihedralGrouping::diff(atom2.x, atom1.x);
        Vec r2 = DihedralGrouping::diff(atom3.x, atom2.x);
        Vec r3 = DihedralGrouping::diff(atom4.x, atom3.x);
        E += it->dih.energy(r1, r2, r3);
    }
    return E;
}

void Dihedrals::set_forces(Box &box) {
    vector<DihedralGrouping>::iterator it;
    for (it = groups.begin(); it < groups.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Atom &atom3 = *it->a3;
        Atom &atom4 = *it->a4;
        Vec r1 = DihedralGrouping::diff(atom2.x, atom1.x);
        Vec r2 = DihedralGrouping::diff(atom3.x, atom2.x);
        Vec r3 = DihedralGrouping::diff(atom4.x, atom3.x);
        barray<Vec, 4> f = it->dih.forces(r1, r2, r3);
        atom1.f += f[0];
        atom2.f += f[1];
        atom3.f += f[2];
        atom4.f += f[3];
        //~ flt maxf = 1000000;
        //~ if(f[0].squaredNorm() > maxf or f[1].squaredNorm() > maxf or
        // f[2].squaredNorm() > maxf
        //~ or f[3].squaredNorm() > maxf){
        //~ cout << "Dihedral overload: " << r1 << r2 << r3 << " :: " <<
        //~ f[0] << f[1] << f[2] << f[3] << "\n";
        //~ cout << "Dihedral overload energy: " <<
        // Dihedral(it->nums).energy(r1,r2,r3) << "\n";
        //~ }
    }
};

flt Dihedrals::mean_dists() const {
    flt dist = 0;
    uint N = 0;
    vector<DihedralGrouping>::const_iterator it;
    for (it = groups.begin(); it < groups.end(); ++it) {
        Atom &atom1 = *it->a1;
        Atom &atom2 = *it->a2;
        Atom &atom3 = *it->a3;
        Atom &atom4 = *it->a4;
        Vec r1 = DihedralGrouping::diff(atom1.x, atom2.x);
        Vec r2 = DihedralGrouping::diff(atom2.x, atom3.x);
        Vec r3 = DihedralGrouping::diff(atom3.x, atom4.x);
        flt cosine = Dihedral::get_cos(r1, r2, r3);
        dist += cosine;
        N++;
    }
    return dist / N;
};
//~
//~ flt Dihedrals::std_dists() const{
//~ flt stds=0;
//~ uint N=0;
//~ vector<AngleGrouping>::const_iterator it;
//~ for(it = triples.begin(); it < triples.end(); ++it){
//~ Vec r1 = diff(it->a1->x, it->a2->x);
//~ Vec r2 = diff(it->a3->x, it->a2->x);
//~ flt theta = acos(r1.dot(r2) / r1.norm() / r2.norm());
//~ flt curdist = theta - (it->x0);
// cout << "angles: " << theta << " x0: " << it->x0
// << " cur: " << curdist << "\n";
//~ stds += curdist*curdist;
//~ N++;
//~ }
//~ return sqrt(stds/N);
//~ };
#endif

Charges::Charges(flt screen, flt k, vector<Charged> atms)
    : atoms(atms), screen(screen), k(k){};

AtomID Charges::get_id(Atom *a) {
    for (vector<Charged>::iterator it = atoms.begin(); it != atoms.end(); ++it)
        if ((*it) == a) return *it;
    return AtomID();
};

flt Charges::energy(Box &box) {
    flt E = 0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for (it = atoms.begin(); it != atoms.end(); ++it)
        for (it2 = atoms.begin(); it2 != it; ++it2) {
            if (ignorepairs.has_pair(*it, *it2)) continue;
            Vec dist = box.diff((*it)->x, (*it2)->x);
            E += k * ElectricScreened::energy(dist.norm(), (it->q) * (it2->q),
                                              screen);
        }
    return E;
};

void Charges::set_forces(Box &box) {
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for (it = atoms.begin(); it != atoms.end(); ++it)
        for (it2 = atoms.begin(); it2 != it; ++it2) {
            if (ignorepairs.has_pair(*it, *it2)) continue;
            Vec r = box.diff((*it)->x, (*it2)->x);
            Vec f = ElectricScreened::forces(r, (it->q) * (it2->q), screen) * k;
            (*it)->f += f;
            (*it2)->f -= f;
        }
};

flt Charges::pressure(Box &box) {
    flt P = 0;
    vector<Charged>::iterator it;
    vector<Charged>::iterator it2;
    for (it = atoms.begin(); it != atoms.end(); ++it)
        for (it2 = atoms.begin(); it2 != it; ++it2) {
            if (ignorepairs.has_pair(*it, *it2)) continue;
            Vec r = box.diff((*it)->x, (*it2)->x);
            Vec f = ElectricScreened::forces(r, (it->q) * (it2->q), screen) * k;
            P += r.dot(f);
        }
    return P;
};

flt SoftWall::energy(Box &box) {
    flt E = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm) * 2;
        if (dist > it->sigma) continue;
        E += it->epsilon * pow(1 - (dist / (it->sigma)), expt) / expt / 2.0;
        // Note that normally you have ε(1-r/σ)^n for 2 particles.
        // We divide by 2 because now there is only one particle, pushing
        // on its mirror image; the force should be the same as if the
        // mirror image was there, so the energy needs to be half
    }
    return E;
};

void SoftWall::set_forces(Box &box) {
    lastf = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm) * 2;
        if (dist > it->sigma) continue;
        flt f = it->epsilon * pow(1 - (dist / (it->sigma)), expt - 1.0);
        lastf += -f;
        a.f += norm * f;
    }
};

flt SoftWall::set_forces_get_pressure(Box &box) {
    flt p = 0;
    lastf = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm) * 2;
        if (dist > it->sigma) continue;
        flt f = it->epsilon * pow(1 - (dist / (it->sigma)), expt - 1.0);
        p += dist * f / 2;
        lastf += -f;
        a.f += norm * f;
    }
    return p;
};

flt SoftWall::pressure(Box &box) {
    flt p = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        flt dist = r.dot(norm) * 2;
        if (dist > it->sigma) continue;
        flt f = it->epsilon * pow(1 - (dist / (it->sigma)), expt - 1.0);
        p += dist * f / 2;
    }
    return p;
};
flt SoftWallCylinder::energy(Box &box) {
    flt E = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        r -= axis * (r.dot(axis));
        flt dist = (radius - r.norm()) * 2;
        if (dist > it->sigma) continue;
        E += it->epsilon * pow(1 - (dist / (it->sigma)), expt) / expt / 2.0;
        // Note that normally you have ε(1-r/σ)^n for 2 particles.
        // We divide by 2 because now there is only one particle, pushing
        // on its mirror image; the force should be the same as if the
        // mirror image was there, so the energy needs to be half
    }
    return E;
};

void SoftWallCylinder::set_forces(Box &box) {
    lastf = 0;
    vector<WallAtom>::iterator it;
    for (it = group.begin(); it != group.end(); ++it) {
        Atom &a = **it;
        Vec r = box.diff(a.x, loc);
        r -= axis * (r.dot(axis));
        flt rmag = r.norm();
        flt dist = (radius - rmag) * 2;
        if (dist > it->sigma) continue;
        flt f = it->epsilon * pow(1 - (dist / (it->sigma)), expt - 1.0);
        a.f -= r * (f / rmag);  // equal to a.f += (-r.normalized()) * f;
        lastf += f;
    }
};

flt SoftWallCylinder::set_forces_get_pressure(Box &box) {
    throw std::runtime_error(
        "SoftWallCylinder::set_forces_get_pressure not implemented");
    return 0;
};

flt SoftWallCylinder::pressure(Box &box) {
    throw std::runtime_error("SoftWallCylinder::pressure not implemented");
    return 0;
};

flt SCAtomVec::volume(flt diameter, flt length, uint dim) {
    flt cap_V = M_PI * pow(diameter, (flt)dim) / (2 * dim);
    flt shaft_V;
    if (dim == 2) {
        shaft_V = diameter * (length - diameter);
    } else if (dim == 3) {
        shaft_V = diameter * diameter * M_PI / 4 * length;
    } else {
        throw std::invalid_argument("SCAtomVec::volume: needs dim=2 or dim=3");
        return 0;
    }
    return cap_V + shaft_V;
};

SpheroCylinderDiff SCPair::nearest_location(Box &box) {
    // see Abreu, Charlles RA and Tavares, Frederico W. and Castier, Marcelo,
    // "Influence of particle shape on the packing and on the segregation of
    // spherocylinders via Monte Carlo simulations", Powder Technology 134, 1
    // (2003), pp. 167–180.
    // Uses that notation, just i -> 1, j -> 2, adds s1,s2
    SpheroCylinderDiff diff;

    Atom &a1 = *p1.first();
    Atom &a1p = *p1.last();
    Atom &a2 = *p2.first();
    Atom &a2p = *p2.last();
    Vec r1 = (a1.x + a1p.x) / 2, r2 = (a2.x + a2p.x) / 2;
    Vec s1 = (a1p.x - a1.x), s2 = (a2p.x - a2.x);
    // flt myl1 = s1.norm(), myl2 = s2.norm();
    Vec u1 = s1.normalized(), u2 = s2.normalized();
    diff.r = box.diff(r2, r1);

    flt u1u2 = u1.dot(u2);
    //~ cout << "u1: " << u1 << "  u2: " << u2 << "  u1u2:" << u1u2 << "\n";

    flt u1u2sq = u1u2 * u1u2;
    flt u1r12 = u1.dot(diff.r), u2r12 = u2.dot(diff.r);
    //~ cout << "r: " << diff.r << "  u1r12: " << u1r12 << "  u2r12: " << u2r12
    //<< "\n";

    // Where the two lines would intersect
    flt lambda1p, lambda2p;

    if (abs(1 - u1u2sq) < 1e-8) {
        // They are too close to parallel, so we just say the "middle"
        // of the two spherocylinders (r12*l1/(l1+l2)) projected onto their axes
        // (u1)
        lambda1p = u1r12 * l1 / (l1 + l2);
        // symmetry would be u2r21/2, but r21 = -r12
        lambda2p = -u2r12 * l2 / (l1 + l2);
        //~ cout << "Too small, adjusted." << '\n';
    } else {
        // Where the two lines actually would intersect
        lambda1p = (u1r12 - (u1u2 * u2r12)) / (1 - u1u2sq);
        lambda2p = ((u1u2 * u1r12) - u2r12) / (1 - u1u2sq);
    }

    //~ cout << "l1p: " << lambda1p << "  l2p: " << lambda2p << "\n";

    flt lambda1s = lambda1p, lambda2s = lambda2p;

    flt L1 = abs(lambda1p) - (l1 / 2);
    flt L2 = abs(lambda2p) - (l2 / 2);
    if (L1 > 0 or L2 > 0) {
        if (L2 > L1) {
            lambda2s = copysign(l2 / 2, lambda2p);
            lambda1s = u1r12 + (lambda2s * u1u2);
            //~ cout << "new l2s: " << lambda2s << "  l1s: " << lambda1s;
            lambda1s = confine_range(-l1 / 2, lambda1s, l1 / 2);
            //~ cout << " -> " << lambda1s << "\n";
        } else {
            lambda1s = copysign(l1 / 2, lambda1p);
            lambda2s = -u2r12 + (lambda1s * u1u2);
            //~ cout << "new l1s: " << lambda1s << "  l2s: " << lambda2s;
            lambda2s = confine_range(-l2 / 2, lambda2s, l2 / 2);
            //~ cout << " -> " << lambda2s << "\n";
        }
    }

    diff.lambda1 = lambda1s;
    diff.lambda2 = lambda2s;
    diff.delta = box.diff(r2 + (u2 * lambda2s), r1 + (u1 * lambda1s));

    return diff;
};

void SCPair::apply_force(Box &box, Vec f, SpheroCylinderDiff diff, flt IoverM1,
                         flt IoverM2) {
    Atom &a1 = *p1.first();
    Atom &a1p = *p1.last();
    Atom &a2 = *p2.first();
    Atom &a2p = *p2.last();
    Vec s1 = (a1.x - a1p.x), s2 = (a2.x - a2p.x);
    flt M1 = a1.m + a1p.m;
    flt M2 = a2.m + a2p.m;

    a1.f -= f / 2;  // note that the force on a1 is half the total force, this
                    // carries through to atau1
    a1p.f -= f / 2;
    a2.f += f / 2;
    a2p.f += f / 2;

    Vec t1 = s1 * (diff.lambda1 / l1);
    Vec atau1 = cross(s1, cross(t1, f)) / (-2 * IoverM1 * M1);
    //~ cout << "t1: " << t1 << "  atau1: " << atau1 << endl;
    // Formula says (t1×f)×s1 / 2I
    // -I because it should be (t1×f)×s1, but we wrote s1×(t1×f)
    // 4 and not 2 because f is the force on the whole thing, we only want half
    a1.f += atau1 * a1.m;
    a1p.f -= atau1 * a1p.m;

    Vec t2 = s2 * (diff.lambda2 / l2);
    Vec atau2 = cross(s2, cross(t2, f)) /
                (2 * IoverM2 * M2);  // 2 because it should be -f
    a2.f += atau2 * a2.m;
    a2p.f -= atau2 * a2p.m;

    //~ cout << "t2: " << t2 << "  atau2: " << atau1 << endl;
};

flt SCSpringList::energy(Box &box) {
    flt E = 0;
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, ls[i], ls[j]);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            //~ cout << "SCSpringList diff delta: " << diff.delta << '\n';
            //~ cout << "SCSpringList diff lambdas: " << diff.lambda1 << "    "
            //<< diff.lambda2 << '\n';
            E += scp.energy(box, diff);
            //~ cout << "SCSpringList energy: " << E << '\n';
        }
    }
    return E;
};

void SCSpringList::set_forces(Box &box) {
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        flt l1 = ls[i];
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            Vec f = scp.forces(box, diff);
            scp.apply_force(box, f, diff, l1 * l1 / 4, l2 * l2 / 4);
        }
    }
};

flt SCSpringList::set_forces_get_pressure(Box &box) {
    flt P = 0;
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        flt l1 = ls[i];
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            Vec f = scp.forces(box, diff);
            scp.apply_force(box, f, diff, l1 * l1 / 4, l2 * l2 / 4);
            flt rdotf = diff.r.dot(f);
            P += rdotf;
        }
    }
    return P;
};

flt SCSpringList::pressure(Box &box) {
    flt P = 0;
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        flt l1 = ls[i];
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            Vec f = scp.forces(box, diff);
            flt rdotf = diff.r.dot(f);
            P += rdotf;
        }
    }
    return P;
};

Matrix SCSpringList::set_forces_get_stress(Box &box) {
    Matrix stress = Matrix::Zero();
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        flt l1 = ls[i];
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            Vec f = scp.forces(box, diff);
            scp.apply_force(box, f, diff, l1 * l1 / 4, l2 * l2 / 4);
            stress += diff.r * f.transpose();
        }
    }
    return stress;
};

Matrix SCSpringList::stress(Box &box) {
    Matrix stress = Matrix::Zero();
    barray<uint, 2> pair;
    for (uint i = 0; i < scs->pairs() - 1; ++i) {
        IDPair pi = scs->pair(i);
        flt l1 = ls[i];
        pair[0] = i;
        for (uint j = i + 1; j < scs->pairs(); ++j) {
            pair[1] = j;
            if (ignore_list.count(pair) > 0) continue;
            IDPair pj = scs->pair(j);
            flt l2 = ls[j];
            SCSpringPair scp = SCSpringPair(pi, pj, eps, sig, l1, l2);
            SpheroCylinderDiff diff = scp.nearest_location(box);
            Vec f = scp.forces(box, diff);
            stress += diff.r * f.transpose();
        }
    }
    return stress;
};

flt SCSpringList::volume() {
    flt V = 0;
    for (uint i = 0; i < scs->pairs(); ++i) {
        V += SCAtomVec::volume(sig, ls[i]);
    };
    return V;
};
