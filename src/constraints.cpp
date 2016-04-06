#include "constraints.hpp"

void CoordCOMConstraint::apply_positions(Box& box) {
    Vec com = a->com() - loc;

    for (uint i = 0; i < a->size(); i++) {
        Atom& atm = (*a)[i];
        for (uint j = 0; j < NDIM; j++) {
            if (not fixed[j]) continue;
            atm.x[j] -= com[j];
        }
    }
}

void CoordCOMConstraint::apply_velocities(Box& box) {
    Vec com_velocity = a->com_velocity();

    for (uint i = 0; i < a->size(); i++) {
        Atom& atm = (*a)[i];
        for (uint j = 0; j < NDIM; j++) {
            if (not fixed[j]) continue;
            atm.v[j] -= com_velocity[j];
        }
    }
}

void CoordCOMConstraint::apply_forces(Box& box) {
    Vec totf = Vec::Zero();
    for (uint i = 0; i < a->size(); i++) {
        totf += (*a)[i].f;
    }
    Vec tota = totf / a->mass();

    for (uint i = 0; i < a->size(); i++) {
        Atom& atm = (*a)[i];
        Vec df = (tota * (atm.m));
        for (uint j = 0; j < NDIM; j++) {
            if (not fixed[j]) continue;
            atm.f[j] -= df[j];
        }
    }
}

void RelativeConstraint::apply_positions(Box& box) {
    Vec dx = a2->x - a1->x;
    for (uint i = 0; i < NDIM; i++) {
        if (not fixed[i]) continue;
        a1->x[i] += dx[i] / 2;
        a2->x[i] -= dx[i] / 2;
        assert(abs(a2->x[i] - a1->x[i]) < 1e-5);
    }
}

void RelativeConstraint::apply_velocities(Box& box) {
    Vec dv = a2->v - a1->v;
    for (uint i = 0; i < NDIM; i++) {
        if (not fixed[i]) continue;
        a1->v[i] += dv[i] / 2;
        a2->v[i] -= dv[i] / 2;
        assert(abs(a2->v[i] - a1->v[i]) < 1e-5);
    }
}

void RelativeConstraint::apply_forces(Box& box) {
    flt mratio1 = a1->m / (a1->m + a2->m);
    flt mratio2 = a2->m / (a1->m + a2->m);
    Vec totf = a2->f + a1->f;
    for (uint i = 0; i < NDIM; i++) {
        if (not fixed[i]) continue;
        a1->f[i] = totf[i] * mratio1;
        a2->f[i] = totf[i] * mratio2;
    }
}

void DistConstraint::apply_positions(Box& box) {
    flt M = (a1->m + a2->m);
    flt mratio1 = a1->m / M;
    flt mratio2 = a2->m / M;

    Vec dx = a2->x - a1->x;
    flt dxmag = dx.norm();

    a1->x += dx * ((1 - dist / dxmag) * mratio2);
    a2->x -= dx * ((1 - dist / dxmag) * mratio1);
}

void DistConstraint::apply_velocities(Box& box) {
    Vec dx = a2->x - a1->x;
    flt dxmag = dx.norm();
    Vec dxnorm = dx / dxmag;

    Vec baddv = dxnorm * ((a2->v - a1->v).dot(dxnorm) / 2);
    a1->v += baddv;
    a2->v -= baddv;
    if ((a2->v - a1->v).dot(dxnorm) > 1e-8) {
        throw std::overflow_error("Velocities are not minimal.");
    }
}

void DistConstraint::apply_forces(Box& box) {
    Vec dx = a2->x - a1->x;
    flt dxmag = dx.norm();
    Vec dxnorm = dx / dxmag;

// TODO: Fix mass ratio stuff
#ifdef VEC2D
    flt omega = (a1->v.dot(dxnorm) - a2->v.dot(dxnorm)) / (2 * dxmag);
    Vec omega_cross_omega_cross_r = -omega * omega * dx;
#else
    // TODO: test this
    Vec omega = (a1->v.cross(dx) + a2->v.cross(dx)) / (2 * dxmag * dxmag);
    Vec omega_cross_omega_cross_r = omega.cross(omega.cross(dx));
#endif

    Vec baddf =
        dxnorm * ((a2->f - a1->f).dot(dxnorm) / 2) + omega_cross_omega_cross_r;
    a1->f += baddf;
    a2->f -= baddf;
    // assert((a2->f - a1->f).dot(dxnorm) < 1e-8);
    // if((a2->f - a1->f).dot(dxnorm) > 1e-8){
    //     throw std::overflow_error("Forces are not minimal.");
    // }
}

void LinearConstraint::set_lvec_com() {
    com = atms->com();
    uint sz = atms->size();
    lvec = Vec::Zero();

    for (uint i = 0; i < sz; i++) {
        flt chaindist = i * dist - lincom;
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;
        lvec += dx.normalized() * chaindist;
    }

    lvec.normalize();
}

void LinearConstraint::apply_positions(Box& box) {
    set_lvec_com();
    uint sz = atms->size();

    for (uint i = 0; i < sz; i++) {
        flt chaindist = i * dist - lincom;
        Vec dx = lvec * chaindist;
        Atom& ai = (*atms)[i];
        ai.x = com + dx;
    }
}

void LinearConstraint::apply_velocities(Box& box) {
    Vec com_velocity = atms->com_velocity();

    uint sz = atms->size();
#ifdef VEC3D
    Vec L = Vec::Zero();
    Vec omega = Vec::Zero();
#else
    flt L = 0;
    flt omega = 0;
#endif

    for (uint i = 0; i < sz; i++) {
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;
        L += cross(dx, ai.v) * ai.m;
    }

    lvec.normalize();
    omega = L / I;

    for (uint i = 0; i < sz; i++) {
        flt chaindist = i * dist - lincom;
        Vec dx = lvec * chaindist;
        Atom& ai = (*atms)[i];
        ai.v = com_velocity + cross(dx, omega);
    }
}

void LinearConstraint::apply_forces(Box& box) {
    Vec com_force = Vec::Zero();

    uint sz = atms->size();
#ifdef VEC3D
    Vec tau = Vec::Zero();
    Vec alpha = Vec::Zero();
#else
    flt tau = 0;
    flt alpha = 0;
#endif

    for (uint i = 0; i < sz; i++) {
        flt chaindist = i * dist - lincom;
        Vec dx = lvec * chaindist;
        Atom& ai = (*atms)[i];
        com_force += ai.f;
        tau += cross(dx, ai.f);
    }
    alpha = tau / I;

    for (uint i = 0; i < sz; i++) {
        flt chaindist = i * dist - lincom;
        Vec dx = lvec * chaindist;
        Atom& ai = (*atms)[i];
        ai.f = com_force + (cross(dx, alpha) * ai.m);
    }
}

#ifdef VEC3D
RigidConstraint::RigidConstraint(sptr<Box> box, sptr<AtomGroup> atms)
    : atms(atms),
      M(atms->mass()),
      MoI(atms->moment()),
      MoI_solver(MoI, Eigen::ComputeFullU | Eigen::ComputeFullV),
      expected(atms->size(), NDIM),
      rot(Matrix::Identity()),
      com(atms->com()),
      omega(Vec::Zero()),
      alpha(Vec::Zero()) {
    finite_or_throw(MoI);
    for (uint i = 0; i < atms->size(); i++) {
        expected.row(i) = ((*atms)[i].x - com);
    };
};

Matrix RigidConstraint::get_rotation() {
    com = atms->com();
    uint sz = atms->size();
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs(sz, NDIM);
    for (uint i = 0; i < sz; i++) {
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;
        locs.row(i) = dx;
    }
    finite_or_throw(expected);
    finite_or_throw(locs);
    return best_rotation_matrix(expected, locs);
}

void RigidConstraint::apply_positions(Box& box) {
    com = atms->com();

    uint sz = atms->size();
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs(sz, NDIM);

    for (uint i = 0; i < sz; i++) {
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;
        locs.row(i) = dx;
    }

    rot = best_rotation_matrix(expected, locs);

    for (uint i = 0; i < sz; i++) {
        Vec loc(expected.row(i));
        Vec dx(rot * loc);
        Atom& ai = (*atms)[i];
        ai.x = com + dx;
    }
}

void RigidConstraint::apply_velocities(Box& box) {
    Vec com_velocity = atms->com_velocity();
    uint sz = atms->size();
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs(sz, NDIM);

    Vec L = Vec::Zero();

    for (uint i = 0; i < sz; i++) {
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;

        // L = sum((r - R) × m v)
        L += cross(dx, ai.v) * ai.m;  // L, L/T, M -> L²M/T
    }

    omega = Vec(rot * MoI_solver.solve(rot.adjoint() * L));
    // Vec omega((I_inv * L).transpose());             // 1/ML², L²M/T -> 1/T

    for (uint i = 0; i < sz; i++) {
        Vec loc(expected.row(i));
        Vec dx(rot * loc);
        Atom& ai = (*atms)[i];
        finite_or_throw(ai.x);
        Vec omega_cross_r = cross(omega, dx);
        ai.v = com_velocity + omega_cross_r;  // L/T + L/T; V + ω × r
    }
}

void RigidConstraint::apply_forces(Box& box) {
    Vec com_force = Vec::Zero();

    uint sz = atms->size();
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs(sz, NDIM);

    Vec tau = Vec::Zero();

    for (uint i = 0; i < sz; i++) {
        Atom& ai = (*atms)[i];
        Vec dx = ai.x - com;
        com_force += ai.f;

        // L = sum((r - R) × m v)
        tau += cross(dx, ai.f);  // L, ML/T² -> ML²/T²
    }

    alpha = Vec(rot * MoI_solver.solve(rot.adjoint() * tau));
    // Vec alpha((I_inv * tau).transpose());           // 1/ML², ML²/T² -> 1/T²
    Vec comf_M = com_force / M;  // L/T²

    for (uint i = 0; i < sz; i++) {
        Vec loc(expected.row(i));
        Vec dx(rot * loc);
        Atom& ai = (*atms)[i];
        ai.x = com + dx;
        Vec omega_cross_r = cross(omega, dx);
        ai.f = (comf_M + cross(alpha, dx) + cross(omega, omega_cross_r)) *
               ai.m;  // (L/T² + L/T²) M; V + α × r
    }
}
#endif

ContactTracker::ContactTracker(sptr<Box> box, sptr<AtomGroup> atoms,
                               vector<flt> dists)
    : atoms(atoms),
      dists(dists),
      contacts(),
      breaks(0),
      formations(0),
      incontact(0) {
    //~ cout << "Making contact tracker." << endl;
    uint N = atoms->size();
    dists.resize(N);
    contacts.resize(N);
    for (uint i = 0; i < N; ++i) {
        contacts[i].resize(i, false);
    }
    update(*box);
    breaks = 0;
    formations = 0;
};

void ContactTracker::update(Box& box) {
    //~ cout << "Contact tracker update." << endl;
    uint N = atoms->size();
    incontact = 0;
    for (uint i = 0; i < N; ++i) {
        Vec ri = atoms->get(i).x;
        for (uint j = 0; j < i; ++j) {
            Vec rj = atoms->get(j).x;
            Vec dr = box.diff(ri, rj);
            bool curcontact = (dr.norm() <= ((dists[i] + dists[j]) / 2));
            if (curcontact && (!contacts[i][j]))
                formations++;
            else if (!curcontact && contacts[i][j])
                breaks++;
            if (curcontact) incontact++;

            contacts[i][j] = curcontact;
        }
    }
    //~ cout << "Contact tracker update done." << endl;
};

void EnergyTracker::update(Box& box) {
    if (n_skipped + 1 < n_skip) {
        n_skipped += 1;
        return;
    }

    n_skipped = 0;
    uint Natoms = atoms->size();
    flt curU = 0, curK = 0;
    for (uint i = 0; i < Natoms; ++i) {
        Atom& curatom = atoms->get(i);
        curK += curatom.v.squaredNorm() * curatom.m / 2;
    }

    vector<sptr<Interaction> >::iterator it;
    for (it = interactions.begin(); it != interactions.end(); ++it) {
        curU += (*it)->energy(box);
    }

    curU -= U0;
    Ks += curK;
    Us += curU;
    Es += curK + curU;
    Ksq += curK * curK;
    Usq += curU * curU;
    Esq += (curK + curU) * (curK + curU);
    N++;
};

void EnergyTracker::set_U0(Box& box) {
    flt curU = 0;
    vector<sptr<Interaction> >::iterator it;
    for (it = interactions.begin(); it != interactions.end(); ++it) {
        curU += (*it)->energy(box);
    }
    set_U0(curU);
};

RsqTracker1::RsqTracker1(AtomGroup& atoms, unsigned long skip, Vec com)
    : pastlocs(atoms.size(), NDIM),
      xyz2sums(atoms.size(), NDIM),
      xyz4sums(atoms.size(), NDIM),
      r4sums(atoms.size(), 0),
      skip(skip),
      count(0) {
    xyz2sums.setZero();
    xyz4sums.setZero();
    for (uint i = 0; i < atoms.size(); ++i) {
        pastlocs.row(i) = atoms[i].x - com;
    };
};

void RsqTracker1::reset(AtomGroup& atoms, Vec com) {
    pastlocs.resize(atoms.size(), NDIM);
    xyz2sums.setZero();
    xyz4sums.setZero();
    r4sums.assign(atoms.size(), 0);
    for (uint i = 0; i < atoms.size(); ++i) {
        pastlocs.row(i) = atoms[i].x - com;
    };
    count = 0;
};

bool RsqTracker1::update(Box& box, AtomGroup& atoms, unsigned long t, Vec com) {
    if (t % skip != 0) return false;

    for (uint i = 0; i < atoms.size(); ++i) {
        // flt dist = box.diff(atoms[i].x, pastlocs[i]).squaredNorm();
        Vec r = atoms[i].x - com;
        // We don't want the boxed distance - we want the actual distance moved!
        Vec pastr = pastlocs.row(i);
        Vec distance_squared = r - pastr;
        Vec distsqsq = Vec::Zero();
        flt dist4 = 0;
        for (uint j = 0; j < NDIM; ++j) {
            distance_squared[j] *= distance_squared[j];
            distsqsq[j] = distance_squared[j] * distance_squared[j];
            dist4 += distance_squared[j];
        };

        dist4 *= dist4;

        xyz2sums.row(i) += distance_squared;
        xyz4sums.row(i) += distsqsq;
        r4sums[i] += dist4;
        pastlocs.row(i) = r;
    };
    count += 1;
    return true;
};

Eigen::Matrix<flt, Eigen::Dynamic, NDIM> RsqTracker1::xyz2() {
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> means(xyz2sums.rows(), NDIM);
    means.setZero();
    for (uint i = 0; i < xyz2sums.rows(); ++i) {
        means.row(i) = xyz2sums.row(i) / ((flt)count);
    }
    return means;
};

Eigen::Matrix<flt, Eigen::Dynamic, NDIM> RsqTracker1::xyz4() {
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> means(xyz4sums.rows(), NDIM);
    means.setZero();
    for (uint i = 0; i < xyz4sums.rows(); ++i) {
        means.row(i) = xyz4sums.row(i) / ((flt)count);
    }
    return means;
};

vector<flt> RsqTracker1::r4() {
    vector<flt> means(r4sums.size(), 0);
    for (uint i = 0; i < r4sums.size(); ++i) {
        means[i] = r4sums[i] / ((flt)count);
    }
    return means;
};

RsqTracker::RsqTracker(sptr<AtomGroup> atoms, vector<unsigned long> ns,
                       bool usecom)
    : atoms(atoms), curt(0), usecom(usecom) {
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<unsigned long>::iterator n = ns.begin(); n != ns.end(); ++n) {
        singles.push_back(RsqTracker1(*atoms, *n, com));
    }
};

void RsqTracker::update(Box& box) {
    curt++;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        it->update(box, *atoms, curt, com);
    }
};

void RsqTracker::reset() {
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        it->reset(*atoms, com);
    }
};

vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > RsqTracker::xyz2() {
    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > vals;
    vals.reserve(singles.size());
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back(it->xyz2());
    }
    return vals;
};

vector<vector<flt> > RsqTracker::r2() {
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vector<flt> val;
        Eigen::Matrix<flt, Eigen::Dynamic, NDIM> xyz2 = it->xyz2();
        val.reserve(xyz2.size());

        for (uint i = 0; i < xyz2.rows(); i++) {
            val.push_back(xyz2.row(i).sum());
        }
        vals.push_back(val);
    }
    return vals;
};

vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > RsqTracker::xyz4() {
    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > vals;
    vals.reserve(singles.size());
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back(it->xyz4());
    }
    return vals;
};

vector<vector<flt> > RsqTracker::r4() {
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back(it->r4());
    }
    return vals;
};
vector<flt> RsqTracker::counts() {
    vector<flt> vals;
    vals.reserve(singles.size());
    for (vector<RsqTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back((flt)it->get_count());
    }
    return vals;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
ISFTracker1::ISFTracker1(AtomGroup& atoms, unsigned long skip, vector<flt> ks,
                         Vec com)
    : pastlocs(atoms.size(), NDIM),
      ISFsums(ks.size(), vector<barray<cmplx, NDIM> >()),
      ks(ks),
      skip(skip),
      count(0) {
    for (uint i = 0; i < atoms.size(); ++i) {
        pastlocs.row(i) = atoms[i].x - com;
    };

    for (uint ki = 0; ki < ks.size(); ++ki) {
        ISFsums[ki] =
            vector<barray<cmplx, NDIM> >(atoms.size(), barray<cmplx, NDIM>());
    };
};

void ISFTracker1::reset(AtomGroup& atoms, Vec com) {
    pastlocs.resize(atoms.size(), Eigen::NoChange);
    for (uint ki = 0; ki < ks.size(); ki++) {
        ISFsums[ki].assign(atoms.size(), barray<cmplx, NDIM>());
    }
    for (uint i = 0; i < atoms.size(); ++i) {
        pastlocs.row(i) = atoms[i].x - com;
    };
    count = 0;
};

bool ISFTracker1::update(Box& box, AtomGroup& atoms, unsigned long t, Vec com) {
    if (t % skip != 0) return false;

    for (uint i = 0; i < atoms.size(); ++i) {
        // flt dist = box.diff(atoms[i].x, pastlocs[i]).squaredNorm();
        Vec r = atoms[i].x - com;
        // We don't want the boxed distance - we want the actual distance moved!
        Vec pastr = pastlocs.row(i);
        Vec dr = r - pastr;
        for (uint ki = 0; ki < ks.size(); ++ki) {
            for (uint j = 0; j < NDIM; ++j) {
                //~ cout << "i: " << i << "  ki: " << ki << "  j: " << j;
                //~ cout << "  dr[j]: " << dr[j];
                //~ cout << "  ks[ki]: " << ks[ki];
                //~ cout << endl;

                ISFsums[ki][i][j] += exp(cmplx(0, ks[ki] * dr[j]));
            }
        }

        pastlocs.row(i) = r;
    };
    count += 1;
    return true;
};

vector<vector<cmplx> > ISFTracker1::ISFs() {
    vector<vector<cmplx> > means(ks.size(), vector<cmplx>());
    for (uint ki = 0; ki < ks.size(); ++ki) {
        means[ki].assign(ISFsums[ki].size(), 0);
        for (uint i = 0; i < ISFsums[ki].size(); ++i) {
            for (uint j = 0; j < NDIM; ++j) {
                means[ki][i] += ISFsums[ki][i][j];
            }
            means[ki][i] /= NDIM;
            means[ki][i] /= (flt)count;
        }
    }

    return means;
};

vector<vector<barray<cmplx, NDIM> > > ISFTracker1::ISFxyz() {
    vector<vector<barray<cmplx, NDIM> > > means(ks.size(),
                                                vector<barray<cmplx, NDIM> >());
    for (uint ki = 0; ki < ks.size(); ++ki) {
        means[ki].assign(ISFsums[ki].size(), barray<cmplx, NDIM>());
        for (uint i = 0; i < ISFsums[ki].size(); ++i) {
            for (uint j = 0; j < NDIM; ++j)
                means[ki][i][j] = ISFsums[ki][i][j] / (cmplx((flt)count, 0));
        }
    }

    return means;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

ISFTracker::ISFTracker(sptr<AtomGroup> atoms, vector<flt> ks,
                       vector<unsigned long> ns, bool usecom)
    : atoms(atoms), curt(0), usecom(usecom) {
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<unsigned long>::iterator n = ns.begin(); n != ns.end(); ++n) {
        singles.push_back(ISFTracker1(*atoms, *n, ks, com));
    }
};

void ISFTracker::update(Box& box) {
    curt++;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<ISFTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        it->update(box, *atoms, curt, com);
    }
};

void ISFTracker::reset() {
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (vector<ISFTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        it->reset(*atoms, com);
    }
};

vector<vector<vector<cmplx> > > ISFTracker::ISFs() {
    vector<vector<vector<cmplx> > > vals;
    vals.reserve(singles.size());
    for (vector<ISFTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back(it->ISFs());
    }
    return vals;
};

vector<vector<vector<barray<cmplx, NDIM> > > > ISFTracker::ISFxyz() {
    vector<vector<vector<barray<cmplx, NDIM> > > > vals;
    vals.reserve(singles.size());
    for (vector<ISFTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back(it->ISFxyz());
    }
    return vals;
};

vector<flt> ISFTracker::counts() {
    vector<flt> vals;
    vals.reserve(singles.size());
    for (vector<ISFTracker1>::iterator it = singles.begin();
         it != singles.end(); ++it) {
        vals.push_back((flt)it->get_count());
    }
    return vals;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

SmoothLocs::SmoothLocs(sptr<AtomGroup> atoms, Box& box, uint smoothn,
                       uint skipn, bool usecom)
    : atoms(atoms),
      smoothn(smoothn),
      skipn(skipn),
      curlocs(atoms->size(), NDIM),
      numincur(0),
      locs(),
      curt(0),
      usecom(usecom) {
    curlocs.setZero();
    update(box);
};

void SmoothLocs::reset() {
    curlocs.setZero();
    numincur = 0;
    curt = 0;
    locs.clear();
};

void SmoothLocs::update(Box& box) {
    if (curt % skipn != 0) {
        curt++;
        return;
    };

    numincur++;
    bool smooth_time = (numincur >= smoothn);

    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (uint i = 0; i < atoms->size(); ++i) {
        Vec curloc = (*atoms)[i].x - com;
        curlocs.row(i) += curloc;
        if (smooth_time) {
            curlocs.row(i) /= numincur;
        }
    };

    if (smooth_time) {
        locs.push_back(curlocs);
        curlocs.resize(atoms->size(), Eigen::NoChange);
        curlocs.setZero();
        numincur = 0;
    };

    curt++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
RDiffs::RDiffs(sptr<AtomGroup> atoms, unsigned long skip, bool usecom)
    : atoms(atoms),
      pastlocs(atoms->size(), NDIM),
      dists(),
      skip(skip),
      curt(0),
      usecom(usecom) {
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (uint i = 0; i < atoms->size(); ++i) {
        pastlocs.row(i) = (*atoms)[i].x - com;
    }
}

void RDiffs::reset() {
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    for (uint i = 0; i < atoms->size(); ++i) {
        pastlocs.row(i) = (*atoms)[i].x - com;
    }
    dists.clear();
}

void RDiffs::update(Box& box) {
    curt++;
    if (curt < skip) return;
    Vec com = usecom ? atoms->com() : Vec::Zero();
    vector<flt> rdiff = vector<flt>(atoms->size(), 0.0);
    for (uint i = 0; i < atoms->size(); ++i) {
        Vec loc = (*atoms)[i].x - com;
        Vec pastr = pastlocs.row(i);
        rdiff[i] = (loc - pastr).norm();
        pastlocs.row(i) = loc;
    }

    dists.push_back(rdiff);

    curt = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool JammingList::operator<(const JammingList& other) const {
    // return distance_squared < other.distance_squared;
    if (other.distance_squared - distance_squared > 1e-8) return true;
    if (distance_squared - other.distance_squared > 1e-8) return false;
    //~ cout << "\nWithin 1e-8\n";

    uint sz = size();
    uint osz = other.size();
    if (sz < osz) return true;
    if (sz > osz) return false;
    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for (uint i = 0; i < sz; ++i) {
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing sizes...\n";
    //~ cout << "Equal!\n";
    return false;  // consider them equal
};

/* There are two ways of looking at the different arrangements.
 * In both cases, we leave A the same as it was, and rotate / flip / translate
 B.
 * Also in both cases, we wrap A, then subtract off its COM (in an infinite
 box).
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
   * Only try different rotations, but measure distances as Σ (\vec r_{ij} -
 \vec s_{ij})²

 * We'll go with method 3.
*/

bool JammingListRot::operator<(const JammingListRot& other) const {
    // return distance_squared < other.distance_squared;
    if (other.distance_squared - distance_squared > 1e-8) return true;
    if (distance_squared - other.distance_squared > 1e-8) return false;

    uint sz = size();
    uint osz = other.size();
    if (sz < osz) return true;
    if (sz > osz) return false;
    //~ cout << "\nWithin 1e-8. ";

    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for (uint i = 0; i < sz; ++i) {
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing rotations...";
    if (rotation < other.rotation) return true;
    if (rotation > other.rotation) return false;

    //~ cout << " Comparing sizes...";
    if (sz < osz) return true;
    if (sz > osz) return false;
    //~ cout << "Equal!\n";
    return false;  // consider them equal
};

JammingTreeRot::JammingTreeRot(sptr<Box> box,
                               Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& locsA0,
                               Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& locsB0,
                               bool use_rotations, bool use_inversions)
    : box(box), jlists(), A(locsA0), Bs(DIMROTATIONS, locsB0) {
    for (uint rot = 0; rot < DIMROTATIONS; ++rot) {
        if (!use_rotations and (rot % DIMROTATIONS != 0)) continue;
        if (!use_inversions and (rot >= DIMROTATIONS)) continue;
        for (uint i = 0; i < locsB0.rows(); ++i) {
            Vec loc = locsB0.row(i);
            Bs[rot].row(i) = rotate_flip(loc, rot);
        }
        if (locsA0.rows() <= locsB0.rows())
            jlists.push_back(JammingListRot(rot));
        //~ cout << "Created, now size " << jlists.size() << endl;
    }
};

flt JammingTreeRot::distance(JammingListRot& jlist) {
    flt dist = 0;
    uint rot = jlist.rotation;
    for (uint i = 1; i < jlist.size(); ++i) {
        uint si = jlist.assigned[i];
        for (uint j = 0; j < i; ++j) {
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A.row(i), A.row(j));
            Vec sij = box->diff(Bs[rot].row(si), Bs[rot].row(sj));
            dist += box->diff(rij, sij).squaredNorm();
        }
    }
    return dist / ((flt)jlist.assigned.size());
};

list<JammingListRot> JammingTreeRot::expand(JammingListRot curjlist) {
    vector<uint>& curlist = curjlist.assigned;
    list<JammingListRot> newlists = list<JammingListRot>();
    if (curlist.size() >= (uint)A.rows()) {
        return newlists;
    }

    uint N = (uint)Bs[curjlist.rotation].rows();
    for (uint i = 0; i < N; ++i) {
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        // if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;

        JammingListRot newjlist = JammingListRot(curjlist, i, 0);
        newjlist.distance_squared = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool JammingTreeRot::expand() {
    JammingListRot curjlist = jlists.front();
    list<JammingListRot> newlists = expand(curjlist);

    if (newlists.empty()) {
        //~ cout << "No lists made\n";
        return false;
    }
    //~ cout << "Have " << newlists.size() << "\n";
    newlists.sort();
    //~ cout << "Sorted.\n";
    jlists.pop_front();
    //~ cout << "Popped.\n";
    jlists.merge(newlists);
    //~ cout << "Merged to size " << jlists.size() << "best dist now " <<
    // jlists.front().distance_squared << "\n";
    return true;
};

JammingTreeBD::JammingTreeBD(sptr<Box> box,
                             Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                             Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B,
                             uint cutoffA, uint cutoffB, bool use_rotations,
                             bool use_inversions)
    : JammingTreeRot(box, A, B, use_rotations, use_inversions),
      cutoff1(cutoffA),
      cutoff2(cutoffB) {
    if (cutoffA > cutoffB) {
        jlists.clear();
    }
    if (A.rows() - cutoffA > B.rows() - cutoffB) {
        jlists.clear();
    }
};

list<JammingListRot> JammingTreeBD::expand(JammingListRot curjlist) {
    vector<uint>& curlist = curjlist.assigned;
    list<JammingListRot> newlists = list<JammingListRot>();
    if (curlist.size() >= (uint)A.rows()) {
        return newlists;
    }

    uint N = (uint)Bs[curjlist.rotation].rows();
    uint start = 0, end = cutoff2;
    if (curlist.size() >= cutoff1) {
        start = cutoff2;
        end = N;
    }
    for (uint i = start; i < end; ++i) {
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        // if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;

        JammingListRot newjlist = JammingListRot(curjlist, i, 0);
        newjlist.distance_squared = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool JammingTreeBD::expand() {
    JammingListRot curjlist = jlists.front();
    list<JammingListRot> newlists = expand(curjlist);

    if (newlists.empty()) return false;
    newlists.sort();
    jlists.pop_front();
    jlists.merge(newlists);
    return true;
};

Eigen::Matrix<flt, Eigen::Dynamic, NDIM> JammingTreeRot::locations_B(
    JammingListRot jlist) {
    uint rot = jlist.rotation;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs =
        Eigen::Matrix<flt, Eigen::Dynamic, NDIM>(jlist.size(), NDIM);

    uint N = jlist.size();
    for (uint i = 0; i < N; ++i) {
        uint si = jlist.assigned[i];
        locs.row(i) = A.row(i);
        for (uint j = 0; j < N; ++j) {
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A.row(i), A.row(j));
            Vec sij = box->diff(Bs[rot].row(si), Bs[rot].row(sj));
            locs.row(i) -= box->diff(rij, sij) / N;
        }
    }
    return locs;
};

Eigen::Matrix<flt, Eigen::Dynamic, NDIM> JammingTreeRot::locations_A(
    JammingListRot jlist) {
    uint rot = jlist.rotation;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs =
        Eigen::Matrix<flt, Eigen::Dynamic, NDIM>(Bs[rot].rows(), NDIM);

    uint N = jlist.size();
    for (uint i = 0; i < N; ++i) {
        uint si = jlist.assigned[i];
        locs.row(si) = Bs[rot].row(si);
        for (uint j = 0; j < N; ++j) {
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A.row(i), A.row(j));
            Vec sij = box->diff(Bs[rot].row(si), Bs[rot].row(sj));
            locs.row(si) += box->diff(rij, sij) / N;
        }

        // this is an inverse rotateflip

        Vec loc = locs.row(si);
        locs.row(si) = rotate_flip_inv(loc, rot);
    }
    return locs;
};

Vec JammingTreeRot::straight_diff(
    Box& bx, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& As,
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& Bs) {
    uint N = (uint)As.rows();
    if (Bs.rows() != N) {
        throw std::runtime_error("As and Bs are not of the same shape");
    }

    Vec loc = Vec::Zero();
    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < N; ++j) {
            Vec rij = bx.diff(As.row(i), As.row(j));
            Vec sij = bx.diff(Bs.row(i), Bs.row(j));
            loc += bx.diff(rij, sij);
        }
    }
    return loc / N;
};

flt JammingTreeRot::straight_distsq(
    Box& bx, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& As,
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& Bs) {
    long int N = As.rows();
    if (Bs.rows() != N) return NAN;

    flt dist = 0;
    for (uint i = 0; i < N; ++i) {
        for (uint j = 0; j < N; ++j) {
            Vec rij = bx.diff(As.row(i), As.row(j));
            Vec sij = bx.diff(Bs.row(i), Bs.row(j));
            dist += bx.diff(rij, sij).squaredNorm();
        }
    }
    return dist / ((flt)N);
};

////////////////////////////////////////////////////////////////////////////////////////////////////

void CNodePath::add(CNode node, OriginBox& box) {
    if (size() > 0) {
        distance += box.diff(node.x, nodes.back().x);
    }
    nodes.push_back(node);
};

void Connectivity::add_edge(CNode node1, CNode node2) {
    if (node1.n < 0 || node2.n < 0) {
        throw invalid_argument("Got node with negative n");
    }
    nodes.insert(node1);
    nodes.insert(node2);

    neighbors[node1.n].push_back(node2);
    std::sort(neighbors[node1.n].begin(), neighbors[node1.n].end());
    neighbors[node2.n].push_back(node1);
    std::sort(neighbors[node2.n].begin(), neighbors[node2.n].end());
};

barray<bool, NDIM> Connectivity::nonzero(Vec diff_vec) {
    barray<bool, NDIM> nonzeros;
    Vec half_shape = box->box_shape() / 2;
    for (uint i = 0; i < NDIM; i++) {
        nonzeros[i] = (abs(diff_vec[i]) > half_shape[i]);
    }
    return nonzeros;
};

void Connectivity::add(Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs,
                       vector<flt> diameters) {
    vector<CNode> cnodes;

    uint n = 0;
    for (uint i = 0; i < locs.rows(); i++) {
        CNode cn = CNode((int)n, locs.row(i));

        uint n2 = 0;
        for (vector<CNode>::const_iterator cit = cnodes.begin();
             cit != cnodes.end(); cit++) {
            Vec dx = box->diff(cn.x, cit->x);
            if (dx.norm() <= (diameters[n] + diameters[n2]) / 2.0) {
                neighbors[n].push_back(*cit);
                std::sort(neighbors[n].begin(), neighbors[n].end());
                neighbors[n2].push_back(cn);
                std::sort(neighbors[n2].begin(), neighbors[n2].end());
            }
            n2++;
        }
        nodes.insert(cn);
        cnodes.push_back(cn);

        n++;
    }
};

CNodePath Connectivity::make_cycle(CNodePath forward, CNodePath backward) {
    // forward and backward should both start and end in the same place
    // and vice versa

    if (forward.size() == 0)
        throw invalid_argument(
            "make_cycle got an empty forward path. This should not happen.");
    else if (backward.size() == 0)
        throw invalid_argument(
            "make_cycle got an empty forward path. This should not happen.");

    if (forward.nodes.front() != backward.nodes.front())
        throw invalid_argument(
            "make_cycle: forward.front() != backward.front()");
    if (backward.nodes.back() != forward.nodes.back())
        throw invalid_argument("make_cycle: backward.back() != forward.back()");

    CNodePath cycle = forward;

    for (vector<CNode>::reverse_iterator it = ++(backward.nodes.rbegin());
         it != backward.nodes.rend(); it++) {
        cycle.add(*it, *box);
    }
    return cycle;
};

map<uint, CNodePath> Connectivity::circular_from(CNode node, set<uint>& visited,
                                                 bool check_all) {
    map<uint, CNodePath>
        full_paths;  // complete roots around the box. key is DIMENSION
    map<uint, CNodePath>
        prev_paths;  // other ways we've found to get to any node
    queue<CNodePath> paths;

    CNodePath path0 = CNodePath(node);
    paths.push(path0);
    prev_paths[node.n] = path0;
    visited.insert(node.n);

    while (!paths.empty()) {
        CNodePath path = paths.front();
        paths.pop();
        assert(path.size() > 0);
        CNode lastnode = path.nodes.back();
        vector<CNode>& nextnodes = neighbors[lastnode.n];
        for (vector<CNode>::iterator it = nextnodes.begin();
             it < nextnodes.end(); ++it) {
            if (it->n < node.n) continue;
            CNode nextnode = *it;
            CNodePath newpath = CNodePath(path, nextnode, *box);
            map<uint, CNodePath>::iterator found_path_it =
                prev_paths.find(nextnode.n);
            if (found_path_it != prev_paths.end()) {
                // We've been to this node before.
                // But have we found a circle (around the box), or just another
                // route?
                bool found_nonzero = false;
                CNodePath& found_path = found_path_it->second;

                Vec pathdiff = found_path.distance - newpath.distance;
                barray<bool, NDIM> nonzeros = nonzero(pathdiff);
                for (uint i = 0; i < NDIM; i++) {
                    if (nonzeros[i]) {
                        found_nonzero = true;
                        // We've found a cycle around the whole box!
                        map<uint, CNodePath>::iterator found_cycle =
                            full_paths.find(i);
                        CNodePath cycle_path = make_cycle(found_path, newpath);
                        if (found_cycle == full_paths.end()) {
                            // And its along a dimension we've never found
                            // before
                            // Note that now both *found_path and newpath
                            // begin with CNode "node" and end with "nextnode"

                            full_paths[i] = cycle_path;
                            // Have we found enough paths?
                            // If so, we're done. We're not trying to find the
                            // "best" cycle,
                            // just any cycle (or any NDIM cycles if check_all)
                            if (!check_all)
                                return full_paths;
                            else if (full_paths.size() >= NDIM)
                                return full_paths;
                        } else {
                            // Already found a cycle in this dimension, maybe
                            // replace it
                            if (full_paths[i].nodes.size() >
                                cycle_path.nodes.size())
                                full_paths[i] = cycle_path;
                        }
                    }
                }

                if (!found_nonzero) {
                    // Different route to the same node, but through the same
                    // box.
                    // If this is a shorter route, might as well hold onto it...
                    if (newpath.size() < found_path.size()) {
                        prev_paths[nextnode.n] = newpath;
                    }
                    continue;
                }
            } else {
                // we've never been to this node before
                visited.insert(nextnode.n);
                prev_paths[nextnode.n] = newpath;
                paths.push(newpath);
                continue;
            }
            // Now same path, new neighbor to add to it
        }
        // We're done, move on to the next incomplete path
    };

    return full_paths;
};

map<uint, CNodePath> Connectivity::find_percolation(bool check_all) {
    map<uint, CNodePath> full_paths;
    set<uint> visited;
    for (set<CNode>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
        if (visited.count(it->n) > 0) continue;
        // we've already been to this node, no need to start there

        map<uint, CNodePath> new_paths = circular_from(*it, visited, check_all);
        for (map<uint, CNodePath>::iterator fit = new_paths.begin();
             fit != new_paths.end(); ++fit) {
            uint dim = fit->first;
            CNodePath& new_path = fit->second;
            if (full_paths.find(dim) == full_paths.end()) {
                // never been along this dimension before
                full_paths[dim] = new_path;
            } else {
                // found one along this dimension before. If the new one is
                // shorter, might as
                // well take it instead
                CNodePath& old_path = full_paths[dim];
                if (new_path.nodes.size() < old_path.nodes.size()) {
                    full_paths[dim] = new_path;
                }
            }
        }
        if (!check_all and !full_paths.empty())
            return full_paths;
        else if (full_paths.size() >= NDIM)
            return full_paths;
    }
    return full_paths;
};
