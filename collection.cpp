#include "collection.hpp"

collection::collection(Box *box, vector<atomgroup*> gs, vector<interaction*> is,
              vector<statetracker*> ts, vector<constraint*> cs)
        : box(box), groups(gs), interactions(is), trackers(ts), constraints(cs),
                atoms(gs){
    update_trackers();
    setForces(true);
    update_constraints();
    update_trackers();
}

void collection::scaleVs(flt scaleby){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.size(); i++){
            g[i].v *= scaleby;
        }
    }
}

void collection::scaleVelocitiesT(flt T){
    flt t = temp();
    flt scaleby = sqrt(T/t);
    scaleVs(scaleby);
}

void collection::scaleVelocitiesE(flt E){
    flt E0 = energy();
    flt k0 = kinetic();
    flt goalkinetic = k0 + (E - E0);
    flt scaleby = sqrt(goalkinetic/k0);
    scaleVs(scaleby);
}

void collection::update_trackers(){
    vector<statetracker*>::iterator git;
    for(git = trackers.begin(); git<trackers.end(); git++){
        (*git)->update(box);
    }
}

void collection::update_constraints(){
    vector<constraint*>::iterator git;
    for(git = constraints.begin(); git<constraints.end(); git++){
        (*git)->apply(box);
    }
}

flt collection::kinetic(){
    flt E=0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        E+= group.kinetic();
    }
    return E;
}

flt collection::virial(){
    flt E = 0;
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        E += (*it)->pressure(box);
    }
    return E;
}

flt collection::pressure(){
    flt V = box->V();
    //~ if(isnan(V) or isinf(V)) return NAN;
    //~ int ndof = 0, nc = 0;
    //~ vector<constraint*>::iterator cit;
    //~ for(cit = constraints.begin(); cit<constraints.end(); cit++){
        //~ nc += (*cit)->ndof();
    //~ }
    //~ vector<atomgroup*>::iterator git;
    //~ for(git = groups.begin(); git<groups.end(); git++){
        //~ atomgroup &group = **git;
        //~ ndof += 3*(group.size());
    //~ }
    
    
    flt E = 2.0 * kinetic();// * (ndof - nc - 3) / ndof;
    //flt E = (ndof - nc) * temp();
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        E += (*it)->pressure(box);
        assert(not isnan(E));
    }
    return E / V / 3.0;
}

flt collection::potentialenergy(){
    flt E=0;
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        E += inter.energy(box);
        //~ cout << "potential energy: " << E << endl;
        assert(not isnan(E));
    }
    return E;
}

flt collection::energy(){
    flt E = potentialenergy() + kinetic();
    assert(not isnan(E));
    return E;
};

flt collection::dof(){
    int ndof = 0;
    vector<constraint*>::iterator cit;
    for(cit = constraints.begin(); cit<constraints.end(); cit++){
        ndof -= (*cit)->ndof();
    }
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        ndof += 3*(group.size());
    }
    
    return ndof;
}

flt collection::temp(){
    flt totatoms = 0;
    flt totkinetic = 0;
    Vec v = comv();
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        totkinetic += group.kinetic(v);
        //~ cout << "group K: " << group.kinetic(v) << ", totatoms: " << group.size() << "\n";
        totatoms += group.size();
    }
    
    int ndof = dof();
    //~ cout << "K: " << totkinetic << ", totatoms: " << totatoms << "\n";
    return totkinetic * 2 / (ndof-3);
}

//~ Vec collection::com(){
    //~ vector<atomgroup*>::iterator git;
    //~ flt mass = 0;
    //~ Vec totcom = Vec();
    //~ for(git = groups.begin(); git<groups.end(); git++){
        //~ atomgroup &g = **git;
        //~ for(uint i = 0; i<g.size(); i++){
            //~ flt curmass = g.getmass(i);
            //~ mass += curmass;
            //~ totcom += g[i].x * curmass;
        //~ }
    //~ }
    //~ return totcom / mass;
//~ }
//~ 
//~ Vec collection::comv(){
    //~ vector<atomgroup*>::iterator git;
    //~ flt mass = 0;
    //~ Vec totcom = Vec();
    //~ for(git = groups.begin(); git<groups.end(); git++){
        //~ atomgroup &g = **git;
        //~ for(uint i = 0; i<g.size(); i++){
            //~ flt curmass = g.getmass(i);
            //~ mass += curmass;
            //~ totcom += g[i].v * curmass;
        //~ }
    //~ }
    //~ return totcom / mass;
//~ }

flt collection::gyradius(){
    vector<atomgroup*>::iterator git;
    Vec avgr = Vec();
    flt N = 0;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.size(); i++){
            avgr += g[i].x;
            N++;
        }
    }
    avgr /= N; // now avgr is the average location, akin to c.o.m.
    flt Rgsq = 0;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.size(); i++) Rgsq += (g[i].x - avgr).sq();
    }
    
    return sqrt(Rgsq/N);
}

void collection::setForces(bool seta){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        group.resetForces();
    }
    
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        inter.setForces(box);
    }
    if(!seta) return;
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
        }
    }
}

collectionSol::collectionSol(Box *box, const flt dt, const flt damp,const flt T,
        vector<atomgroup*> groups,vector<interaction*> interactions,
        vector<statetracker*> trackers, vector<constraint*> constraints) :
            collection(box, groups, interactions, trackers, constraints), 
            dt(dt), damping(damp), desT(T){
    setCs();
    update_trackers();
    setForces(true);
    update_constraints();
    update_trackers();
};

void collectionSol::setCs(){
    if(damping == 0){
        c0 = 1; c1 = 1; c2 = .5;
        sigmar = 0; sigmav = 0; corr = 1;
        gauss.set(sigmar, sigmav, corr);
        return;
    }
    // from Allen and Tildesley, 262
    flt dampdt = damping * dt;
    c0 = exp(-dampdt);
    c1 = (-expm1(-dampdt))/dampdt;
    c2 = (1-c1)/dampdt;
    
    // note that for sigmar, sigmav, and corr, we have taken the
    // k_B*T/m out and put it in the timestep place
    // note that these are unitless sigmas and correlations,
    // so sigmar has a dt*sqrt(k T/m) removed, and sigmav has a sqrt(k T/m)
    // removed
    if(dampdt > 1e-4)
        sigmar = sqrt((1/dampdt) * (
            2-(-4*expm1(-dampdt) + expm1(-2*dampdt))/dampdt
            ));
    else
        sigmar = sqrt(2*dampdt/3 - dampdt*dampdt/2 + 7*dampdt*dampdt*dampdt/30);
    //~ cout << "setCs: " << -4*expm1(-dampdt) << ',' << expm1(-2*dampdt)
         //~ << ", (...): " << (-4*expm1(-dampdt) + expm1(-2*dampdt))/dampdt
         //~ << ", sigmar:" << sigmar << '\n';
    sigmav = sqrt(-expm1(-2*dampdt));
    flt exdpdt = (-expm1(-dampdt));
    corr = exdpdt*exdpdt/dampdt/sigmar/sigmav;
    //~ cout << "vals: " << c0 << ',' << c1 << ',' << c2 << "; T=" <<desT
         //~ << ", sigmas: (" << sigmar << ',' << sigmav << ',' << corr << ')'
         //~ << ", dt=" << dt << ", damping=" << damping << "\n";
    gauss.set(sigmar, sigmav, corr);
}

void collectionSol::timestep(){
    vector<atomgroup*>::iterator git;
    // From Allen and Tildesley 263, Verlet-like, our equations are
    // r(t+dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + drG
    // v(t+dt) = c0 v(t) + (c1 - c2) dt a(t) + c2 dt a(t+dt) + dvG
    // I split this:
    // r(t+dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + drG
    // vp(t) = c0 v(t) + (c1 - c2) dt a(t) + dvG
    // (set forces)
    // v(t+dt) = vp(t) + c2 dt a(t+dt)
    
    //Step 1: set atom.x = r(t+dt) and atom.v = vp(t)
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            flt v0 = sqrt(desT/m.getmass(i));
            flt r0 = dt * v0;
            VecPair vecpair;
            if (damping > 0) vecpair = gauss.genVecs();
            else vecpair[0] = vecpair[1] = Vec();
            // vecpair[0] is drG, and vecpair[1] is dvG
            m[i].x += m[i].v * (c1 * dt) + m[i].a * (c2*dt*dt) + vecpair[0]*r0;
            m[i].v = m[i].v*c0 + m[i].a * (dt*(c1-c2)) + vecpair[1]*v0;
            //~ if(i==0 and git == groups.begin()) cout 
                //~ << "drG: " << vecpair[0].mag() 
                //~ << ", dvG: " << vecpair[1].mag()
                //~ << "\n";
        }
    }
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
        }
    }
    
    // And finish m[i].v
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].v += m[i].a * (dt*c2);
        }
    }
    
    update_trackers();
};

collectionSolHT::collectionSolHT(Box *box, const flt dt, const flt damp,const flt T,
        vector<atomgroup*> groups,vector<interaction*> interactions,
        vector<statetracker*> trackers, vector<constraint*> constraints) :
            collection(box, groups, interactions, trackers, constraints), 
            dt(dt), damping(damp), desT(T), gauss(sqrt(2.0*desT*damping/dt)){
    setGauss();
    update_trackers();
    setForces(true);
    update_constraints();
    update_trackers();
};

void collectionSolHT::setGauss(){
    gauss.set(sqrt(2.0*desT*damping/dt)); // TODO: Is that correct?
}

void collectionSolHT::timestep(){
    //Honeycutt and Thirumalai
    // note atom.a is the second derivative of position, d²r/dt², and
    // includes the random term / damping term.
    // atom.f includes only the forces from other particles.
    
    vector<atomgroup*>::iterator git;
    flt keepv = 1-(damping*dt);
    flt xpartfromv = dt - (dt*dt*damping/2);
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            // step 1: get new x
            m[i].x += (m[i].v * xpartfromv) + (m[i].a * (.5 * dt*dt));
            
            // step 2: make v1 (intermediate v)
            m[i].v = m[i].v*keepv + m[i].a*(dt/2);
        }
    }
    
    // step 3: set forces
    setForces();
    
    // step 4: intermediate acceleration
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            Vec g = gauss.generate();
            //~ cout << "g " << g << "\n";
            m[i].a = (m[i].f + g) / m.getmass(i);
            
            // step 5: finish v
            m[i].v += m[i].a * (dt/2);
        }
    }
    update_trackers();
}

void collectionVerlet::timestep(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].x += m[i].v * dt + m[i].a * (dt*dt/2);
            m[i].v += m[i].a * (dt/2);
        }
    }
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
            // And finish m[i].v
            m[i].v += m[i].a * (dt/2);
        }
    }
    
    update_trackers();
    
    return;
}

void collectionNoseHoover::timestep(){
    // From Corey's notes. Note that xi and lns are now coordinates, with
    // lns a variable.
    // K = <sum over i> m v^2_i
    // L = d N   // that's n_{degrees of freedom}
    // Q is our damping coefficient.
    // * r(t+dt) = r(t) + dt v(t) + dt^2 (a(t) - xi(t) v(t))
    // * lns(t+dt) = lns(t) + xi(t) dt + (T(t) - L T0) dt^2/(2 Q)
    // v(t+dt) = v(t) + (a(t) - xi(t+dt) v(t+dt) + F(t) - xi(t) v(t)) dt/2
    // xi(t+dt) = xi(t) + (T(t+dt) - L T0 + T(t) - L T0) dt/(2 Q)
    
    // I modified v -> y:
    // Definition: y(t+dt) = (1 + xi(t+dt) dt/2) v(t) // that's \vec y_i
    // Definition: Ky = <sum over i> m y^2_i
    // * y(t+dt) = v(t) + (a(t+dt) - F(t) - xi(t) v(t)) dt/2 // put this in v
    // * z0 = xi(t) + (T(t) - 2L T)/ (2 Q) // temporary
    // * z1 = (Ty(t+dt)) dt^3/(8 Q) // temporary
    
    // Now we know
    // 0 = xi^3 + xi^2 (dt - z0) + xi (dt^2/4 - dt z0) - (z0 dt^2/4 z0 + z1)
    // so
    // * xi(t+dt) = solveCubic(dt - z0, dt^2/4 - dt z0, z0 dt^2/4 z0 + z1)
    
    flt ndof = dof();
    flt Kt = 2*kinetic();
    
    //Step 1: set atom.x = r(t+dt)
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            Vec Ftilde = (m[i].a - (m[i].v*xi));
            m[i].x += m[i].v * dt + Ftilde * (dt*dt/2);
            // first half of y
            m[i].v += Ftilde * (dt/2);
        }
    }
    
    lns += xi * dt + (Kt - ndof*T) * (dt*dt/2/Q);
    
    
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
            // And finish y
            m[i].v += m[i].a * (dt/2);
        }
    }
    
    // Now we solve for xi
    flt Ky = 2*kinetic();
    flt z0 = xi + (Kt - 2*ndof*T)*(dt/2/Q);
    flt z1 = Ky*2/dt/Q;
    
    //~ flt oldxi = xi;
    xi = solveCubic(4/dt - z0, 4/dt/dt - 4*z0/dt, -(z0*4/dt/dt) - z1, xi);
    
    
    // And finish v
    flt ytov = 1 + xi*dt/2;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].v /= ytov;
        }
    }
    
    //~ flt Kt2 = 2*kinetic();
    //~ flt xicheck = oldxi + (Kt2 - ndof*T + Kt - ndof*T)*(dt/2/Q);
    //~ flt xicheck2 = z0 + z1 / ((2/dt) + xi) / ((2/dt) + xi);
    //~ if (abs(xi-xicheck) > 1e-4){
        //~ printf("new xi: %10.5f  xicheck: %10.5f  zs: %10.5f\n", xi, xicheck, xicheck2);
    //~ }
    
    update_trackers();
}

flt collectionNoseHoover::Hamiltonian(){
    flt H = (kinetic() + potentialenergy() + 
                (xi*xi*Q/2) + (dof() * lns * T));
    assert(not isnan(H));
    return H;
}

flt collectionGaussianT::setxi(){
    flt xiNumerator = 0, xiDenominator = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            flt m = g.getmass(i);
            xiNumerator += g[i].f.dot(g[i].v) * m;
            xiDenominator += (g[i].v.sq()) * m*m;
        }
    }
    xi = xiNumerator / xiDenominator;
    return xi;
}

void collectionGaussianT::setForces(bool seta, bool shouldwesetxi){
    collection::setForces(seta);
    if(shouldwesetxi) setxi();
}

void collectionGaussianT::timestep(){
    //~ flt oldT = temp();
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].x += m[i].v * dt + m[i].a * (dt*dt/2);
            m[i].v += (m[i].a - (m[i].v*xi)) * (dt/2);
        }
    }
    // Now we set forces and accelerations
    setForces(false, false);
    update_constraints();
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
            // And finish y
            m[i].v += m[i].a * (dt/2);
        }
    }
    
    // now we get z and set xi
    //~ flt oldxi = xi;
    flt z = setxi();
    xi = z/(1 - z*dt/2);
    
    // And finish v
    flt ytov = 1 + xi*dt/2;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].v /= ytov;
        }
    }
    
    //~ flt midxi = xi;
    //~ setxi();
    //~ if(abs(midxi-xi) > 1e-9){
        //~ printf("xi change of %7.3g (%7.3f -> %7.3f)\n", 
            //~ abs(midxi-xi), midxi, xi);
    //~ }
    
    //~ flt T = temp();
    //~ printf("Temperature (%6.3f) %7.3g -> %7.3g, xi (%6.3f) %7.3g -> %7.3g\n",
            //~ T/oldT, oldT, T, xi/oldxi,oldxi,xi);
    
    update_trackers();
    return;
}


void collectionGear3A::timestep(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].x += m[i].v * dt + m[i].a * (dt*dt/2);
            m[i].v += m[i].a * dt;
        }
    }
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    
    // Make correction
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            Vec a = m[i].f / m.getmass(i);
            Vec correction = a - m[i].a;
            m[i].v += correction * (dt/2);
            m[i].a = a;
        }
    }
    
    update_trackers();
    
    return;
}

void collectionGear4A::timestep(){
    uint n=0;
    //~ cout << "Gear 4A timestep, ncorrec = " << ncorrec << "\n";
    
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        //~ cout << "Gear 4A timestep loop 1, g.size() = " << g.size() << "\n";
        for(uint i=0; i<g.size(); i++){
            //~ cout << "Gear 4A timestep loop 2, i = " << i << "\n";
            g[i].x += g[i].v * dt + g[i].a * (dt*dt/2) + bs[n] * (dt*dt*dt/6);
            g[i].v += g[i].a * dt + bs[n] * (dt*dt/2);
            g[i].a += bs[n] * dt;
            n++;
        }
    }
    
    //~ cout << "Gear 4A timestep 2, ncorrec = " << ncorrec << "\n";
    // Now we set forces and accelerations
    for(uint m=0; m<ncorrec; m++){
        //~ cout << "Gear 4A setting forces, m = " << m << ", ncorrec = " << ncorrec << "\n";
        setForces(false);
        update_constraints();
    
        // Make correction
        n=0;
        for(git = groups.begin(); git<groups.end(); git++){
            atomgroup &g = **git;
            for(uint i=0; i<g.size(); i++){
                Vec a = g[i].f / g.getmass(i);
                Vec correction = a - g[i].a;
                g[i].x += correction * (dt*dt/12);
                g[i].v += correction * (5*dt/12);
                g[i].a = a;
                bs[n] += correction / dt;
                n++;
            }
        }
    }
    update_trackers();
    return;
}

void collectionGear5A::timestep(){
    uint n=0;
    flt dt2 = dt*dt/2;
    flt dt3 = dt*dt2/3;
    flt dt4 = dt*dt3/4;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            g[i].x += g[i].v * dt + g[i].a * dt2 + bs[n] * dt3 + cs[n] * dt4;
            g[i].v += g[i].a * dt + bs[n] * dt2 + cs[n] * dt3;
            g[i].a += bs[n] * dt + cs[n] * dt2;
            bs[n] += cs[n] * dt;
            n++;
        }
    }
    flt c0 = 19*dt*dt/240, c1 = 3*dt/8;
    //flt c2 = 1
    flt c3 = 3 / (2 * dt), c4 = 1/(dt*dt);
    
    for(uint m=0; m<ncorrec; m++){
        // Now we set forces and accelerations
        setForces(false);
        update_constraints();
        
        // Make correction
        n=0;
        for(git = groups.begin(); git<groups.end(); git++){
            atomgroup &g = **git;
            for(uint i=0; i<g.size(); i++){
                Vec a = g[i].f / g.getmass(i);
                Vec correction = a - g[i].a;
                g[i].x += correction * c0;
                g[i].v += correction * c1;
                g[i].a = a;
                bs[n] += correction * c3;
                cs[n] += correction * c4;
                n++;
            }
        }
    }
    
    update_trackers();
    return;
}

void collectionGear6A::timestep(){
    //~ cout << "Timestep\n";
    flt dt2 = dt*dt/2;
    flt dt3 = dt*dt2/3;
    flt dt4 = dt*dt3/4;
    flt dt5 = dt*dt4/5;
    
    uint n=0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            //~ assert(n < bs.size());
            //~ assert(n < cs.size());
            //~ assert(n < ds.size());
            g[i].x += g[i].v * dt + g[i].a * dt2 + bs[n] * dt3 + cs[n] * dt4 + ds[n] * dt5;
            g[i].v += g[i].a * dt + bs[n] * dt2 + cs[n] * dt3 + ds[n] * dt4;
            g[i].a += bs[n] * dt + cs[n] * dt2 + ds[n] * dt3;
            bs[n] += cs[n] * dt + ds[n] * dt2;
            cs[n] += ds[n] * dt;
            n++;
        }
    }
    
    //~ atomgroup &gbegin = **(groups.begin());
    //~ printf("predic = {%.15f, %.15f, %.15f, %.15f, %.15f, %.15f}\n",
        //~ gbegin[0].x[0], gbegin[0].v[0], gbegin[0].a[0],
            //~ bs[0][0], cs[0][0], ds[0][0]);
    //~ printf("predicscaled = {%.15f, %.15f, %.15f, %.15f, %.15f, %.15f}\n",
        //~ gbegin[0].x[0], gbegin[0].v[0]*dt, gbegin[0].a[0]*dt2,
            //~ bs[0][0]*dt3, cs[0][0]*dt4, ds[0][0]*dt5);
    //~ cout << "Prediction: " 
         //~ << gbegin[0].x[0] << ", "
         //~ << gbegin[0].v[0] << ", "
         //~ << gbegin[0].a[0] << ", "
         //~ << bs[0][0] << ", "
         //~ << cs[0][0] << ", "
         //~ << ds[0][0] << "\n";
    
    
    flt c0 = 3*dt*dt/40, c1 = 251*dt/720;
    //flt c2 = 1
    flt c3 = 11 / (6*dt), c4 = 2/(dt*dt), c5 = 1/(dt*dt*dt);
    
    for(uint m=0; m<ncorrec; m++){
        // Now we set forces and accelerations
        setForces(false);
        update_constraints();
        
        //~ flt correctax = (gbegin[0].f / gbegin.getmass(0))[0];
        //~ printf("correcta = %.15f\n", correctax);
        //~ printf("correctascaled = %.15f\n", correctax*dt2);
    
        // Make correction
        n=0;
        for(git = groups.begin(); git<groups.end(); git++){
            atomgroup &g = **git;
            for(uint i=0; i<g.size(); i++){
                //~ assert(n < bs.size());
                //~ assert(n < cs.size());
                //~ assert(n < ds.size());                
                Vec a = g[i].f / g.getmass(i);
                Vec correction = a - g[i].a;
                g[i].x += correction * c0;
                g[i].v += correction * c1;
                g[i].a = a;
                bs[n] += correction * c3;
                cs[n] += correction * c4;
                ds[n] += correction * c5;
                //~ if(correction.mag()*dt2 > .05) cout << "|c[" << i << "]| = " << correction.mag()*dt2 << '\n';
                //~ if(a.mag()*dt2 > .1) cout << "|a[" << i << "]| = " << a.mag()*dt2 << '\n';
                //~ if(bs[n].mag()*dt3 > .1) cout << "|bs[" << n << "]| = " << bs[n].mag()*dt3 << '\n';
                //~ if(cs[n].mag()*dt4 > .1) cout << "|cs[" << n << "]| = " << cs[n].mag()*dt4 << '\n';
                //~ if(ds[n].mag()*dt5 > .1) cout << "|ds[" << n << "]| = " << ds[n].mag()*dt5 << '\n';
                //~ assert(correction.mag()*dt2 < 1);
                //~ assert(a.mag()*dt2 < 3);
                //~ assert(bs[n].mag()*dt3 < 3);
                //~ assert(cs[n].mag()*dt4 < 3);
                //~ assert(ds[n].mag()*dt5 < 3);
                n++;
            }
        }
        //~ printf("final = {%.15f, %.15f, %.15f, %.15f, %.15f, %.15f}\n",
        //~ gbegin[0].x[0], gbegin[0].v[0], gbegin[0].a[0],
            //~ bs[0][0], cs[0][0], ds[0][0]);
        //~ printf("finalscaled = {%.15f, %.15f, %.15f, %.15f, %.15f, %.15f}\n",
        //~ gbegin[0].x[0], gbegin[0].v[0]*dt, gbegin[0].a[0]*dt2,
            //~ bs[0][0]*dt3, cs[0][0]*dt4, ds[0][0]*dt5);
    
        //~ cout << "Final:      " 
         //~ << gbegin[0].x[0] << ", "
         //~ << gbegin[0].v[0] << ", "
         //~ << gbegin[0].a[0] << ", "
         //~ << bs[0][0] << ", "
         //~ << cs[0][0] << ", "
         //~ << ds[0][0] << "\n";
    
    }
    
    update_trackers();
    return;
}

atomid atomvecRK4::get_id(atom* atm){
    atomRK4* a = (atomRK4*) atm;
    uint n = a - atoms;
    if (n >= sz or a < atoms) return atomid();
    return atomid(atoms + n, n);
};

void collectionRK4::timestep(){
    //~ atom &atm0 = *((*(groups.begin()))->get(0));
    //~ atomRK4 &a0 = (atomRK4 &) atm0;
    //~ Vec xold = a0.x, vold = a0.v;
    
    
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        //atomgroup *agroup = *git;
        //atomvecRK4 &m = *aRK;
        //atomvecRK4 &m = *((atomvecRK4*) *git);
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            atomRK4 &a = *((atomRK4 *) m.get(i));
            a.Kxa = a.v * dt;
            a.Kva = a.f * (dt / m.getmass(i));
            a.x += a.Kxa / 2;
            a.v += a.Kva / 2;
            //~ if(i == 0) cout << "f " << a.f;
        }
    }
    //~ cout << (xold + (a0.Kxa / 2) - a0.x) << "    "
         //~ << (vold + (a0.Kva / 2) - a0.v) << '\n';
    
    // x = x0 + Kxa/2
    // v = v0 + Kva/2
    
    setForces(false);
    update_constraints();
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            atomRK4 &a = *((atomRK4 *) m.get(i));
            a.Kxb = a.v * dt;
            a.Kvb = a.f * (dt / m.getmass(i));
            a.x += a.Kxb / 2 - a.Kxa / 2;
            a.v += a.Kvb / 2 - a.Kva / 2;
            //~ if(i == 0) cout << " " << a.f;
        }
    }
    
    //~ cout << (xold + (a0.Kxb / 2) - a0.x) << "    "
         //~ << (vold + (a0.Kvb / 2) - a0.v) << '\n';
    
    // x = x0 + Kxb/2
    // v = v0 + Kvb/2
    
    setForces(false);
    update_constraints();
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            atomRK4 &a = *((atomRK4 *) m.get(i));
            a.Kxc = a.v * dt;
            a.Kvc = a.f * (dt / m.getmass(i));
            a.x += a.Kxc - a.Kxb / 2;
            a.v += a.Kvc - a.Kvb / 2;
            //~ if(i == 0) cout << " " << a.f;
        }
    }
    // x = x0 + Kxc
    // v = v0 + Kvc
    //~ cout << (xold + (a0.Kxc) - a0.x) << "    "
         //~ << (vold + (a0.Kvc) - a0.v) << '\n';
    
    setForces(false);
    update_constraints();
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            atomRK4 &a = *((atomRK4 *) m.get(i));
            a.a = a.f / m.getmass(i);
            a.Kxd = a.v * dt;
            a.Kvd = a.a * dt;
            a.x += (a.Kxa + a.Kxd)/6 + (a.Kxb)/3 - (a.Kxc * (2.0/3.0));
            a.v += (a.Kva + a.Kvd)/6 + (a.Kvb)/3 - (a.Kvc * (2.0/3.0));
            //~ if(i == 0){
                
                //~ cout << " " << a.f << '\n';
                //~ cout << "x " << a.x << a.Kxa << a.Kxb << a.Kxc << a.Kxd << '\n' 
                     //~ << "v " << a.v << a.Kva << a.Kvb << a.Kvc << a.Kvd << '\n';
                 //~ }
        }
    }
    //~ Vec x1 = xold + (a0.Kxa/6) + (a0.Kxb/3) + (a0.Kxc/3) + (a0.Kxd/6);
    //~ Vec v1 = xold + (a0.Kxa/6) + (a0.Kxb/3) + (a0.Kxc/3) + (a0.Kxd/6);
    //~ cout << x1 << a0.x << (a0.x - x1) << "\n"
         //~ << v1 << a0.v << (a0.v - v1) << '\n';
    
    setForces(false);
    update_constraints();
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
        }
    }
    
    update_trackers();
    
    return;
}

flt collectionGear4NPH::setForcesGetPressure(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        group.resetForces();
    }
    
    flt p=0;
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        p += inter.setForcesGetPressure(box);
    }
    return p;
}

flt collectionGear4NPH::kinetic(){
    flt E=0;
    flt Vfac = dV/box->V()/3.0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            Vec v = g[i].v - (g[i].x * Vfac);
            E += v.sq() * g.getmass(i);
        }
    }
    return E/2.0;
}

flt collectionGear4NPH::temp(){
    flt totatoms = 0;
    flt totkinetic2 = 0; // kinetic * 2
    Vec cv = comv();
    flt Vfac = dV/box->V()/3.0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;for(uint i=0; i<g.size(); i++){
            Vec v = g[i].v - cv - (g[i].x * Vfac);
            totkinetic2 += v.sq() * g.getmass(i);
        }
        totatoms += g.size();
    }
    
    int ndof = dof();
    return totkinetic2 / (ndof-3);
}

void collectionGear4NPH::timestep(){
    uint n=0;
    OriginBox* obox = (OriginBox*) box;
    
    /// Prediction
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            g[i].x += g[i].v * dt + g[i].a * (dt*dt/2) + bs[n] * (dt*dt*dt/6);
            g[i].v += g[i].a * dt + bs[n] * (dt*dt/2);
            g[i].a += bs[n] * dt;
            n++;
        }
    }
    flt V = obox->V();
    flt newV = V + dV * dt + ddV * (dt*dt/2) + dddV * (dt*dt*dt/6);
    dV += ddV * dt + dddV * (dt*dt/2);
    ddV += dddV * dt;
    V = obox->resizeV(newV);
    //~ if (abs(V-newV)/abs(newV) > 1e-8){
        //~ printf("V: %.4g; newV: %.4g\n", V, newV);
        //~ assert(abs(V-newV)/abs(newV) < 1e-8);
    //~ }
    
    // Now we set forces and accelerations
    for(uint m=0; m<ncorrec; m++){
        flt interacP = setForcesGetPressure();
        flt newP = (interacP + (kinetic()*2.0))/3.0/V;
        //~ for(git = groups.begin(); git<groups.end(); git++){
            //~ atomgroup &g = **git;
            //~ for(uint i=0; i<g.size(); i++){
                //~ Vec a = g[i].f / g.getmass(i);
            //~ }
        //~ }
        //~ newP = pressure();
        update_constraints();
        
        /// Correction
        flt Vfac = (-(2.0/9.0)*(dV/V)*(dV/V)) + (ddV/V/3.0);
        flt newddV = (newP - P) / Q;
        flt correctionV = newddV - ddV;
        newV = V + correctionV * (dt*dt/12);
        dV += correctionV * (5*dt/12);
        ddV = newddV;
        dddV += correctionV / dt;
        
        V = obox->resizeV(newV);
        //~ if (abs(V-newV)/abs(newV) > 1e-8){
            //~ printf("V: %.4g; newV: %.4g\n", V, newV);
            //~ assert(abs(V-newV)/abs(newV) < 1e-8);
        //~ }
        
        
        n=0;
        for(git = groups.begin(); git<groups.end(); git++){
            atomgroup &g = **git;
            for(uint i=0; i<g.size(); i++){
                Vec a = g[i].f / g.getmass(i) + g[i].x * Vfac;
                Vec correction = a - g[i].a;
                g[i].x += correction * (dt*dt/12);
                g[i].v += correction * (5*dt/12);
                g[i].a = a;
                bs[n] += correction / dt;
                n++;
            }
        }
    }
    update_trackers();
    return;
}

vector<interaction*> collectionGear4NPT::tointerpair(vector<interactionpairsx*> &v){
    vector<interaction*> newv = vector<interaction*>();
    vector<interactionpairsx*>::iterator it;
    for(it = v.begin(); it<v.end(); it++){
        newv.push_back((interaction*) *it);
    }
    return newv;
}

void xrpsummer::run(forcepairx *fpairx){
    Vec rij = box->diff(fpairx->a1->x, fpairx->a2->x);
    Vec vij = fpairx->a1->v - fpairx->a2->v;
    rpxsum += rij.dot(vij) * fpairx->xij / rij.sq();
    xsum += fpairx->xij;
    vfsum += vij.dot(fpairx->fij);
    rfsum += rij.dot(fpairx->fij);
}

void collectionGear4NPT::setForces(bool seta){
    xrpsums.reset();
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        group.resetForces();
    }
    
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interactionpairsx* inter = (interactionpairsx*) *it;
        inter->setForces(box, &xrpsums);
    }
    if(!seta) return;
    
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].a = m[i].f / m.getmass(i);
        }
    }
}

void collectionGear4NPT::timestep(){
    uint n=0;
    OriginBox* obox = (OriginBox*) box;
    
    const flt c1 = dt, c2 = dt*dt/2, c3 = dt*dt*dt/6;
    
    /// Prediction
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            g[i].x += xs1[n]*c1 + xs2[n]*c2 + xs3[n]*c3;
            xs1[n] += xs2[n]*c1 + xs3[n]*c2;
            xs2[n] += xs3[n]*c1;
            
            g[i].v += g[i].a*c1 + vs2[n]*c2 + vs3[n]*c3;
            g[i].a += vs2[n]*c1 + vs3[n]*c2;
            vs2[n] += vs3[n]*c1;
            n++;
        }
    }
    
    flt V = obox->V();
    flt newV = V + V1*c1 + V2*c2 + V3*c3;
    V1 += V2*c1 + V3*c2;
    V2 += V3*c1;
    V = obox->resizeV(newV);
    //~ if (abs(V-newV)/abs(newV) > 1e-8){
        //~ printf("V: %.4g; newV: %.4g\n", V, newV);
        //~ assert(abs(V-newV)/abs(newV) < 1e-8);
    //~ }
    
    // Now we set forces and accelerations
    
    const flt d0 = 3.0*dt/8.0, d2 = 3.0/2.0/dt, d3 = 1.0/dt/dt;
    
    for(uint m=0; m<ncorrec; m++){
        setForces(false);
        update_constraints();
        
        /// Correction
        flt K2 = 2.0 * kinetic();
        flt PV = (xrpsums.rfsum + K2)/3.0;
        flt rpx = xrpsums.rpxsum;
        flt xx = xrpsums.xsum;
        chi = (-rpx) / (9*PV + xx);
        chixi = (xrpsums.vfsum) / K2;
        
        flt Vcorr = 3*V*chi - V1;
        V = obox->resizeV(V + Vcorr*d0);
        V1 += Vcorr;
        V2 += Vcorr*d2;
        V3 += Vcorr*d3;
        
        
        n=0;
        for(git = groups.begin(); git<groups.end(); git++){
            atomgroup &g = **git;
            for(uint i=0; i<g.size(); i++){
                Vec newx1 = g[i].v + (g[i].x * chi);
                Vec xcorr = newx1 - xs1[n];
                g[i].x += xcorr*d0;
                xs1[n] = newx1;
                xs2[n] += xcorr*d2;
                xs3[n] += xcorr*d3;
                
                Vec newa = g[i].f / g.getmass(i) - g[i].v*chixi;
                Vec vcorr = newa - g[i].a;
                g[i].v += vcorr*d0;
                g[i].a = newa;
                vs2[n] += vcorr*d2;
                vs3[n] += vcorr*d3;
                n++;
            }
        }
    }
    update_trackers();
    return;
}

void collectionVerletNPT::resetvhalf(){
    uint Natoms = 0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        Natoms += (*git)->size();
    };
    vhalf.resize(Natoms, Vec());
    uint n=0;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i=0; i<g.size(); i++){
            vhalf[n] = g[i].v;
            n++;
        }
    }
}

void collectionVerletNPT::timestep(){
    OriginBox* obox = (OriginBox*) box;
    setForces(false);
    update_constraints();
    
    //~ vector<Vec> lastvhalf = vhalf; // This does copy, right?
    uint ndof = dof();
    flt xieta = (xidot + eta) * dt/2.0;
    flt Ktot2 = 0, lastKtot2 = 0;
    
    uint n=0;
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            flt mi = m.getmass(i);
            m[i].a = m[i].f / mi;
            lastKtot2 += vhalf[n].sq() / mi;
            Vec lastvhalf = vhalf[n];
            vhalf[n] = (vhalf[n]*(1-xieta) + m[i].a*dt)/(1+xieta);
            //~ if(n == 0){
                //~ cout << "v: " << lastvhalf << " -> " << vhalf[n]
                     //~ << ", xieta: " << xieta << ", a: " << m[i].a
                     //~ << "\n";
                //~ }
            //~ 
            //~ if(vhalf[n].mag() > 0.6){
                //~ cout << "n: " << n
                     //~ << ", v: " << lastvhalf << " -> " << vhalf[n]
                     //~ << ", xieta: " << xieta << ", a: " << m[i].a
                     //~ << "\n";
                //~ }
            m[i].v = (lastvhalf + vhalf[n])/2.0;
            Ktot2 += vhalf[n].sq() / mi;
            n++;
        }
    }
    
    eta += dt*(Ktot2 - ndof*T)/QT;
    
    flt V = box->V();
    flt newV = obox->resizeV(lastV + (dt*xidot*V));
    flt Verr = newV / (lastV + (dt*xidot*V));
    
    lastV = V;
    curP = ((Ktot2 + lastKtot2)/2 + virial()) / 3.0 / V;
    flt xidott = xidot;
    xidot = lastxidot + (2*dt*(curP - P)*V)/QP;
    lastxidot = xidott;
    
    if((not (Verr < 1.001)) or (not (Verr > 0.999))){
        cout << "lastV: " << V << ", Set to: " <<  newV / Verr 
             << ", got: " << newV << ", Verr: " << Verr << '\n';
    }
    
    assert(Verr < 1.001);
    assert(Verr > 0.999);
    
    n=0;
    flt Vfac1=pow(newV/V, 1.0/3.0), Vfac2 = pow(2*newV/(newV+V), 1.0/3.0);
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            m[i].x = m[i].x*Vfac1 + vhalf[n]*dt*Vfac2;
            n++;
        }
    }
    
    update_trackers();
    
    return;
}
