#include "collection.hpp"

collection::collection(vector<atomgroup*> gs, vector<interaction*> is,
              vector<statetracker*> ts, vector<constraint*> cs)
        : groups(gs), interactions(is), trackers(ts), constraints(cs),
                atoms(gs){
    update_trackers();
    setForces();
    update_constraints();
    update_trackers();
}

void collection::scaleVelocities(flt T){
    flt t = temp();
    flt scaleby = sqrt(T/t);
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.size(); i++){
            g[i].v *= scaleby;
        }
    }
}

void collection::update_trackers(){
    vector<statetracker*>::iterator git;
    for(git = trackers.begin(); git<trackers.end(); git++){
        (*git)->update();
    }
}

void collection::update_constraints(){
    vector<constraint*>::iterator git;
    for(git = constraints.begin(); git<constraints.end(); git++){
        (*git)->apply();
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

flt collection::potentialenergy(){
    flt E=0;
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        E += inter.energy();
    }
    assert(not isnan(E));
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
    //~ Vec totcom = Vec(0,0,0);
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
    //~ Vec totcom = Vec(0,0,0);
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
    Vec avgr = Vec(0,0,0);
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

void collection::setForces(){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        group.resetForces();
    }
    
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        inter.setForces();
    }
}

collectionSol::collectionSol(const flt dt, const flt damp,const flt T, 
        vector<atomgroup*> groups,vector<interaction*> interactions,
        vector<statetracker*> trackers, vector<constraint*> constraints) :
            collection(groups, interactions, trackers, constraints), 
            dt(dt), damping(damp), desT(T){
    setCs();
    update_trackers();
    setForces();
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
            else vecpair[0] = vecpair[1] = Vec(0,0,0);
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
    setForces();
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
    
    return;
    //Honeycutt and Thirumalai
    /*
    cout << "damping " << damping << "\n";
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            // step 1: get new x
            m[i].x += (m[i].v * dt) + (m[i].a * (.5 * dt*dt));
            
            //(useful values)
            flt mass = m.getmass(i);
            flt curdamp = damping / mass;
            flt curdamp2 = dt * curdamp / 2;
            
            // step 2: make v1 (intermediate v)
            m[i].v += m[i].v * (-curdamp2 + curdamp2*curdamp2)
                        + m[i].f * (dt / 2 / mass);
        }
    }
    
    // step 3: set forces
    setForces();
    
    // step 4: intermediate acceleration
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.size(); i++){
            Vec g = gauss.genVec();
            // cout << "g " << g << "\n";
            m[i].a = (m[i].f + g) / m.getmass(i);
            
            // step 5: finish v and a
            m[i].v += m[i].a * (dt/2);
            m[i].a -= m[i].v * damping;
        }
    }
    update_trackers();
    */
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
    setForces();
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
    setForces();
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

void collectionGaussianT::setForces(bool shouldwesetxi){
    collection::setForces();
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
