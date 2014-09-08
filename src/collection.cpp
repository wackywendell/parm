#include "collection.hpp"

collection::collection(sptr<Box> box, sptr<atomgroup> atoms, vector<sptr<interaction> > is,
              vector<sptr<statetracker> > ts, vector<sptr<constraint> > cs)
        : box(box), atoms(atoms), interactions(is), trackers(ts), constraints(cs),
                E0(0){
    update_trackers();
    setForces(true);
    update_constraints();
    update_trackers();
}

void collection::scaleVs(flt scaleby){
    for(uint i = 0; i<atoms->size(); i++){
        (*atoms)[i].v *= scaleby;
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
    vector<sptr<statetracker> >::iterator git;
    for(git = trackers.begin(); git<trackers.end(); ++git){
        (*git)->update(*box);
    }
}

void collection::update_constraints(){
    vector<sptr<constraint> >::iterator git;
    for(git = constraints.begin(); git<constraints.end(); ++git){
        (*git)->apply(*box);
    }
}


flt collection::virial(){
    flt E = 0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        E += (*it)->pressure(*box);
    }
    return E;
}

flt collection::pressure(){
    // Returns (1/d V) Σ ri dot fi 
    // where d is number of dimensions
    // note that interaction->pressure just returns Σ ri dot fi
    flt V = box->V();

    flt E = 2.0 * kinetic();// * (ndof - nc - 3) / ndof;
    //flt E = (ndof - nc) * temp();
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        E += (*it)->pressure(*box);
        assert(not isnan(E));
    }
    return E / V / flt(NDIM);
}

flt collection::potentialenergy(){
    flt E=0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        interaction &inter = **it;
        E += inter.energy(*box);
        //~ cout << "potential energy: " << E << endl;
        assert(not isnan(E));
    }
    return E;
}

flt collection::energy(){
    flt E = potentialenergy() - E0 + kinetic();
    assert(not isnan(E));
    return E;
};

flt collection::dof(){
    int ndof = 0;
    vector<sptr<constraint> >::iterator cit;
    for(cit = constraints.begin(); cit<constraints.end(); ++cit){
        ndof -= (*cit)->ndof();
    }
    
    return ndof + NDIM*(atoms->size());
}

flt collection::temp(bool minuscomv){
    Vec v = Vec();
    if(minuscomv) v = comv();
    
    int ndof = (int) dof();
    if (minuscomv) ndof -= NDIM;
    return atoms->kinetic(v) * 2 / ndof;
}

flt collection::gyradius(){
    Vec avgr = Vec();
    flt N = atoms->size();
    for(uint i = 0; i<N; i++){
        avgr += (*atoms)[i].x;
    }
    avgr /= N; // now avgr is the average location, akin to c.o.m.
    flt Rgsq = 0;
    for(uint i = 0; i<atoms->size(); i++){
        Rgsq += ((*atoms)[i].x - avgr).sq();
    }
    
    return sqrt(Rgsq/N);
}

void collection::setForces(bool seta){
    atoms->resetForces();
    
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        interaction &inter = **it;
        inter.setForces(*box);
    }
    if(!seta) return;
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].a = (*atoms)[i].f / (*atoms)[i].m;
    }
}

flt collection::setForcesGetPressure(bool seta){
    atoms->resetForces();
    
    flt p=0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        interaction &inter = **it;
        p += inter.setForcesGetPressure(*box);
    }
    assert(!isnan(p));
    if(!seta) return p;
    
    for(uint i=0; i<atoms->size(); i++){
            (*atoms)[i].a = (*atoms)[i].f / (*atoms)[i].m;
    }
    assert(!isnan(p));
    return p;
}

collectionSol::collectionSol(sptr<Box> box, sptr<atomgroup> atoms,
        const flt dt, const flt damp, const flt T,
        vector<sptr<interaction> > interactions,
        vector<sptr<statetracker> > trackers, vector<sptr<constraint> > constraints) :
            collection(box, atoms, interactions, trackers, constraints), 
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
    // From Allen and Tildesley 263, Verlet-like, our equations are
    // r(t+dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + drG
    // v(t+dt) = c0 v(t) + (c1 - c2) dt a(t) + c2 dt a(t+dt) + dvG
    // I split this:
    // r(t+dt) = r(t) + c1 dt v(t) + c2 dt^2 a(t) + drG
    // vp(t) = c0 v(t) + (c1 - c2) dt a(t) + dvG
    // (set forces)
    // v(t+dt) = vp(t) + c2 dt a(t+dt)
    
    //Step 1: set atom.x = r(t+dt) and atom.v = vp(t)
        atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        flt v0 = sqrt(desT/m[i].m);
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
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(uint i=0; i<m.size(); i++){
        m[i].a = m[i].f / m[i].m;
    }
    
    // And finish m[i].v
    for(uint i=0; i<m.size(); i++){
    m[i].v += m[i].a * (dt*c2);
    }
    
    update_trackers();
};

collectionSolHT::collectionSolHT(sptr<Box> box, sptr<atomgroup> atoms,
        const flt dt, const flt damp,const flt T,
        vector<sptr<interaction> > interactions,
        vector<sptr<statetracker> > trackers, vector<sptr<constraint> > constraints) :
            collection(box, atoms, interactions, trackers, constraints), 
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
    
    flt keepv = 1-(damping*dt);
    flt xpartfromv = dt - (dt*dt*damping/2);
    
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        // step 1: get new x
        m[i].x += (m[i].v * xpartfromv) + (m[i].a * (.5 * dt*dt));
        
        // step 2: make v1 (intermediate v)
        m[i].v = m[i].v*keepv + m[i].a*(dt/2);
    }
    
    // step 3: set forces
    setForces();
    
    // step 4: intermediate acceleration
    for(uint i=0; i<m.size(); i++){
        Vec g = gauss.generate();
        //~ cout << "g " << g << "\n";
        m[i].a = (m[i].f + g) / m[i].m;
        
        // step 5: finish v
        m[i].v += m[i].a * (dt/2);
    }
    update_trackers();
}

void collectionVerlet::timestep(){
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].x += m[i].v * dt + m[i].a * (dt*dt/2);
        m[i].v += m[i].a * (dt/2);
    }
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(uint i=0; i<m.size(); i++){
        m[i].a = m[i].f / m[i].m;
        // And finish m[i].v
        m[i].v += m[i].a * (dt/2);
    }
    
    update_trackers();
}

void collectionOverdamped::timestep(){
    setForces(false);
    update_constraints();
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].a = m[i].f / m[i].m;
        m[i].v = m[i].a * gamma;
        m[i].x += m[i].v * dt;
    }
    
    update_trackers();
}

void collectionConjGradient::timestep(){
    // Algorithm from https://en.wikipedia.org/wiki/Energy_minimization#Nonlinear_conjugate_gradient_method
    // Where dt = κ, v = h, F = F
    // Note that in the middle of this loop, a = a(t-1), and newa = a(t)
    // and while its called 'a', its actually 'v'
    setForces(false);
    update_constraints();
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].x += m[i].v * dt;
        Vec newa = m[i].f / m[i].m;
        flt asq = m[i].a.sq();
        flt gamma = 0;
        if(asq > 0){// if a was 0, then just set γ = 0
            gamma = (newa - m[i].a).dot(newa) / asq; // Polak-Ribière
            // gamma = newa.sq() / asq; // Fletcher-Reeves
        }
        //if(gamma > 100) gamma = 0;
        if(gamma < 0) gamma = 0;
        m[i].v = newa + (m[i].v * gamma);
        m[i].a = newa;
        //~ if(git == groups.begin() and i == 0){
            //~ cout << ' ' << i << " :: Gamma: " << gamma << " asq:" << asq << " newasq:" << newa.sq() << endl;
            //~ //cout << "        " << 
        //~ }
        
    }
    
    update_trackers();
}

void collectionConjGradient::timestepNewton(){
    setForces(false);
    update_constraints();
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = (*atoms)[i].f / (*atoms)[i].m;
        (*atoms)[i].x += (*atoms)[i].v * dt;
    }
    
    update_trackers();
}

void collectionConjGradient::reset(){
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = Vec();
    }
}

flt collectionConjGradientBox::kinetic(){
    flt E=0;
    flt Vfac = dV/box->V()/flt(NDIM);
    for(uint i=0; i<atoms->size(); i++){
        //Vec v = g[i].v - (g[i].x * Vfac);
        Vec v = (*atoms)[i].v + ((*atoms)[i].x * Vfac);
        E += v.sq() * (*atoms)[i].m;
    }
    return E/2.0;
}

void collectionConjGradientBox::reset(){
    dV = 0;
    hV = 0;
    FV = 0;
    lastFV = 0;
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = Vec();
    }
}


void collectionConjGradientBox::resize(flt V){
    dV = 0;
    hV = 0;
    FV = 0;
    lastFV = 0;
    
    OriginBox& obox = (OriginBox&) *box;
    flt oldV = obox.V();
    #ifdef VEC3D
    flt Vfac = pow(V/oldV, 1.0/3.0);
    #endif
    #ifdef VEC2D
    flt Vfac = sqrt(V/oldV);
    #endif
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x *= Vfac;
        (*atoms)[i].v = Vec();
    }
    
    obox.resizeV(V);
}

void collectionConjGradientBox::timestep(){
    // Algorithm from https://en.wikipedia.org/wiki/Energy_minimization#Nonlinear_conjugate_gradient_method
    // Where dt = κ, v = h, F = F
    // Note that in the middle of this loop, a = a(t-1), and newa = a(t)
    // and while its called 'a', its actually 'v'
    
    OriginBox& obox = (OriginBox&) *box;
    atomgroup& g = (*atoms);
    update_constraints();
    
    flt oldV = obox.V();
    dV = hV * oldV;
    if(dV*dt > 1/kappaV) dV = 1/kappaV/dt;
    if(dV < 0) dV *= kappaV;
    if(maxdV >= 0 and dV > maxdV*oldV) dV = maxdV*oldV;
    if(maxdV >= 0 and dV < -maxdV*oldV) dV = -maxdV*oldV;
    if(maxdV == 0) dV = 0;

    flt newV = oldV + (dV*dt);
    //flt Vfac = dV/((oldV + newV)/2)/NDIM;
    #ifdef VEC3D
    flt Vfac = pow(newV/oldV, 1.0/3.0);
    #endif
    #ifdef VEC2D
    flt Vfac = sqrt(newV/oldV);
    #endif
    
    for(uint i=0; i<g.size(); i++){
        g[i].x *= Vfac;
        g[i].x += g[i].v * dt;// + (m[i].x * Vfac);
    }
    
    obox.resizeV(newV);
    flt interacP = setForcesGetPressure();
    flt newFV = (interacP/NDIM) - (P0*newV);
    
    flt gammaV = 0;
    if(FV*FV > 0) gammaV = newFV * (newFV-FV) / (FV*FV);
    if(gammaV < 0 or isnan(gammaV) or isinf(gammaV)) gammaV = 0;
    if(gammaV > 1) gammaV = 1;
    
    if(isnan(oldV) or isnan(newV) or isnan(dV)){
        cout << "P: " << (interacP/NDIM) << " - " << P0 << " V: " << dV << " / (" << oldV << " -> " << newV << ")\n";
        cout << "FV: " << FV << " -> " << newFV << " gammaV: " << gammaV << endl;
        assert(!isnan(oldV));
    }
    hV = newFV + (gammaV*hV);
    FV = newFV;
    
    for(uint i=0; i<g.size(); i++){
        Vec newa = g[i].f / g[i].m;
        flt asq = g[i].a.sq();
        flt gamma = 0;
        if(asq > 0){// if a was 0, then just set γ = 0
            gamma = (newa - g[i].a).dot(newa) / asq; // Polak-Ribière
            // gamma = newa.sq() / asq; // Fletcher-Reeves
        }
        //if(gamma > 100) gamma = 0;
        if(gamma < 0 or isnan(gamma) or isinf(gamma)) gamma = 0;
        g[i].v = newa + (g[i].v * gamma);
        g[i].a = newa;
        
    }
    
    update_trackers();
}

void collectionConjGradientBox::timestepBox(){
    OriginBox& obox = (OriginBox&) *box;
    update_constraints();
    
    flt oldV = obox.V();
    dV = hV * oldV;
    if(dV*dt > 1/kappaV) dV = 1/kappaV/dt;
    if(dV < 0) dV *= kappaV;
    if(maxdV >= 0 and dV > maxdV*oldV) dV = maxdV*oldV;
    if(maxdV >= 0 and dV < -maxdV*oldV) dV = -maxdV*oldV;
    flt newV = oldV + (dV*dt);
    
    #ifdef VEC3D
    flt Vfac = pow(newV/oldV, 1.0/3.0);
    #endif
    #ifdef VEC2D
    flt Vfac = sqrt(newV/oldV);
    #endif
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x *= Vfac;
        (*atoms)[i].v = Vec();
    }
    
    obox.resizeV(newV);
    flt interacP = setForcesGetPressure();
    flt newFV = (interacP/NDIM) - (P0*newV);
    
    flt gammaV = 0;
    if(FV*FV > 0) gammaV = newFV * (newFV) / (FV*FV);
    if(gammaV < 0 or isnan(gammaV) or isinf(gammaV)) gammaV = 0;
    if(gammaV > 1) gammaV = 1;
    
    if(isnan(oldV) or isnan(newV) or isnan(dV)){
        cout << "P: " << (interacP/NDIM) << " - " << P0 << " V: " << dV << " / (" << oldV << " -> " << newV << ")\n";
        cout << "FV: " << FV << " -> " << newFV << " gammaV: " << gammaV << endl;
        assert(!isnan(oldV));
    }
    hV = newFV + (gammaV*hV);
    FV = newFV;
    
    update_trackers();
    return;
}

void collectionConjGradientBox::timestepAtoms(){
    // Algorithm from https://en.wikipedia.org/wiki/Energy_minimization#Nonlinear_conjugate_gradient_method
    // Where dt = κ, v = h, F = F
    // Note that in the middle of this loop, a = a(t-1), and newa = a(t)
    // and while its called 'a', its actually 'v'
    
    setForces(false);
    update_constraints();
    
    dV = 0;
    hV = 0;
    FV = 0;
    lastFV = 0;
    
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].x += m[i].v * dt;// + (m[i].x * Vfac);
        
        Vec newa = m[i].f / m[i].m;
        flt asq = m[i].a.sq();
        flt gamma = 0;
        if(asq > 0){// if a was 0, then just set γ = 0
            gamma = (newa - m[i].a).dot(newa) / asq; // Polak-Ribière
            // gamma = newa.sq() / asq; // Fletcher-Reeves
        }
        //if(gamma > 100) gamma = 0;
        if(gamma < 0 or isnan(gamma) or isinf(gamma)) gamma = 0;
        m[i].v = newa + (m[i].v * gamma);
        m[i].a = newa;            
    }
    
    update_trackers();
}

collectionNLCG::collectionNLCG(sptr<OriginBox> box, sptr<atomgroup> atoms,
                const flt dt, const flt P0,
                vector<sptr<interaction> > interactions,
                vector<sptr<statetracker> > trackers,
                vector<sptr<constraint> > constraints,
                const flt kappa, const flt kmax,
                const uint secmax, const flt seceps) :
            collection(box, atoms, interactions, trackers, 
                constraints), dt(dt), seceps(seceps), secmax(secmax),
                kappa(kappa), alphamax(2.0), afrac(0), dxmax(100),
                stepmax(1e-3), kmax(kmax), 
                P0(P0), Knew(0), k(0), vl(0), fl(0), al(0), alpha(0), 
                beta(0), betaused(0), dxsum(0), alphavmax(0), maxdV(0),
                sec(0){
};

void collectionNLCG::reset(){
    k = 0;
    setForces(true, true);
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = (*atoms)[i].a;
    }
    vl = al;
}

void collectionNLCG::setForces(bool seta, bool setV){
    flt V = box->V();
    if(setV){
        flt interacP = collection::setForcesGetPressure(false);
        fl = ((interacP/NDIM) - (P0*V))/kappa;
        assert(kappa > 0);
        assert(!isnan(interacP));
        assert(!isnan(V));
        assert(!isnan(P0));
        assert(!isnan(fl));
        if(seta){
            al = fl;
            vl = fl;
        }
        //~ cout << "NLCG::setForces: fl = " << fl
             //~ << " -- P/d = " << (interacP/NDIM) << " -- V*P0 = " << (P0*V)
             //~ << '\n';
    } else {
        collection::setForces(false);
    }
    if(seta){
        for(uint i=0; i<atoms->size(); i++){
            (*atoms)[i].a = (*atoms)[i].f;
            (*atoms)[i].v = (*atoms)[i].f;
        }
        Knew = fdota();
    }
}


flt collectionNLCG::Hamiltonian(){
    return potentialenergy() + P0*(box->V());
};

void collectionNLCG::resize(flt V){
    reset();
    OriginBox& obox = (OriginBox&) *box;
    flt oldV = obox.V();
    #ifdef VEC3D
    flt Lfac = pow(V/oldV, 1.0/3.0);
    #endif
    #ifdef VEC2D
    flt Lfac = sqrt(V/oldV);
    #endif
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x *= Lfac;
        (*atoms)[i].v = Vec();
    }
    
    obox.resizeV(V);
}

flt collectionNLCG::kinetic(){
    flt E=0;
    flt Lfac = exp(vl / (kappa*NDIM));
    for(uint i=0; i<atoms->size(); i++){
        //Vec v = g[i].v - (g[i].x * Vfac);
        Vec v = (*atoms)[i].v + ((*atoms)[i].x * Lfac);
        E += v.sq();// * g[i].m;
    }
    
    return E/2.0;
}

flt collectionNLCG::pressure(){
    flt V = box->V();
    
    flt E = 0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        E += (*it)->pressure(*box);
        assert(not isnan(E));
    }
    return E / V / flt(NDIM);
}

void collectionNLCG::stepx(flt dx){
    OriginBox& obox = (OriginBox&) *box;
    flt Lfac = exp(dx * vl / (kappa*NDIM));
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x *= Lfac;
        (*atoms)[i].x += (*atoms)[i].v * dx;
    }
    
    //~ cout << "NLCG::stepx: V0 = " << obox->V() << " -- dx = " << dx 
         //~ << " -- Vfac = " << Vfac;
    obox.resize(Lfac);
    //~ cout << " -- V = " << obox->V() << '\n';
}

flt collectionNLCG::getLsq(){
    OriginBox& obox = (OriginBox&) *box;
    #ifdef VEC3D
    return pow(obox.V(), 2.0/3.0);
    #endif
    #ifdef VEC2D
    return obox.V();
    #endif
};

flt collectionNLCG::fdotf(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.sq();
    }
    
    return returnvalue / getLsq() + fl * fl;
}

flt collectionNLCG::fdota(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.dot((*atoms)[i].a);
    }
    
    return returnvalue / getLsq() + fl * al;
}

flt collectionNLCG::fdotv(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.dot((*atoms)[i].v);
    }
    
    return returnvalue / getLsq() + fl * vl;
}

flt collectionNLCG::vdotv(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].v.sq();
    }
    
    return returnvalue / getLsq() + vl * vl;
}

void collectionNLCG::timestep(){
    // Algorithm from
    // An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
    // Edition 1 1/4
    // Jonathan Richard Shewchuk
    // August 4, 1994
    // An electronic copy of this article is available by anonymous FTP to WARP.CS.CMU.EDU (IP address 128.2.209.103) under the filename quake-papers/painless-conjugate-gradient.ps. 
    
    // Uses a combination of
    // B4. Nonlinear Conjugate Gradients with Newton-Raphson and Fletcher-Reeves
    // B5. Preconditioned Nonlinear Conjugate Gradients with Secant and Polak-Ribi`ere
    
    // where I use Secant and Polak-Ribi`ere, but no preconditioning.
    
    // Where m[i].a/L = -f', m[i].v/L = d, σ0 = dt
    // Note that in the middle of this loop, a = a(t-1), and newa = a(t)
    // and while its called 'a', its actually 'v'
    
    //~ flt V0 = box->V();
    stepx(dt);
    setForces(false, true); // sets both atom.f and fl, but not atom.a or al
    update_constraints();
    flt eta0 = -fdotv(); // slope at x0 - dt
    flt eta = eta0;
    stepx(-dt);
    update_trackers();
    alpha = -dt;
    setForces(false, true); // sets both atom.f and fl, but not atom.a or al
    update_constraints();
    dxsum = 0;
    flt vdv = vdotv();
    
    /// Secant stepping
    // note that v does not change here; we simply step in the direction
    // of v until the force along v (eta, fdotv) gets very small...
    for(sec=0; sec < secmax; sec++){
        eta = -fdotv(); // slope at x
                            // eta0 is the slope at x - α
        
        // Go in the direction of -slope (-eta), approximating a quadratic.
        // abs added to ensure we go in the direction of the slope, and don't
        // maximize a quadratic
        flt alphafac = -eta / abs(eta0-eta);
        if(abs(eta0 - eta) <= 1e-12 * abs(eta)){
            alphafac = alphamax > 0 ? alphamax : 1.1;
            sec = sec > 0 ? sec*2 - 1 : 1;
        }
        
        if(alphamax > 0 and alphafac > alphamax) alphafac = alphamax;
        //if(alphamax > 0 and alphafac > 0 and alphafac < alphamin) alphafac = alphamin;
        if(alphamax > 0 and alphafac < -alphamax) alphafac = -alphamax;
        //if(alphamax > 0 and alphafac < 0 and alphafac > -alphamin) alphafac = -alphamin;
        
        // abs added here, same reason
        alpha = abs(alpha) * alphafac;
        
        //~ if(abs(alpha) <= 1e-32 and abs(eta) >= 1e-8){
            //~ cout << "NLCGV::timestep 2: secant " << sec
                 //~ << " eta = " << eta << " -- alpha = " << alpha;
            //~ cout << " -- alpha0 = " << oldalpha;
            //~ cout << " -- alphafac = " << alphafac;
            //~ cout << " -- E = " << energy();
            //~ cout << " -- K = " << kinetic();
            //~ cout << " -- fdv = " << fdotv() << endl;
        //~ }
        
        //~ if(sec >= secmax - 9){
            //~ cout << "sec: " << sec << " eta: " << eta << " eta0: " << eta0 
                //~ << " diff: " << (eta0 - eta) << " alpha: " << alpha
                //~ << " alphafac: "  << alphafac << " alphamax: " << alphamax
                //~ << " dxsum: " << (dxsum + alpha) << endl;
        //~ }
        
        flt newdxsum = abs(dxsum + alpha);
        if(dxmax > 0 and newdxsum > dxmax){
            //~ cout << "dxmaxed: " << newdxsum << "  step: " << sqrt(dxsum*dxsum*vdv) << '\n';
            k = 0;
            break;
        }
        flt dVoverV = expm1(abs(dxsum + alpha) * vl / (kappa*NDIM));
        if(maxdV > 0 and dVoverV > maxdV){
            k = 0;
            //~ flt V1 = box->V();
            //~ cout << "V overflow, " << V0 << " -> " << V1 << '\n';
            break;
        }
        dxsum += alpha;
        stepx(alpha);
        setForces(false, true);
        update_constraints();
        eta0 = eta;
        if(alpha*alpha*vdv < seceps*seceps) break;
        if((sec > 1)
                and (afrac > 0)
                and (abs(alpha) < abs(dxsum))
                and (abs(alpha) / abs(dxsum) < afrac)) break;
        if(stepmax > 0 and dxsum*dxsum*vdv > stepmax*stepmax){
            k = 0;
            //~ cout << "stepmaxed: " << sqrt(dxsum*dxsum*vdv)
                //~ << "  alpha: " << alpha
                //~ << "  dx: " << dxsum
                //~ << "  eta: " << eta
                //~ << "  beta: " << beta << " -> " << betaused << '\n';
            break;
        }
    }
    
    alphavmax = sqrt(alpha*alpha*vdv);
    
    flt Kold = Knew;
    flt Kmid = fdota();
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].a = (*atoms)[i].f;
    }
    al = fl;
    
    Knew = fdota();
    beta = (Knew - Kmid) / Kold;
    betaused = beta;
    k++;
    if(k >= kmax or isinf(betaused) or isnan(betaused) or betaused <= 0){
        k = 0;
        betaused = 0;
    } else if(betaused > 1){
        betaused = 1;
    }
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = (*atoms)[i].a + (*atoms)[i].v*betaused;
    }
    vl = al + betaused * vl;
}

void collectionNLCG::descend(){
    setForces(false, true);
    update_constraints();
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = (*atoms)[i].f;
        (*atoms)[i].a = (*atoms)[i].f;
    }
    al = fl;
    vl = fl;
    
    stepx(dt);
    update_trackers();
}

collectionNLCGV::collectionNLCGV(sptr<Box> box, sptr<atomgroup> atoms,
                const flt dt,
                vector<sptr<interaction> > interactions,
                vector<sptr<statetracker> > trackers,
                vector<sptr<constraint> > constraints,
                const flt kmax,
                const uint secmax, const flt seceps) :
            collection(box, atoms, interactions, trackers, 
                constraints), dt(dt), seceps(seceps), secmax(secmax),
                alphamax(0), afrac(0), dxmax(0), stepmax(0),
                kmax(kmax), Knew(0), k(0), vl(0), fl(0), al(0), 
                alpha(0), beta(0), betaused(0), dxsum(0), alphavmax(0),
                sec(0){
};

void collectionNLCGV::stepx(flt dx){
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x += (*atoms)[i].v * dx;
    }
};

flt collectionNLCGV::fdotf(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.sq();
    }
    
    return returnvalue;
}

flt collectionNLCGV::fdota(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.dot((*atoms)[i].a);
    };
    
    return returnvalue;
}

flt collectionNLCGV::fdotv(){
    flt returnvalue = 0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].f.dot((*atoms)[i].v);
        //~ if(isnan(returnvalue)){
            //~ cout << "collectionNLCGV::fdotv : Got nan on atom " << i;
            //~ cout << " f: " << (*atoms)[i].f;
            //~ cout << " v: " << (*atoms)[i].v;
            //~ cout << '\n';
            //~ assert(false);
        //~ }
    }
    
    return returnvalue;
}

flt collectionNLCGV::vdotv(){
    flt returnvalue=0;
    for(uint i=0; i<atoms->size(); i++){
        returnvalue += (*atoms)[i].v.sq();
    }
    
    return returnvalue;
}

void collectionNLCGV::reset(){
    k = 0;
    setForces(false);
    for(uint i=0; i<atoms->size(); i++){
            (*atoms)[i].v = (*atoms)[i].f;
            (*atoms)[i].a = (*atoms)[i].f;
    }
}

void collectionNLCGV::descend(){
    setForces(false);
    update_constraints();
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v = (*atoms)[i].f;
        (*atoms)[i].a = (*atoms)[i].f;
    }

    stepx(dt);
    update_trackers();
}

flt collectionNLCGV::pressure(){
    flt V = box->V();
    
    flt E = 0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        E += (*it)->pressure(*box);
        assert(not isnan(E));
    }
    return E / V / flt(NDIM);
}



void collectionNLCGV::timestep(){
    // Algorithm from
    // An Introduction to the Conjugate Gradient Method Without the Agonizing Pain
    // Edition 1 1/4
    // Jonathan Richard Shewchuk
    // August 4, 1994
    // An electronic copy of this article is available by anonymous FTP to WARP.CS.CMU.EDU (IP address 128.2.209.103) under the filename quake-papers/painless-conjugate-gradient.ps. 
    
    // Uses a combination of
    // B4. Nonlinear Conjugate Gradients with Newton-Raphson and Fletcher-Reeves
    // B5. Preconditioned Nonlinear Conjugate Gradients with Secant and Polak-Ribi`ere
    
    // where I use Secant and Polak-Ribi`ere, but no preconditioning.
    
    // Where m[i].a = -f', m[i].v = d, σ0 = dt
    
    stepx(dt);
    setForces(false);
    update_constraints();
    flt eta0 = -fdotv();
    flt eta = eta0;
    //~ cout << "NLCGV::timestep 1:"
             //~ << " eta0 = " << eta0;
        //~ cout << " -- E = " << energy();
        //~ cout << " -- K = " << kinetic();
        //~ cout << " -- fdv = " << fdotv() << endl;
    
    stepx(-dt);
    update_trackers();
    alpha = -dt;
    setForces(false); // sets both atom.f and fl, but not atom.a or al
    update_constraints();
    dxsum = 0;
    flt vdv = vdotv();
    
    // flt alphamin = 1.0 / alphamax;
    
    /// Secant stepping
    // note that v does not change here; we simply step in the direction
    // of v until the force along v (eta, fdotv) gets very small...
    for(sec=0; sec < secmax; sec++){
        eta = -fdotv(); // slope at x
                            // eta0 is the slope at x - α
        
        // Go in the direction of -slope (-eta), approximating a quadratic.
        // abs added to ensure we go in the direction of the slope, and don't
        // maximize a quadratic
        flt alphafac = -eta / abs(eta0-eta);
        if(abs(eta0 - eta) <= 1e-12 * abs(eta)){
            alphafac = alphamax > 0 ? alphamax : 1.1;
            sec = sec > 0 ? sec*2 - 1 : 1;
        }
        
        if(alphamax > 0 and alphafac > alphamax) alphafac = alphamax;
        //if(alphamax > 0 and alphafac > 0 and alphafac < alphamin) alphafac = alphamin;
        if(alphamax > 0 and alphafac < -alphamax) alphafac = -alphamax;
        //if(alphamax > 0 and alphafac < 0 and alphafac > -alphamin) alphafac = -alphamin;
        
        // abs added here, same reason
        //~ flt oldalpha = alpha;
        alpha = abs(alpha) * alphafac;
        
        //~ if(abs(alpha) <= 1e-32 and abs(eta) >= 1e-8){
            //~ cout << "NLCGV::timestep 2: secant " << sec
                 //~ << " eta = " << eta << " -- alpha = " << alpha;
            //~ cout << " -- alpha0 = " << oldalpha;
            //~ cout << " -- alphafac = " << alphafac;
            //~ cout << " -- E = " << energy();
            //~ cout << " -- K = " << kinetic();
            //~ cout << " -- fdv = " << fdotv() << endl;
        //~ }
        
        //~ if(sec >= secmax - 9){
            //~ cout << "sec: " << sec << " eta: " << eta << " eta0: " << eta0 
                //~ << " diff: " << (eta0 - eta) << " alpha: " << alpha
                //~ << " alphafac: "  << alphafac << " alphamax: " << alphamax
                //~ << " dxsum: " << (dxsum + alpha) << endl;
        //~ }
        
        flt newdxsum = abs(dxsum + alpha);
        if(dxmax > 0 and newdxsum > dxmax){
            //~ cout << "dxmaxed: " << newdxsum << "  step: " << sqrt(dxsum*dxsum*vdv) << '\n';
            k = 0;
            break;
        }
        dxsum += alpha;
        stepx(alpha);
        
        //~ cout << "NLCGV::timestep 3: secant " << sec
             //~ << " eta = " << eta << " -- alpha = " << alpha;
        //~ cout << " -- E = " << energy();
        //~ cout << " -- K = " << kinetic();
        //~ cout << " -- fdv = " << fdotv() << endl;
        
        update_trackers();
        
        setForces(false);
        update_constraints();
        eta0 = eta;
        
        if(alpha*alpha*vdv < seceps*seceps) break;
        if((sec > 1)
                and (afrac > 0)
                and (abs(alpha) < abs(dxsum))
                and (abs(alpha) / abs(dxsum) < afrac)) break;
        if(stepmax > 0 and dxsum*dxsum*vdv > stepmax*stepmax){
            k = 0;
            //~ cout << "stepmaxed: " << sqrt(dxsum*dxsum*vdv)
                //~ << "  alpha: " << alpha
                //~ << "  dx: " << dxsum
                //~ << "  eta: " << eta
                //~ << "  beta: " << beta << " -> " << betaused << '\n';
            break;
        }
    }
    
    alphavmax = sqrt(alpha*alpha*vdv);
    
    flt Kold = Knew;
    flt Kmid = fdota();
    
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].a = m[i].f;
    }
    al = fl;
    
    Knew = fdota();
    beta = (Knew - Kmid) / Kold;
    betaused = beta;
    k++;
    if(k >= kmax or isinf(betaused) or isnan(betaused) or betaused <= 0){
        k = 0;
        betaused = 0;
    } else if(betaused > 1){
        betaused = 1;
    }
    
    for(uint i=0; i<m.size(); i++){
        m[i].v = m[i].a + m[i].v*betaused;
    }
    
    //~ cout << "NLCGV::timestep 4: beta " << beta;
        //~ cout << " -- E = " << energy();
        //~ cout << " -- K = " << kinetic();
        //~ cout << " -- fdv = " << fdotv() << endl;
        
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
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        Vec Ftilde = (m[i].a - (m[i].v*xi));
        m[i].x += m[i].v * dt + Ftilde * (dt*dt/2);
        // first half of y
        m[i].v += Ftilde * (dt/2);
    }
    
    lns += xi * dt + (Kt - ndof*T) * (dt*dt/2/Q);
    
    
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    for(uint i=0; i<m.size(); i++){
        m[i].a = m[i].f / m[i].m;
        // And finish y
        m[i].v += m[i].a * (dt/2);
    }
    
    // Now we solve for xi
    flt Ky = 2*kinetic();
    flt z0 = xi + (Kt - 2*ndof*T)*(dt/2/Q);
    flt z1 = Ky*2/dt/Q;
    
    //~ flt oldxi = xi;
    xi = solveCubic(4/dt - z0, 4/dt/dt - 4*z0/dt, -(z0*4/dt/dt) - z1, xi);
    
    
    // And finish v
    flt ytov = 1 + xi*dt/2;
    for(uint i=0; i<m.size(); i++){
        m[i].v /= ytov;
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
    for(uint i=0; i<atoms->size(); i++){
            flt m = (*atoms)[i].m;
            xiNumerator += (*atoms)[i].f.dot((*atoms)[i].v) * m;
            xiDenominator += ((*atoms)[i].v.sq()) * m*m;
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
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x += (*atoms)[i].v * dt + (*atoms)[i].a * (dt*dt/2);
        (*atoms)[i].v += ((*atoms)[i].a - ((*atoms)[i].v*xi)) * (dt/2);
    }
    // Now we set forces and accelerations
    setForces(false, false);
    update_constraints();
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].a = (*atoms)[i].f / (*atoms)[i].m;
        // And finish y
        (*atoms)[i].v += (*atoms)[i].a * (dt/2);
    }
    
    // now we get z and set xi
    //~ flt oldxi = xi;
    flt z = setxi();
    xi = z/(1 - z*dt/2);
    
    // And finish v
    flt ytov = 1 + xi*dt/2;
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].v /= ytov;
    }
    
    update_trackers();
}


void collectionGear3A::timestep(){
    atomgroup &m = *atoms;
    for(uint i=0; i<m.size(); i++){
        m[i].x += m[i].v * dt + m[i].a * (dt*dt/2);
        m[i].v += m[i].a * dt;
    }
    // Now we set forces and accelerations
    setForces(false);
    update_constraints();
    
    // Make correction
    for(uint i=0; i<m.size(); i++){
        Vec a = m[i].f / m[i].m;
        Vec correction = a - m[i].a;
        m[i].v += correction * (dt/2);
        m[i].a = a;
    }
    
    update_trackers();
}

void collectionGear4A::timestep(){
    uint n=0;
    //~ cout << "Gear 4A timestep, ncorrec = " << ncorrec << "\n";
    
    atomgroup &g = *atoms;
    //~ cout << "Gear 4A timestep loop 1, g.size() = " << g.size() << "\n";
    for(uint i=0; i<g.size(); i++){
        //~ cout << "Gear 4A timestep loop 2, i = " << i << "\n";
        g[i].x += g[i].v * dt + g[i].a * (dt*dt/2) + bs[n] * (dt*dt*dt/6);
        g[i].v += g[i].a * dt + bs[n] * (dt*dt/2);
        g[i].a += bs[n] * dt;
        n++;
    }
    
    //~ cout << "Gear 4A timestep 2, ncorrec = " << ncorrec << "\n";
    // Now we set forces and accelerations
    for(uint m=0; m<ncorrec; m++){
        //~ cout << "Gear 4A setting forces, m = " << m << ", ncorrec = " << ncorrec << "\n";
        setForces(false);
        update_constraints();
    
        // Make correction
        n=0;
        for(uint i=0; i<g.size(); i++){
            Vec a = g[i].f / g[i].m;
            Vec correction = a - g[i].a;
            g[i].x += correction * (dt*dt/12);
            g[i].v += correction * (5*dt/12);
            g[i].a = a;
            bs[n] += correction / dt;
            n++;
        }
    }
    update_trackers();
}

void collectionGear5A::timestep(){
    atomgroup &g = *atoms;
    uint n=0;
    flt dt2 = dt*dt/2;
    flt dt3 = dt*dt2/3;
    flt dt4 = dt*dt3/4;
    for(uint i=0; i<g.size(); i++){
        g[i].x += g[i].v * dt + g[i].a * dt2 + bs[n] * dt3 + cs[n] * dt4;
        g[i].v += g[i].a * dt + bs[n] * dt2 + cs[n] * dt3;
        g[i].a += bs[n] * dt + cs[n] * dt2;
        bs[n] += cs[n] * dt;
        n++;
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
        for(uint i=0; i<g.size(); i++){
            Vec a = g[i].f / g[i].m;
            Vec correction = a - g[i].a;
            g[i].x += correction * c0;
            g[i].v += correction * c1;
            g[i].a = a;
            bs[n] += correction * c3;
            cs[n] += correction * c4;
            n++;
        }
    }
    
    update_trackers();
}

void collectionGear6A::timestep(){
    atomgroup &g = *atoms;
    //~ cout << "Timestep\n";
    flt dt2 = dt*dt/2;
    flt dt3 = dt*dt2/3;
    flt dt4 = dt*dt3/4;
    flt dt5 = dt*dt4/5;
    
    uint n=0;
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
    
    flt c0 = 3*dt*dt/40, c1 = 251*dt/720;
    //flt c2 = 1
    flt c3 = 11 / (6*dt), c4 = 2/(dt*dt), c5 = 1/(dt*dt*dt);
    
    for(uint m=0; m<ncorrec; m++){
        // Now we set forces and accelerations
        setForces(false);
        update_constraints();
        
        //~ flt correctax = (gbegin[0].f / gbegin[0].m)[0];
        //~ printf("correcta = %.15f\n", correctax);
        //~ printf("correctascaled = %.15f\n", correctax*dt2);
    
        // Make correction
        n=0;
        for(uint i=0; i<g.size(); i++){
            //~ assert(n < bs.size());
            //~ assert(n < cs.size());
            //~ assert(n < ds.size());                
            Vec a = g[i].f / g[i].m;
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
}

void collectionRK4::timestep(){
    atomgroup &g = *atoms;
    //~ atom &atm0 = *((*(groups.begin()))->get(0));
    //~ atomRK4 &a0 = (atomRK4 &) atm0;
    //~ Vec xold = a0.x, vold = a0.v;
    
    
    for(uint i=0; i<g.size(); i++){
        atomid a = g.get_id(i);
        RK4data adat = data[a.n()];
        adat.Kxa = a->v * dt;
        adat.Kva = a->f * (dt / g[i].m);
        a->x += adat.Kxa / 2;
        a->v += adat.Kva / 2;
        //~ if(i == 0) cout << "f " << a.f;
    }
    
    setForces(false);
    update_constraints();
    
    for(uint i=0; i<g.size(); i++){
        atomid a = g.get_id(i);
        RK4data adat = data[a.n()];
        adat.Kxb = a->v * dt;
        adat.Kvb = a->f * (dt / g[i].m);
        a->x += adat.Kxb / 2 - adat.Kxa / 2;
        a->v += adat.Kvb / 2 - adat.Kva / 2;
        //~ if(i == 0) cout << " " << a.f;
    }
    
    //~ cout << (xold + (a0.Kxb / 2) - a0.x) << "    "
         //~ << (vold + (a0.Kvb / 2) - a0.v) << '\n';
    
    // x = x0 + Kxb/2
    // v = v0 + Kvb/2
    
    setForces(false);
    update_constraints();
    
    for(uint i=0; i<g.size(); i++){
        atomid a = g.get_id(i);
        RK4data adat = data[a.n()];
        adat.Kxc = a->v * dt;
        adat.Kvc = a->f * (dt / a->m);
        a->x += adat.Kxc - adat.Kxb / 2;
        a->v += adat.Kvc - adat.Kvb / 2;
    }
    
    setForces(false);
    update_constraints();
    
    
    for(uint i=0; i<g.size(); i++){
        atomid a = g.get_id(i);
        RK4data adat = data[a.n()];
        a->a = a->f / a->m;
        adat.Kxd = a->v * dt;
        adat.Kvd = a->a * dt;
        a->x += (adat.Kxa + adat.Kxd)/6 + (adat.Kxb)/3 - (adat.Kxc * (2.0/3.0));
        a->v += (adat.Kva + adat.Kvd)/6 + (adat.Kvb)/3 - (adat.Kvc * (2.0/3.0));
    }
    
    setForces(false);
    update_constraints();
    
    for(uint i=0; i<g.size(); i++){
        g[i].a = g[i].f / g[i].m;
    }
    
    update_trackers();
}

flt collectionGear4NPH::kinetic(){
    flt E=0;
    flt Vfac = dV/box->V()/flt(NDIM);
    for(uint i=0; i<atoms->size(); i++){
        Vec v = (*atoms)[i].v - ((*atoms)[i].x * Vfac);
        E += v.sq() * (*atoms)[i].m;
    }
    return E/2.0;
}

flt collectionGear4NPH::temp(bool minuscomv){
    flt totkinetic2 = 0; // kinetic * 2
    Vec cv = Vec();
    if(minuscomv) cv = comv();
    flt Vfac = dV/box->V()/flt(NDIM);
    for(uint i=0; i<atoms->size(); i++){
            Vec v = (*atoms)[i].v - cv - ((*atoms)[i].x * Vfac);
            totkinetic2 += v.sq() * (*atoms)[i].m;
    }
    
    int ndof = (int) dof();
    if (minuscomv) ndof -= NDIM;
    return totkinetic2 / ndof;
}

void collectionGear4NPH::timestep(){
    uint n=0;
    OriginBox& obox = (OriginBox&) *box;
    atomgroup &g = *atoms;
    
    /// Prediction
    for(uint i=0; i<g.size(); i++){
        g[i].x += g[i].v * dt + g[i].a * (dt*dt/2) + bs[n] * (dt*dt*dt/6);
        g[i].v += g[i].a * dt + bs[n] * (dt*dt/2);
        g[i].a += bs[n] * dt;
        n++;
    }

    flt V = obox.V();
    flt newV = V + dV * dt + ddV * (dt*dt/2) + dddV * (dt*dt*dt/6);
    dV += ddV * dt + dddV * (dt*dt/2);
    ddV += dddV * dt;
    V = obox.resizeV(newV);
    //~ if (abs(V-newV)/abs(newV) > 1e-8){
        //~ printf("V: %.4g; newV: %.4g\n", V, newV);
        //~ assert(abs(V-newV)/abs(newV) < 1e-8);
    //~ }
    
    // Now we set forces and accelerations
    for(uint m=0; m<ncorrec; m++){
        flt interacP = setForcesGetPressure(false);
        flt newP = (interacP + (kinetic()*2.0))/NDIM/V;
        //~ for(git = groups.begin(); git<groups.end(); git++){
            //~ atomgroup &g = **git;
            //~ for(uint i=0; i<g.size(); i++){
                //~ Vec a = g[i].f / g[i].m;
            //~ }
        //~ }
        //~ newP = pressure();
        update_constraints();
        
        /// Correction
        flt Vfac = (-(2.0/NDIM/NDIM)*(dV/V)*(dV/V)) + (ddV/V/NDIM);
        flt newddV = (newP - P) / Q;
        flt correctionV = newddV - ddV;
        newV = V + correctionV * (dt*dt/12);
        dV += correctionV * (5*dt/12);
        ddV = newddV;
        dddV += correctionV / dt;
        
        V = obox.resizeV(newV);
        //~ if (abs(V-newV)/abs(newV) > 1e-8){
            //~ printf("V: %.4g; newV: %.4g\n", V, newV);
            //~ assert(abs(V-newV)/abs(newV) < 1e-8);
        //~ }
        
        
        n=0;
        
        for(uint i=0; i<g.size(); i++){
            Vec a = g[i].f / g[i].m + g[i].x * Vfac;
            Vec correction = a - g[i].a;
            g[i].x += correction * (dt*dt/12);
            g[i].v += correction * (5*dt/12);
            g[i].a = a;
            bs[n] += correction / dt;
            n++;
        }
    }
    update_trackers();
}

vector<sptr<interaction> > collectionGear4NPT::tointerpair(vector<sptr<interactionpairsx> > &v){
    vector<sptr<interaction> > newv = vector<sptr<interaction> >();
    vector<sptr<interactionpairsx> >::iterator it;
    for(it = v.begin(); it<v.end(); ++it){
        newv.push_back((sptr<interaction>) *it);
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
    atoms->resetForces();
    
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it<interactions.end(); ++it){
        sptr<interactionpairsx> inter = boost::dynamic_pointer_cast<interactionpairsx>(*it);
        inter->setForces(*box, &xrpsums);
    }
    if(!seta) return;
    
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].a = (*atoms)[i].f / (*atoms)[i].m;
    }
}

void collectionGear4NPT::timestep(){
    uint n=0;
    OriginBox& obox = (OriginBox&) *box;
    atomgroup &g = *atoms;
    
    const flt c1 = dt, c2 = dt*dt/2, c3 = dt*dt*dt/6;
    
    /// Prediction
    for(uint i=0; i<g.size(); i++){
        g[i].x += xs1[n]*c1 + xs2[n]*c2 + xs3[n]*c3;
        xs1[n] += xs2[n]*c1 + xs3[n]*c2;
        xs2[n] += xs3[n]*c1;
        
        g[i].v += g[i].a*c1 + vs2[n]*c2 + vs3[n]*c3;
        g[i].a += vs2[n]*c1 + vs3[n]*c2;
        vs2[n] += vs3[n]*c1;
        n++;
    }
    
    flt V = obox.V();
    flt newV = V + V1*c1 + V2*c2 + V3*c3;
    V1 += V2*c1 + V3*c2;
    V2 += V3*c1;
    V = obox.resizeV(newV);
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
        V = obox.resizeV(V + Vcorr*d0);
        V1 += Vcorr;
        V2 += Vcorr*d2;
        V3 += Vcorr*d3;
        
        
        n=0;
        for(uint i=0; i<g.size(); i++){
            Vec newx1 = g[i].v + (g[i].x * chi);
            Vec xcorr = newx1 - xs1[n];
            g[i].x += xcorr*d0;
            xs1[n] = newx1;
            xs2[n] += xcorr*d2;
            xs3[n] += xcorr*d3;
            
            Vec newa = g[i].f / g[i].m - g[i].v*chixi;
            Vec vcorr = newa - g[i].a;
            g[i].v += vcorr*d0;
            g[i].a = newa;
            vs2[n] += vcorr*d2;
            vs3[n] += vcorr*d3;
            n++;
        }
    }
    update_trackers();
}

void collectionVerletNPT::resetvhalf(){
    vhalf.resize(atoms->size(), Vec());
    uint n=0;
    for(uint i=0; i<atoms->size(); i++){
        vhalf[n] = (*atoms)[i].v;
        n++;
    }
}

void collectionVerletNPT::timestep(){
    OriginBox& obox = (OriginBox&) *box;
    atomgroup &g = *atoms;
    setForces(false);
    update_constraints();
    
    //~ vector<Vec> lastvhalf = vhalf; // This does copy, right?
    uint ndof = (uint) dof();
    flt xieta = (xidot + eta) * dt/2.0;
    flt Ktot2 = 0, lastKtot2 = 0;
    
    uint n=0;
    for(uint i=0; i<g.size(); i++){
        flt mi = g[i].m;
        g[i].a = g[i].f / mi;
        lastKtot2 += vhalf[n].sq() / mi;
        Vec lastvhalf = vhalf[n];
        vhalf[n] = (vhalf[n]*(1-xieta) + g[i].a*dt)/(1+xieta);
        
        g[i].v = (lastvhalf + vhalf[n])/2.0;
        Ktot2 += vhalf[n].sq() / mi;
        n++;
    }
    
    eta += dt*(Ktot2 - ndof*T)/QT;
    
    flt V = box->V();
    flt newV = obox.resizeV(lastV + (dt*xidot*V));
    flt Verr = newV / (lastV + (dt*xidot*V));
    
    lastV = V;
    curP = ((Ktot2 + lastKtot2)/2 + virial()) / NDIM / V;
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
    flt Vfac1=pow(newV/V, OVERNDIM), Vfac2 = pow(2*newV/(newV+V), OVERNDIM);
    for(uint i=0; i<g.size(); i++){
        g[i].x = g[i].x*Vfac1 + vhalf[n]*dt*Vfac2;
        n++;
    }
    
    update_trackers();
};

void collectionCDBDgrid::reset_velocities(){
    for(uint i=0; i<atoms->size(); i++){
        flt mi = (*atoms)[i].m;
        (*atoms)[i].v = randVec() * sqrt(T/mi);
    }
    reset_events();
};

void collectionCDBDgrid::line_advance(flt deltat){
    if(deltat == 0) return;
    
    //~ assert(deltat > 0);
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x += (*atoms)[i].v * deltat;
    }
    
    curt += deltat;
    update_grid();
};

bool make_event(Box &box, event& e, atomid a, atomid b, flt sigma, flt curt){
    Vec rij = box.diff(a->x, b->x);
    Vec vij = a->v - b->v;
    flt bij = rij.dot(vij);
    if(bij >= 0) return false;
    
    flt bijsq = bij*bij;
    flt vijsq = vij.sq();
    flt underroot = bijsq - vijsq*(rij.sq() - sigma*sigma);
    if(underroot < 0) return false;
    
    flt tau_ij = -(bij + sqrt(underroot))/vijsq;
    // collide immediately if the event already happened (its overlapping)
    if(tau_ij < 0) tau_ij = 0;
    
    e.t = curt + tau_ij;
    if(a.n() > b.n()){
        atomid c = a;
        a = b;
        b = c;
    }
    e.a = a;
    e.b = b;
    return true;
}

void collectionCDBDgrid::update_grid(bool force){
    if(!force && (gridt == curt)) return;
    grid.optimize_widths();
    grid.make_grid();
};

event collectionCDBDgrid::next_event(atomid a){
    event e;
    flt vmag = a->v.mag();
    if(!isfinite(vmag) or vmag <= 0) vmag = sqrt(T/a->m) * edge_epsilon;
    flt epsilon_t = atomsizes[a.n()]*edge_epsilon/vmag;
    assert(epsilon_t > 0);
    
    e.t = curt + grid.time_to_edge(*a) + epsilon_t;
    e.a = a;
    e.b = a;
        
    for(Grid::pair_iter j=grid.pairs(a); j!=j.end(); ++j){
        event e2;
        flt sigma = (atomsizes[a.n()] + atomsizes[(*j).n()])/2.0;
        if(make_event(*box, e2, *j, a, sigma, curt)
                and e2 < e){
            e = e2;
        };
    }
    assert(e.t >= curt);
    return e;
}

void collectionCDBDgrid::reset_events(bool force){
    //~ std::cerr << "CDBD::reset_events...\n";
    events.clear();
    //~ std::cerr << "CDBD::reset_events cleared...\n";
    update_grid(force);
    //~ std::cerr << "CDBD::reset_events updated grid...\n";
    
    for(uint i=0; i<atoms->size(); i++){
        //~ std::cerr << "CDBD::reset_events" << i << "\n";
        events.insert(next_event(atoms->get_id(i)));
    };
    //~ std::cerr << "CDBD::reset_events done.\n";
};

inline void collide(Box &box, atom& a, atom &b){
    assert(&a != &b);
    Vec rij = box.diff(a.x, b.x);
    flt rmag = rij.mag();
    Vec rhat = rij / rmag;
    
    flt rvi = rhat.dot(a.v);
    flt rvj = rhat.dot(b.v);
    if(rvi >= rvj) return; // they're already moving away
    Vec p = rhat * (2 * a.m * b.m * (rvj - rvi) / (a.m + b.m));
    a.v += p / a.m;
    b.v -= p / b.m;
};

bool collectionCDBDgrid::take_step(flt tlim){
    // TODO: more efficient looping and merging
    
    // If we have a limit, and we've already passed it, stop.
    if((tlim > 0) && (tlim <= curt)) return false;
    
    if(events.empty()) reset_events();
    
    assert(atoms->size() > 0);
    assert(atomsizes.size() > 0);
    if(events.empty()){
        assert(!events.empty());
        //~ // TODO: this should check time to leave box, and only go that far
        //~ line_advance(tlim - curt);
        //~ curt = tlim;
        //~ return false;
    }
    event e = *(events.begin());
    if ((tlim > 0) & (e.t > tlim)){
        // if we have a limit, and the next event is farther in the 
        // future than that limit, then we just move everyone forward and
        // no collisions happen
        line_advance(tlim - curt);
        return false;
    };
    events.erase(e);
    
    // move everyone forward
    line_advance(e.t - curt);
    
    // Remove all "bad" scheduled events (involving these two atoms),
    // and mark which atoms need rescheduling
    vector<atomid> badatoms(1);
    badatoms[0] = e.a;
    
    bool fakecollision = (e.a == e.b);
    
    // if e.a == e.b, that just means that e.a is entering a new box,
    // and we don't need to worry about any other atoms
    if(!fakecollision){
        // collide our two atoms (i.e. have them bounce)
        collide(*box, *(e.a), *(e.b));
        badatoms.push_back(e.b);
        
        set<event>::iterator eit, eit2;
        eit = events.begin();
        while (eit != events.end()){
            if((eit->a == e.a) || (eit->a == e.b)){
                badatoms.push_back(eit->b);
            } else if((eit->b == e.a) || (eit->b == e.b)){
                badatoms.push_back(eit->a);
            } else {
                // no matches, carry on
                ++eit;
                continue;
            }
            
            // need to remove that event
            eit2 = eit;
            ++eit;
            events.erase(eit2);
        };
    }
    
    for(uint i=0; i<badatoms.size(); i++){
        event e = next_event(badatoms[i]);
        events.insert(e);
    };
    return true;
};

void collectionCDBDgrid::timestep() {
    reset_velocities();
    flt newt = curt + dt;
    while(take_step(newt)){};
    update_trackers();
};

void collectionCDBD::reset_velocities(){
    for(uint i=0; i<atoms->size(); i++){
        flt mi = (*atoms)[i].m;
        (*atoms)[i].v = randVec() * sqrt(T/mi);
    }
    reset_events();
};

void collectionCDBD::line_advance(flt deltat){
    for(uint i=0; i<atoms->size(); i++){
        (*atoms)[i].x += (*atoms)[i].v * deltat;
    }
    
    curt += deltat;
};

void collectionCDBD::reset_events(){
    events.clear();
    
    for(uint i=1; i<atoms->size(); i++){
        for(uint j=0; j<i; j++){
            flt sigma = (atomsizes[i] + atomsizes[j]) / 2;
            
            //~ std::cerr << i << "-" << j;
            event e;
            if(make_event(*box, e, atoms->get_id(j), atoms->get_id(i), sigma, curt)){
                events.insert(e);
                //~ std::cerr << " MADE!";
            };
            //~ std::cerr << ", ";
        }
    };
    //~ std::cerr << " Done, events " << events.size() << "\n";
};

bool collectionCDBD::take_step(flt tlim){
    // TODO: more efficient looping and merging
    
    // If we have a limit, and we've already passed it, stop.
    if((tlim > 0) && (tlim <= curt)) return false;
    
    //~ std::cerr << "take_step: events " << events.size() << "\n";
    if(events.empty()) reset_events();
    
    assert(atoms->size() > 0);
    assert(atomsizes.size() > 0);
    if(events.empty() == 0){
        // TODO: this should check time to leave box, and only go that far
        line_advance(tlim - curt);
        curt = tlim;
        return false;
    }
    //~ std::cerr << "take_step: started, events " << events.size() << "\n";
    event e = *(events.begin());
    if ((tlim > 0) & (e.t > tlim)){
        // if we have a limit, and the next event is farther in the 
        // future than that limit, then we just move everyone forward and
        // no collisions happen
        line_advance(tlim - curt);
        curt = tlim;
        return false;
    };
    events.erase(e);
    
    // move everyone forward
    line_advance(e.t - curt);
    curt = e.t;
    //~ std::cerr << "take_step: advanced.\n";
    
    // collide our two atoms (i.e. have them bounce)
    //~ std::cerr << "take_step: colliding " << e.a.n() << " - " << e.b.n() << "\n";
    collide(*box, *(e.a), *(e.b));
    //~ std::cerr << "take_step: collided " << e.a.n() << " - " << e.b.n() << "\n";
    
    // Remove all "bad" scheduled events (involving these two atoms),
    // and mark which atoms need rescheduling
    set<atomid> badatoms;
    set<atomid>::iterator ait;
    badatoms.insert(e.a);
    badatoms.insert(e.b);
    
    set<event>::iterator eit, eit2;
    eit = events.begin();
    while (eit != events.end()){
        if((eit->a == e.a) || (eit->a == e.b)){
            badatoms.insert(eit->b);
        } else if((eit->b == e.a) || (eit->b == e.b)){
            badatoms.insert(eit->a);
        } else {
            // no matches, carry on
            ++eit;
            continue;
        }
        
        
        //~ std::cerr << "take_step: removing event...";
        // need to remove that event
        eit2 = eit;
        ++eit;
        events.erase(eit2);
        //~ std::cerr << "Removed.\n";
    };
    
    
    //~ std::cerr << "take_step: new events (" << badatoms.size() << ")...";
    // Make new events, add them to the list
    vector<event> newevents(badatoms.size());
    //~ std::cerr << "made events vector...";
    vector<bool> new_filled(badatoms.size(), false);
    //~ std::cerr << "made new_filled.\n";
    //~ std::cerr << "Starting loop.\n";
    for(uint i=0; i<atoms->size(); i++){
        uint bi=0;
        //~ std::cerr << "Starting badatoms loop.\n";
        for(ait=badatoms.begin(); ait!=badatoms.end(); ++ait){
            //~ cerr << "Testing " << ait->n() << " - " << n  << "(" << bi << "), ";
            flt sigma = (atomsizes[ait->n()] + atomsizes[i])/2.;
            //~ cerr << "sigma " << sigma << "\n";
            event testevent;
            if(make_event(*box, testevent, *ait, atoms->get_id(i), sigma, curt)){
                //~ cerr << "Made good event\n";
                if(!new_filled[bi] || (testevent < newevents[bi])){
                    newevents[bi] = testevent;
                    new_filled[bi] = true;
                    //~ cerr << "Its the best event\n";
                } else {
                    //~ cerr << "event ignored\n";
                }
            } else {
                //~ cerr << "event not made.\n";
            };
            bi++;
        }
    };
    
    //~ std::cerr << "take_step: made new events\n";
    
    for(uint bi=0; bi<newevents.size(); bi++){
        if(new_filled[bi]){
            //~ std::cerr << "take_step: inserting collision between " << 
                //~ newevents[bi].a.n() << " - " << newevents[bi].b.n() << "\n";
            events.insert(newevents[bi]);
        }
    };
    //~ std::cerr << "take_step: Done.\n";
    return true;
};

void collectionCDBD::timestep() {
    reset_velocities();
    flt newt = curt + dt;
    while(take_step(newt)){};
    update_trackers();
};
