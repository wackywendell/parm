#include "collection.hpp"

collection::collection(vector<atomgroup*> gs, vector<interaction*> is)
        : groups(gs), interactions(is){
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

flt collection::Energy(){
    flt E=0;
    vector<interaction*>::iterator it;
    for(it = interactions.begin(); it<interactions.end(); it++){
        interaction &inter = **it;
        E += inter.energy();
    }
    
    E += kinetic();
    return E;
};

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

void collection::timestep(const flt dt){
    vector<atomgroup*>::iterator git;
    
    // this will slow it down a bit, but it would be nice to have AUAlocs
    // always correct
    for(git = groups.begin(); git<groups.end(); git++)
        (**git).vverlet1(dt);
    // clear old forces
    // just check all
    setForces();
    
    for(git = groups.begin(); git<groups.end(); git++)
        (**git).setAccel();
    
    // use force info to update velocities
    // [v(t+dt/2) -> v(t+dt) with f(t+dt)]
    for(git = groups.begin(); git<groups.end(); git++)
        (**git).vverlet2(dt);
    
}

collectionSol::collectionSol(flt damp, flt sigma, 
        vector<atomgroup*> groups,vector<interaction*> interactions) :
        collection(groups, interactions), gauss(sigma), damping(damp){
    if(damping == 0 and sigma == 0){
        //*timestep = *collection::timestep; // use the general form
    } 
};

void collectionSol::timestep(const flt dt){
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &m = **git;
        for(uint i=0; i<m.N(); i++){
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
        for(uint i=0; i<m.N(); i++){
            m[i].a = m[i].f / m.getmass(i) + gauss.generate();
            
            // step 5: finish v and a
            m[i].v += m[i].a * (dt/2);
            m[i].a -= m[i].v * damping;
        }
    }
}
