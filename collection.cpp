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

flt collection::Temp(){
    flt totatoms = 0;
    flt totkinetic = 0;
    Vec v = comv();
    vector<atomgroup*>::iterator git;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &group = **git;
        totkinetic += group.kinetic(v);
        totatoms += group.N();
    }
    return totkinetic * 2 / (3*totatoms-3);
}

Vec collection::com(){
    vector<atomgroup*>::iterator git;
    flt mass = 0;
    Vec totcom = Vec(0,0,0);
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.N(); i++){
            flt curmass = g.getmass(i);
            mass += curmass;
            totcom += g[i].x * curmass;
        }
    }
    return totcom / mass;
}

Vec collection::comv(){
    vector<atomgroup*>::iterator git;
    flt mass = 0;
    Vec totcom = Vec(0,0,0);
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.N(); i++){
            flt curmass = g.getmass(i);
            mass += curmass;
            totcom += g[i].v * curmass;
        }
    }
    return totcom / mass;
}

flt collection::gyradius(){
    vector<atomgroup*>::iterator git;
    Vec avgr = Vec(0,0,0);
    flt N = 0;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.N(); i++){
            avgr += g[i].x;
            N++;
        }
    }
    avgr /= N; // now avgr is the average location, akin to c.o.m.
    flt Rgsq = 0;
    for(git = groups.begin(); git<groups.end(); git++){
        atomgroup &g = **git;
        for(uint i = 0; i<g.N(); i++) Rgsq += (g[i].x - avgr).sq();
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

collectionSol::collectionSol(flt damp, flt sigma, 
        vector<atomgroup*> groups,vector<interaction*> interactions) :
        collection(groups, interactions), gauss(sigma), damping(damp){
    if(damping == 0 and sigma == 0){
        //*timestep = *collection::timestep; // use the general form
    } 
};

void collectionSol::timestep(const flt dt){
    //~ cout << "damping " << damping << "\n";
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
            Vec g = gauss.generate();
            //~ cout << "g " << g << "\n";
            m[i].a = (m[i].f + g) / m.getmass(i);
            
            // step 5: finish v and a
            m[i].v += m[i].a * (dt/2);
            m[i].a -= m[i].v * damping;
        }
    }
}
