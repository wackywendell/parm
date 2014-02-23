#include "constraints.hpp"


ContactTracker::ContactTracker(sptr<Box> box, sptr<atomgroup> atoms, vector<flt> dists) :
    atoms(atoms), dists(dists), contacts(), breaks(0), formations(0),
        incontact(0){
    //~ cout << "Making contact tracker." << endl;
    uint N = atoms->size();
    dists.resize(N);
    contacts.resize(N);
    for(uint i=0; i<N; i++){
        contacts[i].resize(i, false);
    }
    update(*box);
    breaks = 0;
    formations = 0;
};

void ContactTracker::update(Box &box){
    //~ cout << "Contact tracker update." << endl;
    uint N = atoms->size();
    incontact = 0;
    for(uint i=0; i<N; i++){
        Vec ri = atoms->get(i)->x;
        for(uint j=0; j<i; j++){
            Vec rj = atoms->get(j)->x;
            Vec dr = box.diff(ri,rj);
            bool curcontact = (dr.mag() <= ((dists[i] + dists[j])/2));
            if(curcontact && (!contacts[i][j])) formations++;
            else if(!curcontact && contacts[i][j]) breaks++;
            if(curcontact) incontact++;
            
            contacts[i][j] = curcontact;
        }
    }
    //~ cout << "Contact tracker update done." << endl;
};

void EnergyTracker::update(Box &box){
    if(nskipped + 1 < nskip){
        nskipped += 1;
        return;
    }
    
    nskipped = 0;
    uint Natoms = atoms->size();
    flt curU = 0, curK = 0;
    for(uint i=0; i<Natoms; i++){
        atom &curatom = *(atoms->get(i));
        curK += curatom.v.sq() * curatom.m / 2;
    }
    
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it != interactions.end(); it++){
        curU += (*it)->energy(box);
    }
    
    curU -= U0;
    Ks += curK;
    Us += curU;
    Es += curK + curU;
    Ksq += curK*curK;
    Usq += curU*curU;
    Esq += (curK + curU)*(curK + curU);
    N++;
};

void EnergyTracker::setU0(Box &box){
    flt curU = 0;
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it != interactions.end(); it++){
        curU += (*it)->energy(box);
    }
    setU0(curU);
};


RsqTracker1::RsqTracker1(atomgroup& atoms, uint skip, Vec com) :
            pastlocs(atoms.size(), Vec()), rsqsums(atoms.size(), 0.0),
            rsqsqsums(atoms.size(), 0.0), skip(skip), count(0){
    for(uint i = 0; i<atoms.size(); i++){
        pastlocs[i] = atoms[i].x - com;
    };
};
        
void RsqTracker1::reset(atomgroup& atoms, Vec com){
    pastlocs.resize(atoms.size());
    rsqsums.assign(atoms.size(), 0);
    rsqsqsums.assign(atoms.size(), 0);
    for(uint i = 0; i<atoms.size(); i++){
        pastlocs[i] = atoms[i].x - com;
    };
    count = 0;
};

bool RsqTracker1::update(Box& box, atomgroup& atoms, uint t, Vec com){
    if(t % skip != 0) return false;
    
    for(uint i = 0; i<atoms.size(); i++){
        //flt dist = box.diff(atoms[i].x, pastlocs[i]).sq();
        Vec r = atoms[i].x - com;
        flt dist = (r - pastlocs[i]).sq();
        // We don't want the boxed distance - we want the actual distance moved!
        rsqsums[i] += dist;
        rsqsqsums[i] += dist*dist;
        pastlocs[i] = r;
    };
    count += 1;
    return true;
};

vector<flt> RsqTracker1::rsq_mean(){
    vector<flt> means(rsqsums.size(), 0.0);
    for(uint i=0; i<rsqsums.size(); i++){
        means[i] = rsqsums[i] / count;
    }
    return means;
};

vector<flt> RsqTracker1::rsq_var(){
    vector<flt> vars(rsqsums.size(), 0.0);
    for(uint i=0; i<rsqsums.size(); i++){
        flt mn = (rsqsums[i])/count;
        flt sqmn = (rsqsqsums[i])/count;
        vars[i] = sqmn - mn*mn;
    }
    return vars;
};
    

RsqTracker::RsqTracker(sptr<atomgroup> atoms, vector<uint> ns, bool usecom) : 
            atoms(atoms), curt(0), usecom(usecom){
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<uint>::iterator n=ns.begin(); n!=ns.end(); n++){
        singles.push_back(RsqTracker1(*atoms, *n, com));
    }
};

void RsqTracker::update(Box &box){
    curt++;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); it++){
        it->update(box, *atoms, curt, com);
    }
};

void RsqTracker::reset(){
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); it++){
        it->reset(*atoms, com);
    }
};

vector<vector<flt> > RsqTracker::means(){
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); it++){
        vals.push_back(it->rsq_mean());
    }
    return vals;
};

vector<vector<flt> > RsqTracker::vars(){
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); it++){
        vals.push_back(it->rsq_var());
    }
    return vals;
};

vector<flt> RsqTracker::counts(){
    vector<flt> vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); it++){
        vals.push_back(it->get_count());
    }
    return vals;
};
