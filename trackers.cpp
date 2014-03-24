#include "trackers.hpp"

neighborlist::neighborlist(sptr<Box> box, const flt innerradius, const flt outerradius) :
                box(box), critdist(innerradius), skinradius(outerradius),
                atoms(){};

neighborlist::neighborlist(sptr<Box> box, atomgroup &group, const flt innerradius,
            const flt outerradius, pairlist ignore) :
                box(box), critdist(innerradius), skinradius(outerradius),
                atoms(), ignorepairs(ignore), ignorechanged(false){
    lastlocs.resize(group.size());
    for(uint i=0; i<group.size(); i++){
        atomid a = group.get_id(i);
        ignorepairs.ensure(a);
        atoms.add(a.pointer());
        lastlocs[i] = a.x();
    }
    update_list(true);
    updatenum = 1;
};

atomid neighborlist::get_id(atom* a){
    for(uint i=0; i < atoms.size(); i++)
        if (atoms.get(i) == a) return atoms.get_id(i);
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
    
    if(not force and not ignorechanged){ // check if we need to update
        for(uint i=0; i < atoms.size(); i++){
            atom &atm = atoms[i];
            curdist = (atm.x - lastlocs[i]).mag();
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
            //~ if(skinradius - critdist > 5)
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
    for(uint i=0; i<atoms.size(); i++){
        atomid a1=atoms.get_id(i);
        lastlocs[i] = a1.x();
        for(uint j=0; j<i; j++){
            atomid a2=atoms.get_id(j);
            if (ignorepairs.has_pair(a1, a2)) continue;
            if(box->diff(a1.x(), a2.x()).mag() < skinradius)
                curpairs.push_back(idpair(a1, a2));
        }
    }
    //~ cout << "neighborlist::update_list:: done.\n";
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

