#include "box.hpp"

#ifndef TRACKERS_H
#define TRACKERS_H

#include <set>
#include <map>

class statetracker {
    public:
        virtual void update(Box &box) = 0;
        virtual ~statetracker(){};
};

class pairlist {
    protected:
        map<const atomid, set<atomid> > pairs;
    public:
        //~ pairlist(atomgroup *group);
        pairlist(){};
        
        inline void ensure(const atomid a){
            pairs.insert(std::pair<atomid, set<atomid> >(a, set<atomid>()));
        }
        inline void ensure(vector<atomid> ps){
            vector<atomid>::iterator it;
            for(it=ps.begin(); it != ps.end(); it++) ensure(*it);
        }
        inline void ensure(atomgroup &group){
            for(uint i=0; i<group.size(); i++) ensure(group.get_id(i));
        }
        
        inline bool has_pair(atomid a1, atomid a2){
            if(a1 > a2) return pairs[a1].count(a2) > 0;
            else return pairs[a2].count(a1) > 0;
        }
        
        inline void add_pair(atomid a1, atomid a2){
            //~ cout << "pairlist ignore " << a1.n() << '-' << a2.n() << "\n";
            if(a1 > a2){pairs[a1].insert(a2);}
            else{pairs[a2].insert(a1);};
        }
        
        inline void erase_pair(atomid a1, atomid a2){
            if(a1 > a2){pairs[a1].erase(a2);}
            else{pairs[a2].erase(a1);};
        }
        
        inline set<atomid> get_pairs(const atomid a){ensure(a); return pairs[a];};
        
        // for iterating over neighbors
        inline set<atomid>::iterator begin(const atomid a){return pairs[a].begin();};
        inline set<atomid>::iterator end(const atomid a){return pairs[a].end();};
        
        inline uint size() const { uint N=0; for(
            map<const atomid, set<atomid> >::const_iterator it=pairs.begin();
            it != pairs.end(); it++) N+= (uint) it->second.size();
            return N;
        };
        
        void clear();
};

class neighborlist : public statetracker{
    //maintains a Verlet list of "neighbors": molecules within a 
    // 'skin radius' of each other.
    // note that molecules are counted as neighbors if any point within
    // their molecular radius is within a 'skin radius' of any point
    // within another molecule's molecular radius.
    //
    // update(false) should be called frequently; it checks (O(N)) if
    // any two molecules might conceivably overlap by more than a critical
    // distance, and if so, it updates all the neighbor lists.
    
    // the <bool areneighbors()> function returns whether two molecules
    // are neighbors, and the begin(n), end(n) allow for iterating over
    // the neighbor lists
    protected:
        sptr<Box> box;
        flt critdist, skinradius;
        metagroup atoms;
        vector<idpair> curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        atomid get_id(atom* a);
        bool ignorechanged; // if true, forces a full check on next update
        //~ bool checkneighbors(const uint n, const uint m) const;
        // this is a full check
    public:
        neighborlist(sptr<Box> box, const flt innerradius, const flt outerradius);
        neighborlist(sptr<Box> box, atomgroup &vec, const flt innerradius, 
        const flt outerradius, pairlist ignore = pairlist());
        void update(Box &newbox){assert(&newbox == box.get()); update_list(false);};
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first
        
        inline uint which(){return updatenum;};
        inline uint numpairs(){return (uint) curpairs.size();};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b); ignorechanged=true;};
        void ignore(atom*, atom*);
        atomid add(atom* a){
            atomid id = atoms.get_id(a);
            if(id != NULL) return id;
            atoms.add(a);
            assert(lastlocs.size() == atoms.size() - 1);
            lastlocs.push_back(a->x);
            id = atoms.get_id(a);
            ignorechanged = true;
            return id;
        }
        
        inline void changesize(flt inner, flt outer){
            critdist = inner; skinradius = outer; update_list(true);};
        inline void changesize(flt ratio){
            skinradius = critdist*ratio; update_list(true);};

        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint size() const{return atoms.size();};
        inline vector<idpair>::iterator begin(){return curpairs.begin();};
        inline vector<idpair>::iterator end(){return curpairs.end();};
        inline idpair get(uint i){return curpairs[i];};
        //~ inline vector<idpair> getpairs(){return vector<idpair>(curpairs);};
        ~neighborlist(){};
};

inline neighborlist* neighborlistL(sptr<Box> box, const double innerradius, const double outerradius){
    return new neighborlist(box, innerradius, outerradius);
};

#endif
