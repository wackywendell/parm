#include "box.hpp"

#ifndef TRACKERS_H
#define TRACKERS_H

#include <set>
#include <map>
#include "assert.h"

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

class GridIterator;
class GridPairedIterator;

class Grid {
    public:
        sptr<OriginBox> box;
        sptr<atomgroup> atoms;
        // minwidth is minimum size of a box, in real units
        // goalwidth is how many atoms per box (goal)
        flt minwidth, goalwidth;
        uint widths[NDIM];
        vector<set<atomid> > gridlocs;
        
        vector<uint> neighbors(uint i);
        uint get_loc(Vec v, Vec bsize);
        
    public:
        typedef GridIterator iterator;
        typedef GridPairedIterator pair_iter;
        
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, const uint width=1);
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, vector<uint> width);
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, const flt minwidth, const flt goalwidth);
        
        //void add(atom& a){atoms.add(&a);};
        //void add(atomgroup& a){for(uint i=0; i<a.size(); i++) add(a[i]);};
        
        void optimize_widths();
        void make_grid();
        
        friend class GridIterator;
        friend class GridPairedIterator;
        iterator begin();
        iterator end();
        
        pair_iter pairs(atomid a);
        flt time_to_edge(atom& a);
        flt time_to_edge(uint i){return time_to_edge((*atoms)[i]);}
        
        vector<idpair> allpairs();
        vector<atomid> allpairs(atomid a);
        
        uint numcells(uint i){assert(i<NDIM); return widths[i];}
        #ifdef VEC2D
        uint numcells(){ return widths[0] * widths[1];};
        #else
        uint numcells(){ return widths[0] * widths[1] * widths[2];};
        #endif
};

class GridPairedIterator {
    // For iterating over all pairs of a single atom
    protected:
        Grid & grid;
        atomid atom1;
        
        vector<uint> neighbor_cells;
        vector<uint>::iterator cellnum2;
        set<atomid> *cell2;
        set<atomid>::iterator atom2;
        
        // returns "successful increment", e.g. cell1 points to an element
        bool increment_cell2();
        bool increment_atom2();
    
    public:
        typedef set<atomid>::iterator end_type;
    
        GridPairedIterator(Grid & grid, atomid a);
        
        GridPairedIterator& operator++();
        atomid operator*(){return *atom2;};
        bool operator==(const GridPairedIterator &other);
        
        bool operator!=(const GridPairedIterator &other){
            return !(*this == other);
        };
        
        bool operator==(const end_type other){
            return (atom2 == other);
        }
        bool operator!=(const end_type other){
            return !(atom2 == other);
        };
        
        end_type end(){return grid.gridlocs[neighbor_cells.back()].end();};
};

class GridIterator {
    // For iterating over all possible pairs
    protected:
        Grid & grid;
        vector<set<atomid> >::iterator cell1;
        set<atomid>::iterator atom1;
        
        vector<uint> neighbor_cells;
        vector<uint>::iterator cellnum2;
        set<atomid> *cell2;
        set<atomid>::iterator atom2;
    
        // returns "successful increment", e.g. cell1 points to an element
        bool increment_cell1();
        bool increment_atom1();
        bool increment_cell2();
        bool increment_atom2();
    
    public:
        GridIterator(Grid & grid);
        GridIterator(Grid & grid, vector<set<atomid> >::iterator cell1);
        
        GridIterator& operator++();
        idpair operator*(){return idpair(*atom1, *atom2);};
        bool operator==(const GridIterator &other);
        
        bool operator!=(const GridIterator &other){
            return !(*this == other);
        };
};

#endif
