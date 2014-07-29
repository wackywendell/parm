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

class NLPairIter;

class neighborlist : public statetracker{
    //maintains a Verlet list of "neighbors": molecules within a 
    // 'skin radius' of each other.
    // note that molecules are counted as neighbors if the distance
    // between them is less than (r1 + r2 + R), where r1 and r2 are the 
    // radii, and R is the "skin radius".
    //
    // update(false) should be called frequently; it checks (O(N)) if
    // any two molecules might conceivably overlap by more than a critical
    // distance, and if so, it updates all the neighbor lists.
    
    protected:
        sptr<Box> box;
        flt skin;
        SubgroupBool atoms;
        vector<flt> diameters;
        vector<vector<atomid> > curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        bool ignorechanged; // if true, forces a full check on next update
    public:
        typedef NLPairIter iterator;
        friend class NLPairIter;
        
        neighborlist(sptr<Box> box, sptr<atomvec> atoms, const flt skin);
        void update(Box &newbox){assert(&newbox == box.get()); update_list(false);};
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first
        
        atomvec& vec(){return atoms.vec();};
        inline uint which(){return updatenum;};
        //~ inline uint numpairs(){return (uint) curpairs.size();};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b); ignorechanged=true;};
        void add(atomid a, flt diameter){
            atoms.add(a);
            //assert(diameters.size() == atoms.size() - 1);
            uint avsize = atoms.vec().size();
            if(diameters.size() < avsize) diameters.resize(avsize, -1);
            if(lastlocs.size() < avsize) lastlocs.resize(avsize);
            diameters[a.n()] = diameter;
            lastlocs[a.n()] = a->x;
            ignorechanged = true;
        }
        
        inline uint ignore_size() const{return ignorepairs.size();};
        inline uint size() const{return atoms.size();};
        
        flt dist_to_half_edge(atomid a){
            flt curdist = (a->x - lastlocs[a.n()]).mag();
            return skin/2 - curdist;
        };
        
        iterator begin();
        iterator end();
        
        vector<idpair> pair_list();
        vector<vector<atomid> > pair_lists(){return curpairs;};
        
        vector<atomid>& near(atomid a){
			assert(a.n() < curpairs.size());
			return curpairs[a.n()];
		};
        
        ~neighborlist(){};
};

class NLPairIter {
	protected:
		neighborlist& nlist;
		vector<vector<atomid> >::iterator curlist;
        vector<atomid>::iterator curloc;
        uint n1;
        
        // move to the next (including current) location which is not empty
        void find_pair(bool force_advance=false);
        
    public:
		NLPairIter(neighborlist& nlist) :
			nlist(nlist), curlist(nlist.curpairs.begin()), 
			curloc(curlist->begin()), n1(0){find_pair(false);};
		NLPairIter(neighborlist& nlist,
			vector<vector<atomid> >::iterator curlist,
			vector<atomid>::iterator curloc, uint n1) :
			nlist(nlist), curlist(curlist), 
			curloc(curloc), n1(n1){};
		
		NLPairIter& operator++(){find_pair(true); return *this;};
        idpair operator*(){return idpair(nlist.atoms.get_id(n1), *curloc);};
        bool operator==(const NLPairIter &other);
        
        bool operator!=(const NLPairIter &other){
            return !(*this == other);
        };
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
