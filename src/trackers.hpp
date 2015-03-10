#include "box.hpp"

#ifndef TRACKERS_H
#define TRACKERS_H

#include <set>
#include <map>
#include "assert.h"

/**
The general interface for a "tracker", a class that needs to be called every timestep.

statetracker is used as a base for some "helper" classes, like neighbor lists (\ref neighborlist),
as well as for many statistics.
*/
class statetracker {
    public:
        /** This function is called once per timestep, when particles are in their set position.*/
        virtual void update(Box &box) = 0;
        virtual ~statetracker(){};
};

/** A mapping of \ref atom -> [list of \ref atom], used by \ref neighborlist to keep track
of what atoms are near what other atoms.

Pairs are assumed to be symmetric, i.e. either `(a1, a2)` is stored or `(a2, a1)` is stored, and
only one will be returned through the iterator.
*/
class pairlist {
    protected:
        map<const atomid, set<atomid> > pairs;
    public:
        //~ pairlist(atomgroup *group);
        pairlist(){};

        /** Ensure that a given \ref atom is in the overall list. */
        inline void ensure(const atomid a){
            pairs.insert(std::pair<atomid, set<atomid> >(a, set<atomid>()));
        }
        /** Ensure that every \ref atom in input is in the overall list. */
        inline void ensure(vector<atomid> ps){
            vector<atomid>::iterator it;
            for(it=ps.begin(); it != ps.end(); ++it) ensure(*it);
        }
        /** Ensure that every \ref atom in input is in the overall list. */
        inline void ensure(atomgroup &group){
            for(uint i=0; i<group.size(); i++) ensure(group.get_id(i));
        }
        /** Check if a given pair is is the list.

        Note that `has_pair(a1, a2) == has_pair(a2, a1)`.
        */
        inline bool has_pair(atomid a1, atomid a2){
            if(a1 > a2) return pairs[a1].count(a2) > 0;
            else return pairs[a2].count(a1) > 0;
        }

        /**
        Add a pair to the list. If already in the list, nothing changes.
        */
        inline void add_pair(atomid a1, atomid a2){
            //~ cout << "pairlist ignore " << a1.n() << '-' << a2.n() << "\n";
            if(a1 > a2){pairs[a1].insert(a2);}
            else{pairs[a2].insert(a1);};
        }


        /**
        Remove a pair from the list. If not already in the list, nothing changes.
        */
        inline void erase_pair(atomid a1, atomid a2){
            if(a1 > a2){pairs[a1].erase(a2);}
            else{pairs[a2].erase(a1);};
        }

        /**
        Get all atoms neighboring an atom *that are "lesser"* than the atom.

        Because of the way neighbors are stored, there is no easy way to actually get *all*
        neighbors of an atom, but if `get_pairs` is called successively on all atoms, it will return
        each and every pair exactly once (including "false positives").
        */
        inline set<atomid> get_pairs(const atomid a){ensure(a); return pairs[a];};

        /** for iterating over neighbors */
        inline set<atomid>::iterator begin(const atomid a){return pairs[a].begin();};
        /** for iterating over neighbors */
        inline set<atomid>::iterator end(const atomid a){return pairs[a].end();};

        /** for iterating over neighbors */
        inline uint size() const { uint N=0; for(
            map<const atomid, set<atomid> >::const_iterator it=pairs.begin();
            it != pairs.end(); ++it) N+= (uint) it->second.size();
            return N;
        };

        /** Clear the *inner* lists. */
        void clear();
};

/**
Maintains a Verlet list of "neighbors": molecules within a 'skin radius' of each other. When the
list is updated (`update_list(true)`), molecules are counted as neighbors if the distance between
them is \f$d_{ij} < r_i + r_j + R\f$, where \f$r_i\f$ and \f$r_j\f$ are the radii, and \f$R\f$ is
the "skin radius". On "normal" updates (`update(box)` or `update_list(false)`, All atoms fulfilling
the condition \f$d_{ij} < r_i + r_j\f$ are guaranteed to be in the list, but there may be some false
positives.

`update_list(false)` should be called frequently (i.e. every timestep); it checks (\f$O(N)\f$) if
any two molecules might conceivably overlap by more than a critical distance, and if so, it updates
all the neighbor lists.

The neighborlist has an inherent tradeoff, set by the skin radius: If the skin is large, then
full updates (which are slow) will be needed less often, but more particles will be considered
"neighbors" even when they are not within \f$d_{ij} < r_i + r_j\f$, i.e. there will be more false
positives.
*/
class neighborlist : public statetracker{
    protected:
        sptr<Box> box;
        flt skin;
        subgroup atoms;
        vector<flt> diameters;
        vector<idpair> curpairs;
        pairlist ignorepairs;
        vector<Vec> lastlocs;
        uint updatenum;
        bool ignorechanged; /*< if true, forces a full check on next update */
    public:
        neighborlist(sptr<Box> box, sptr<atomvec> atoms, const flt skin);
        /**
        Calls `update_list(false)`; exists as a requirement of \ref statetracker.
        */
        void update(Box &newbox){assert(&newbox == box.get()); update_list(false);};
        /**
        Update the neighbor list based on current positions. If `force`, the list is updated
        immediately \f$O(N^2)\f$; otherwise, neighborlist checks if an update is necessary
        \f$O(N)\f$, and updates the list if so.
        */
        bool update_list(bool force = true);
        // if force = false, we check if updating necessary first

        atomvec& vec(){return atoms.vec();};
        inline uint which(){return updatenum;};
        inline uint numpairs(){return (uint) curpairs.size();};
        inline void ignore(atomid a, atomid b){ignorepairs.add_pair(a,b); ignorechanged=true;};
        void add(atomid a, flt diameter){
            atoms.add(a);
            //assert(diameters.size() == atoms.size() - 1);
            diameters.push_back(diameter);
            //assert(lastlocs.size() == atoms.size() - 1);
            lastlocs.push_back(a->x);
            ignorechanged = true;
        }

        /** Number of pairs currently being ignored */
        inline uint ignore_size() const{return ignorepairs.size();};
        /** Number of atoms in list */
        inline uint size() const{return atoms.size();};
        inline vector<idpair>::iterator begin(){return curpairs.begin();};
        inline vector<idpair>::iterator end(){return curpairs.end();};
        inline idpair get(uint i){
            //assert(i<curpairs.size());
            return curpairs[i];
        };
        //~ inline vector<idpair> getpairs(){return vector<idpair>(curpairs);};
        ~neighborlist(){};
};

class GridIterator;
class GridPairedIterator;

/**
A fast algorithm for finding all pairs of neighboring atoms.

A Grid splits the OriginBox into "cells", where no cell is narrower than any atom. For any given atom,
all atoms within the same cell or any neighboring cell are considered "neighbors".
*/
class Grid {
    public:
        sptr<OriginBox> box;
        sptr<atomgroup> atoms;

        /**
        minwidth is minimum size of a box, in real units.

        This is used to calculate the widths of new boxes in `optimize_widths()`.

        May be 0 or negative, to indicate that `optimize_widths()` should never change the box sizes.
        */
        flt minwidth;
        flt goalwidth; /**< goalwidth is how many atoms per box (goal). See minwidth. */
        uint widths[NDIM]; /**< Number of box divisions per dimension. */
        vector<set<atomid> > gridlocs; /**< Atoms in each grid location. */

        vector<uint> neighbors(uint i);
        uint get_loc(Vec v, Vec bsize);

    public:
        typedef GridIterator iterator;
        typedef GridPairedIterator pair_iter;

        /**
        /param width Number of divisions along all axes
        */
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, const uint width=1);

        /**
        /param widths Number of divisions along each axis
        */
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, vector<uint> width);

        /**
        /param minwidth The minimum size of a box, e.g.\ the widest diameter of a particle.
        /param goalwidth Number of atoms to try and fit per box.
        */
        Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, const flt minwidth, const flt goalwidth);

        /**
        Reshape the grid to optimize number of atom pairs. No-op if `minwidth <= 1` or `goalwidth <=
        0`.
        */
        void optimize_widths();
        /**
        Reshape the grid to optimize number of atom pairs
        */
        void make_grid();

        friend class GridIterator;
        friend class GridPairedIterator;
        /** Iterator over pairs */
        iterator begin();
        iterator end();

        /** Get all pairs of a specific atom */
        pair_iter pairs(atomid a);
        /** Find the amount of time until an atom leaves its current cell. */
        flt time_to_edge(atom& a);
        flt time_to_edge(uint i){return time_to_edge((*atoms)[i]);}

        /** Return a list of all pairs.\ SLOW function, but useful for debugging. */
        vector<idpair> allpairs();
        /** Return a list of all neighbors of a given atom.\ SLOW function, but useful for debugging. */
        vector<atomid> allpairs(atomid a);

        uint numcells(uint i){assert(i<NDIM); return widths[i];}
        #ifdef VEC2D
        uint numcells(){ return widths[0] * widths[1];};
        #else
        uint numcells(){ return widths[0] * widths[1] * widths[2];};
        #endif
};

/** For iterating over all pairs of a single atom */
class GridPairedIterator {
    protected:
        Grid & grid;
        atomid atom1;

        vector<uint> neighbor_cells;
        vector<uint>::iterator cellnum2;
        set<atomid> *cell2;
        set<atomid>::iterator atom2;

        /** returns "successful increment", e.g. cell1 points to an element */
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

/** For iterating over all possible pairs */
class GridIterator {
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
