#include "box.hpp"

#ifndef TRACKERS_H
#define TRACKERS_H

#include <map>
#include <set>
#include "assert.h"

/**
The general interface for a "tracker", a class that needs to be called every
timestep.

StateTracker is used as a base for some "helper" classes, like neighbor lists
(\ref NeighborList),
as well as for many statistics.
*/
class StateTracker {
   public:
    /** This function is called once per timestep, when particles are in their
     * set position.*/
    virtual void update(Box &box) = 0;
    /** for CollectionCD (and heirs). Returns whether or not to call this
    function after every collision, regardless of timestep.
    Most trackers do not need / want this, so the default is false.
    */
    virtual bool every_collision() { return false; };
    virtual void update_collision(Box &box, AtomID a1, AtomID a2, flt time,
                                  Vec delta_p){};
    virtual ~StateTracker(){};
};

/** A mapping of \ref Atom -> [list of \ref Atom], used by \ref NeighborList to
keep track
of what atoms are near what other atoms.

Pairs are assumed to be symmetric, i.e. either `(a1, a2)` is stored or `(a2,
a1)` is stored, and
only one will be returned through the iterator.
*/
class PairList {
   protected:
    map<const AtomID, set<AtomID> > pairs;

   public:
    //~ PairList(AtomGroup *group);
    PairList(){};

    /** Ensure that a given \ref Atom is in the overall list. */
    inline void ensure(const AtomID a) {
        pairs.insert(std::pair<AtomID, set<AtomID> >(a, set<AtomID>()));
    }
    /** Ensure that every \ref Atom in input is in the overall list. */
    inline void ensure(vector<AtomID> ps) {
        vector<AtomID>::iterator it;
        for (it = ps.begin(); it != ps.end(); ++it) ensure(*it);
    }
    /** Ensure that every \ref Atom in input is in the overall list. */
    inline void ensure(AtomGroup &group) {
        for (uint i = 0; i < group.size(); i++) ensure(group.get_id(i));
    }
    /** Check if a given pair is is the list.

    Note that `has_pair(a1, a2) == has_pair(a2, a1)`.
    */
    inline bool has_pair(AtomID a1, AtomID a2) {
        if (a1 > a2)
            return pairs[a1].count(a2) > 0;
        else
            return pairs[a2].count(a1) > 0;
    }

    /**
    Add a pair to the list. If already in the list, nothing changes.
    */
    inline void add_pair(AtomID a1, AtomID a2) {
        //~ cout << "PairList ignore " << a1.n() << '-' << a2.n() << "\n";
        if (a1 > a2) {
            pairs[a1].insert(a2);
        } else {
            pairs[a2].insert(a1);
        };
    }

    /**
    Remove a pair from the list. If not already in the list, nothing changes.
    */
    inline void erase_pair(AtomID a1, AtomID a2) {
        if (a1 > a2) {
            pairs[a1].erase(a2);
        } else {
            pairs[a2].erase(a1);
        };
    }

    /**
    Get all atoms neighboring an Atom *that are "lesser"* than the Atom.

    Because of the way neighbors are stored, there is no easy way to actually
    get *all*
    neighbors of an Atom, but if `get_pairs` is called successively on all
    atoms, it will return
    each and every pair exactly once (including "false positives").
    */
    inline set<AtomID> get_pairs(const AtomID a) {
        ensure(a);
        return pairs[a];
    };

    /** for iterating over neighbors */
    inline set<AtomID>::iterator begin(const AtomID a) {
        return pairs[a].begin();
    };
    /** for iterating over neighbors */
    inline set<AtomID>::iterator end(const AtomID a) { return pairs[a].end(); };

    /** for iterating over neighbors */
    inline uint size() const {
        uint N = 0;
        for (map<const AtomID, set<AtomID> >::const_iterator it = pairs.begin();
             it != pairs.end(); ++it)
            N += (uint)it->second.size();
        return N;
    };

    /** Clear the *inner* lists. */
    void clear();
};

/**
Maintains a Verlet list of "neighbors": molecules within a 'skin radius' of each
other. When the
list is updated (`update_list(true)`), molecules are counted as neighbors if the
distance between
them is \f$d_{ij} < r_i + r_j + R\f$, where \f$r_i\f$ and \f$r_j\f$ are the
radii, and \f$R\f$ is
the "skin radius". On "normal" updates (`update(box)` or `update_list(false)`,
All atoms fulfilling
the condition \f$d_{ij} < r_i + r_j\f$ are guaranteed to be in the list, but
there may be some false
positives.

`update_list(false)` should be called frequently (i.e. every timestep); it
checks (\f$O(N)\f$) if
any two molecules might conceivably overlap by more than a critical distance,
and if so, it updates
all the neighbor lists.

The NeighborList has an inherent tradeoff, set by the skin radius: If the skin
is large, then
full updates (which are slow) will be needed less often, but more particles will
be considered
"neighbors" even when they are not within \f$d_{ij} < r_i + r_j\f$, i.e. there
will be more false
positives.
*/
class NeighborList : public StateTracker {
   protected:
    sptr<Box> box;
    flt skin;
    SubGroup atoms;
    vector<flt> diameters;
    vector<IDPair> curpairs;
    PairList ignorepairs;
    vector<Vec> lastlocs;
    uint updatenum;
    bool ignorechanged; /*< if true, forces a full check on next update */
   public:
    NeighborList(sptr<Box> box, sptr<AtomVec> atoms, const flt skin);
    /**
    Calls `update_list(false)`; exists as a requirement of \ref StateTracker.
    */
    void update(Box &newbox) {
        assert(&newbox == box.get());
        update_list(false);
    };
    /**
    Update the neighbor list based on current positions. If `force`, the list is
    updated
    immediately \f$O(N^2)\f$; otherwise, NeighborList checks if an update is
    necessary
    \f$O(N)\f$, and updates the list if so.
    */
    bool update_list(bool force = true);
    // if force = false, we check if updating necessary first

    AtomVec &vec() { return atoms.vec(); };
    inline uint which() { return updatenum; };
    inline uint numpairs() { return (uint)curpairs.size(); };
    inline void ignore(AtomID a, AtomID b) {
        ignorepairs.add_pair(a, b);
        ignorechanged = true;
    };
    void add(AtomID a, flt diameter) {
        atoms.add(a);
        // assert(diameters.size() == atoms.size() - 1);
        diameters.push_back(diameter);
        // assert(lastlocs.size() == atoms.size() - 1);
        lastlocs.push_back(a->x);
        ignorechanged = true;
    }

    /** Number of pairs currently being ignored */
    inline uint ignore_size() const { return ignorepairs.size(); };
    /** Number of atoms in list */
    inline uint size() const { return atoms.size(); };
    inline vector<IDPair>::iterator begin() { return curpairs.begin(); };
    inline vector<IDPair>::iterator end() { return curpairs.end(); };
    inline IDPair get(uint i) {
        // assert(i<curpairs.size());
        return curpairs[i];
    };
    ~NeighborList(){};
};

class GridIterator;
class GridPairedIterator;

/**
A fast algorithm for finding all pairs of neighboring atoms.

A Grid splits the OriginBox into "cells", where no cell is narrower than any
Atom. For any given Atom,
all atoms within the same cell or any neighboring cell are considered
"neighbors".
*/
class Grid {
   public:
    sptr<OriginBox> box;
    sptr<AtomGroup> atoms;

    /**
    minwidth is minimum size of a box, in real units.

    This is used to calculate the widths of new boxes in `optimize_widths()`.

    May be 0 or negative, to indicate that `optimize_widths()` should never
    change the box sizes.
    */
    flt minwidth;
    flt goalwidth;     /**< goalwidth is how many atoms per box (goal). See
                          minwidth. */
    uint widths[NDIM]; /**< Number of box divisions per dimension. */
    vector<set<AtomID> > gridlocs; /**< Atoms in each grid location. */

    vector<uint> neighbors(uint i);
    uint get_loc(Vec v, Vec bsize);

   public:
    typedef GridIterator iterator;
    typedef GridPairedIterator pair_iter;

    /**
    /param width Number of divisions along all axes
    */
    Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, const uint width = 1);

    /**
    /param widths Number of divisions along each axis
    */
    Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, vector<uint> width);

    /**
    /param minwidth The minimum size of a box, e.g.\ the widest diameter of a
    particle.
    /param goalwidth Number of atoms to try and fit per box.
    */
    Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt minwidth,
         const flt goalwidth);

    /**
    Reshape the grid to optimize number of Atom pairs. No-op if `minwidth <= 1`
    or `goalwidth <=
    0`.
    */
    void optimize_widths();
    /**
    Reshape the grid to optimize number of Atom pairs
    */
    void make_grid();

    friend class GridIterator;
    friend class GridPairedIterator;
    /** Iterator over pairs */
    iterator begin();
    iterator end();

    /** Get all pairs of a specific Atom */
    pair_iter pairs(AtomID a);
    /** Find the amount of time until an Atom leaves its current cell. */
    flt time_to_edge(Atom &a);
    flt time_to_edge(uint i) { return time_to_edge((*atoms)[i]); }

    /** Return a list of all pairs.\ SLOW function, but useful for debugging. */
    vector<IDPair> all_pairs();
    /** Return a list of all neighbors of a given Atom.\ SLOW function, but
     * useful for debugging. */
    vector<AtomID> all_pairs(AtomID a);

    uint numcells(uint i) {
        assert(i < NDIM);
        return widths[i];
    }
#ifdef VEC2D
    uint numcells() { return widths[0] * widths[1]; };
#else
    uint numcells() { return widths[0] * widths[1] * widths[2]; };
#endif
};

/** For iterating over all pairs of a single Atom */
class GridPairedIterator {
   protected:
    Grid &grid;
    AtomID atom1;

    vector<uint> neighbor_cells;
    vector<uint>::iterator cellnum2;
    set<AtomID> *cell2;
    set<AtomID>::iterator atom2;

    /** returns "successful increment", e.g. cell1 points to an element */
    bool increment_cell2();
    bool increment_atom2();

   public:
    typedef set<AtomID>::iterator end_type;

    GridPairedIterator(Grid &grid, AtomID a);

    GridPairedIterator &operator++();
    AtomID operator*() { return *atom2; };
    bool operator==(const GridPairedIterator &other);

    bool operator!=(const GridPairedIterator &other) {
        return !(*this == other);
    };

    bool operator==(const end_type other) { return (atom2 == other); }
    bool operator!=(const end_type other) { return !(atom2 == other); };

    end_type end() { return grid.gridlocs[neighbor_cells.back()].end(); };
};

/** For iterating over all possible pairs */
class GridIterator {
   protected:
    Grid &grid;
    vector<set<AtomID> >::iterator cell1;
    set<AtomID>::iterator atom1;

    vector<uint> neighbor_cells;
    vector<uint>::iterator cellnum2;
    set<AtomID> *cell2;
    set<AtomID>::iterator atom2;

    // returns "successful increment", e.g. cell1 points to an element
    bool increment_cell1();
    bool increment_atom1();
    bool increment_cell2();
    bool increment_atom2();

   public:
    GridIterator(Grid &grid);
    GridIterator(Grid &grid, vector<set<AtomID> >::iterator cell1);

    GridIterator &operator++();
    IDPair operator*() { return IDPair(*atom1, *atom2); };
    bool operator==(const GridIterator &other);

    bool operator!=(const GridIterator &other) { return !(*this == other); };
};

#endif
