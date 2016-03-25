#include "trackers.hpp"

void PairList::clear() {
    map<const AtomID, set<AtomID> >::iterator mapiter;
    for (mapiter = pairs.begin(); mapiter != pairs.end(); ++mapiter) {
        mapiter->second.clear();
    }
};

NeighborList::NeighborList(sptr<Box> box, sptr<AtomVec> atomv, const flt skin)
    : box(box),
      skin(skin),
      atoms(atomv),
      diameters(),
      lastlocs(),
      updatenum(0),
      ignorechanged(true){};

bool NeighborList::update_list(bool force) {
    // biggestdist is the distance the furthest-moving Atom has gone
    // bigdist is the next furthest

    if (not force and not ignorechanged) {  // check if we need to update
        flt bigdist = 0, biggestdist = 0;
        for (uint i = 0; i < atoms.size(); i++) {
            Atom &atm = atoms[i];
            flt curdist = (atm.x - lastlocs[i]).norm();
            if (curdist > biggestdist) {
                bigdist = biggestdist;
                biggestdist = curdist;
            } else if (curdist > bigdist) {
                bigdist = curdist;
            } else
                continue;  // this one hasn't moved enough to worry about

            // if we haven't continued, that means bigdist and/or biggestdist
            // have been changed, so we check
            if (bigdist + biggestdist >= skin) {
                force = true;
                break;  // we do need to update, stop checking
            }
        }
        // if we haven't found anything, then we're done; no need to update.
        if (not force) {
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
    for (uint i = 0; i < atoms.size(); i++) {
        AtomID a1 = atoms.get_id(i);
        lastlocs[i] = a1->x;
        for (uint j = 0; j < i; j++) {
            AtomID a2 = atoms.get_id(j);
            if (ignorepairs.has_pair(a1, a2)) continue;
            flt diam = (diameters[i] + diameters[j]) / 2;
            if (box->diff(a1->x, a2->x).norm() < (diam + skin))
                curpairs.push_back(IDPair(a1, a2));
        }
    }
    //~ cout << "NeighborList::update_list:: done.\n";
    // print stuff about the current update
    //~ set<AtomID> curset = (ignorepairs.get_pairs(atoms.back()));
    //~ cout << "NeighborList | atoms: " << atoms.size() <<  "pairs: " <<
    // curpairs.size() << " ignored -1: "
    //~ << curset.size() << "\n";
    //~ cout << "ignored -1:";
    //~ for(set<AtomID>::iterator it=curset.begin(); it!=curset.end(); it++)
    //~ cout << " " << it->n();
    //~ cout << '\n' << "paired:";
    //~ for(vector<IDPair>::iterator it=begin(); it!=end(); it++)
    //~ if(it->first() == atoms.back()) cout << " " << it->last().n();
    //~ cout << '\n';

    return true;
}

Grid::Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, const uint width)
    : box(box), atoms(atoms), minwidth(-1), goalwidth(-1) {
    widths[0] = width;
    widths[1] = width;
#ifndef VEC2D
    widths[2] = width;
#endif
};

Grid::Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, vector<uint> width)
    : box(box), atoms(atoms), minwidth(-1), goalwidth(-1) {
    assert(width.size() == NDIM);
    widths[0] = width[0];
    widths[1] = width[1];
#ifndef VEC2D
    widths[2] = width[2];
#endif
};

Grid::Grid(sptr<OriginBox> box, sptr<AtomGroup> atoms, const flt minwidth,
           const flt goalwidth)
    : box(box), atoms(atoms), minwidth(minwidth), goalwidth(goalwidth) {
    widths[0] = 0;
    widths[1] = 0;
#ifndef VEC2D
    widths[2] = 0;
#endif
    optimize_widths();
};

vector<uint> Grid::neighbors(uint i) {
    if (widths[0] == 1) {
        return vector<uint>(1, 0);
    };
#ifdef VEC2D
    uint yrow = widths[0];

    uint x = i % yrow + widths[0];
    uint y = (i / yrow) % yrow + widths[1];

    vector<uint> v(9);
    v[0] = ((y - 1) % widths[1]) * yrow + ((x - 1) % yrow);
    v[1] = ((y - 1) % widths[1]) * yrow + x % yrow;
    v[2] = ((y - 1) % widths[1]) * yrow + ((x + 1) % yrow);
    v[3] = (y % widths[1]) * yrow + ((x - 1) % yrow);
    v[4] = (y % widths[1]) * yrow + (x % yrow);
    v[5] = (y % widths[1]) * yrow + ((x + 1) % yrow);
    v[6] = ((y + 1) % widths[1]) * yrow + ((x - 1) % yrow);
    v[7] = ((y + 1) % widths[1]) * yrow + x % yrow;
    v[8] = ((y + 1) % widths[1]) * yrow + ((x + 1) % yrow);
    //~ v[0] = (i + fullsize - yrow - xrow) % fullsize;
    //~ v[1] = (i + fullsize - yrow) % fullsize;
    //~ v[2] = (i + fullsize - yrow + xrow) % fullsize;
    //~ v[3] = (i + fullsize - xrow) % fullsize;
    //~ v[4] = (i + fullsize + xrow) % fullsize;
    //~ v[5] = (i + yrow - xrow) % fullsize;
    //~ v[6] = (i + yrow) % fullsize;
    //~ v[7] = (i + yrow + xrow) % fullsize;

    return v;
#else
    uint yrow = widths[0];
    uint zrow = widths[0] * widths[1];

    uint x = (i % yrow + widths[0]);
    uint y = ((i % zrow) / yrow + widths[1]);
    uint z = ((i / zrow) % zrow + widths[2]);

    vector<uint> v(27);
    uint n = 0;
    for (int dz = -1; dz <= 1; dz++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                assert(n < v.size());
                v[n] = (((uint)((int)z + dz) % widths[2]) * zrow +
                        ((uint)((int)y + dy) % widths[1]) * yrow +
                        ((uint)((int)x + dx) % widths[0]));
                n++;
            }
    return v;
#endif
};

void Grid::optimize_widths() {
    if (minwidth <= 0) return;
    flt width_per_atom = pow(box->V() * goalwidth / atoms->size(), OVERNDIM);
    if (width_per_atom < minwidth) width_per_atom = minwidth;
    Vec bshape = box->box_shape();
    widths[0] = (uint)floor(bshape[0] / width_per_atom);
    widths[1] = (uint)floor(bshape[1] / width_per_atom);
#ifdef VEC2D
    if (widths[0] < 3 or widths[1] < 3) {
        widths[0] = 1;
        widths[1] = 1;
    }
#else
    widths[2] = (uint)floor(bshape[2] / width_per_atom);
    if (widths[0] < 3 or widths[1] < 3 or widths[2] < 3) {
        widths[0] = 1;
        widths[1] = 1;
        widths[2] = 1;
    }
#endif
};

uint Grid::get_loc(Vec v, Vec bsize) {
    v = vec_mod(v - bsize / 2., bsize) + bsize / 2;
    uint x = (uint)floor(v[0] * widths[0] / bsize[0]);
    if (x == widths[0]) x = 0;
    uint y = (uint)floor(v[1] * widths[1] / bsize[1]);
    if (y == widths[1]) y = 0;
#ifdef VEC2D
    return y * widths[0] + x;
#else
    uint z = (uint)floor(v[2] * widths[2] / bsize[2]);
    if (z == widths[2]) z = 0;
    return (z * widths[1] + y) * widths[0] + x;
#endif
};

void Grid::make_grid() {
#ifdef VEC2D
    gridlocs = vector<set<AtomID> >(widths[0] * widths[1]);
#else
    gridlocs = vector<set<AtomID> >(widths[0] * widths[1] * widths[2]);
#endif
    Vec bsize = box->box_shape();
    AtomGroup &g = *atoms;
    for (uint ai = 0; ai < g.size(); ai++) {
        uint i = get_loc(g[ai].x, bsize);
        gridlocs[i].insert(g.get_id(ai));
    }
};

Grid::iterator Grid::begin() { return GridIterator(*this); };

Grid::iterator Grid::end() { return GridIterator(*this, gridlocs.end()); };

Grid::pair_iter Grid::pairs(AtomID a) { return pair_iter(*this, a); };

flt Grid::time_to_edge(Atom &a) {
    Vec bsize = box->box_shape();
    Vec v = vec_mod(a.x - bsize / 2., bsize) + bsize / 2;

    flt t = 0;
    for (uint i = 0; i < NDIM; ++i) {
        if (a.v[i] == 0) continue;
        flt cellwidth = bsize[i] / widths[i];
        uint xd = (uint)floor(v[i] / cellwidth + 1.0);
        if (a.v[i] < 0) --xd;
        flt dist = abs(xd * cellwidth - v[i]);
        if (dist < 1e-14) dist += cellwidth;
        flt newt = dist / abs(a.v[i]);
        if (t == 0 or newt < t) t = newt;
    };

    return t;
};

vector<IDPair> Grid::all_pairs() {
    vector<IDPair> v = vector<IDPair>();
    for (iterator p = begin(); p != end(); ++p) {
        IDPair pr = *p;
        cout << "n1: " << pr.first().n() << "  n2: " << pr.last().n() << '\n';
        assert(pr.first().n() < atoms->size());
        assert(pr.last().n() < atoms->size());
        v.push_back(pr);
    }
    return v;
};

vector<AtomID> Grid::all_pairs(AtomID a) {
    vector<AtomID> v = vector<AtomID>();
    for (pair_iter p = pairs(a); p != p.end(); ++p) v.push_back(*p);
    return v;
};

GridPairedIterator::GridPairedIterator(Grid &grid, AtomID a)
    : grid(grid), atom1(a) {
    uint cellnum = grid.get_loc(a->x, grid.box->box_shape());
    neighbor_cells = grid.neighbors(cellnum);
    cellnum2 = neighbor_cells.begin();

    assert(cellnum2 != neighbor_cells.end());
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while (atom2 == cell2->end()) {
        if (!increment_cell2()) return;
        atom2 = cell2->begin();
    };

    while (*atom2 == atom1) {
        if (!increment_atom2()) return;
    }
    if (cellnum2 != neighbor_cells.end()) {
        assert(atom2->n() < grid.atoms->size());
    };
};

bool GridPairedIterator::increment_cell2() {
    if (cellnum2 == neighbor_cells.end()) return false;
    ++cellnum2;
    if (cellnum2 == neighbor_cells.end()) return false;
    cell2 = &(grid.gridlocs[*cellnum2]);
    return true;
};

bool GridPairedIterator::increment_atom2() {
    ++atom2;
    while (atom2 == cell2->end()) {
        if (!increment_cell2()) return false;
        atom2 = cell2->begin();
    };
    if (*atom2 == atom1) return increment_atom2();
    return true;
};

GridPairedIterator &GridPairedIterator::operator++() {
    increment_atom2();
    return *this;
};

bool GridPairedIterator::operator==(const GridPairedIterator &other) {
    if (atom1 != other.atom1) return false;
    if (cellnum2 != other.cellnum2) return false;
    if (cellnum2 == neighbor_cells.end()) return true;  // both at end
    return (atom2 == other.atom2);
}

GridIterator::GridIterator(Grid &grid)
    : grid(grid), cell1(grid.gridlocs.begin()) {
    if (cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    while (atom1 == cell1->end()) {
        if (!increment_cell1()) return;
        atom1 = cell1->begin();
    };
    uint n = (uint)(cell1 - grid.gridlocs.begin());
    neighbor_cells = grid.neighbors(n);
    assert(neighbor_cells.size() > 0);
    cellnum2 = neighbor_cells.begin();
    while (cellnum2 == neighbor_cells.end()) {
        if (!increment_atom1()) return;
        assert(neighbor_cells.size() > 0);
        cellnum2 = neighbor_cells.begin();
    }
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while (atom2 == cell2->end()) {
        if (!increment_cell2()) return;
        atom2 = cell2->begin();
    };

    while (*atom2 == *atom1) {
        if (!increment_atom2()) return;
    }
    if (cell1 != grid.gridlocs.end()) {
        assert(atom1->n() < grid.atoms->size());
        assert(atom2->n() < grid.atoms->size());
    };
};

GridIterator::GridIterator(Grid &grid, vector<set<AtomID> >::iterator cell1)
    : grid(grid), cell1(cell1) {
    if (cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    while (atom1 == cell1->end()) {
        if (!increment_cell1()) return;
        atom1 = cell1->begin();
    };
    uint n = (uint)(cell1 - grid.gridlocs.begin());
    neighbor_cells = grid.neighbors(n);
    assert(neighbor_cells.size() > 0);
    cellnum2 = neighbor_cells.begin();
    while (cellnum2 == neighbor_cells.end()) {
        if (!increment_atom1()) return;
        assert(neighbor_cells.size() > 0);
        cellnum2 = neighbor_cells.begin();
    }
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while (atom2 == cell2->end()) {
        if (!increment_cell2()) return;
        atom2 = cell2->begin();
    };
    while (*atom2 == *atom1) {
        if (!increment_atom2()) return;
    };

    assert(
        (cell1 == grid.gridlocs.end()) or
        (atom1->n() < grid.atoms->size() and atom2->n() < grid.atoms->size()));
};

bool GridIterator::increment_cell1() {
    if (cell1 == grid.gridlocs.end()) return false;
    ++cell1;
    if (cell1 == grid.gridlocs.end()) return false;
    uint n = (uint)(cell1 - grid.gridlocs.begin());
    neighbor_cells = grid.neighbors(n);
    return true;
};

bool GridIterator::increment_atom1() {
    if (cell1 == grid.gridlocs.end()) return false;
    assert(atom1 != cell1->end());
    ++atom1;
    while (atom1 == cell1->end()) {
        if (!increment_cell1()) return false;
        atom1 = cell1->begin();
    };
    return true;
};

bool GridIterator::increment_cell2() {
    if (cell1 == grid.gridlocs.end()) return false;
    assert(cellnum2 != neighbor_cells.end());
    ++cellnum2;
    while (cellnum2 == neighbor_cells.end()) {
        if (!increment_atom1()) return false;
        cellnum2 = neighbor_cells.begin();
    }
    cell2 = &(grid.gridlocs[*cellnum2]);
    return true;
};

bool GridIterator::increment_atom2() {
    if (cell1 == grid.gridlocs.end()) return false;
    assert(atom2 != cell2->end());
    ++atom2;
    while (atom2 == cell2->end()) {
        if (!increment_cell2()) return false;
        atom2 = cell2->begin();
    }
    if (atom2 == atom1) return increment_atom2();
    return true;
};

GridIterator &GridIterator::operator++() {
    increment_atom2();
    assert(
        (cell1 == grid.gridlocs.end()) or
        (atom1->n() < grid.atoms->size() and atom2->n() < grid.atoms->size()));
    return *this;
};

bool GridIterator::operator==(const GridIterator &other) {
    if (cell1 != other.cell1) return false;
    if (cell1 == grid.gridlocs.end()) return true;  // both at end
    if (atom1 != other.atom1) return false;
    if (cellnum2 != other.cellnum2) return false;
    return (atom2 == other.atom2);
};
