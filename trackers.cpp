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

Grid::Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, const uint width)
        : box(box), atoms(atoms), minwidth(-1), goalwidth(-1){
    widths[0] = width;
    widths[1] = width;
    #ifndef VEC2D
    widths[2] = width;
    #endif
};

Grid::Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, vector<uint> width)
        : box(box), atoms(atoms), minwidth(-1), goalwidth(-1){
    assert(width.size() == NDIM);
    widths[0] = width[0];
    widths[1] = width[1];
    #ifndef VEC2D
    widths[2] = width[2];
    #endif
};

Grid::Grid(sptr<OriginBox> box, sptr<atomgroup> atoms, 
                            const flt goalwidth, const flt minwidth) :
        box(box), atoms(atoms), minwidth(minwidth), goalwidth(goalwidth){
    widths[0] = 0;
    widths[1] = 0;
    #ifndef VEC2D
    widths[2] = 0;
    #endif
    optimize_widths();
};

vector<uint> Grid::neighbors(uint i){
    if(widths[0] == 1){
        return vector<uint>(1,0);
    };
    #ifdef VEC2D
    uint yrow = widths[0];
    
    uint x = i % yrow + widths[0];
    uint y = (i / yrow) % yrow + widths[1];
    
    vector<uint> v(9);
    v[0] = ((y - 1) % widths[1]) * yrow + ((x - 1) % yrow);
    v[1] = ((y - 1) % widths[1]) * yrow + x        % yrow;
    v[2] = ((y - 1) % widths[1]) * yrow + ((x + 1) % yrow);
    v[3] = (y       % widths[1]) * yrow + ((x - 1) % yrow);
    v[4] = (y       % widths[1]) * yrow + (x       % yrow);
    v[5] = (y       % widths[1]) * yrow + ((x + 1) % yrow);
    v[6] = ((y + 1) % widths[1]) * yrow + ((x - 1) % yrow);
    v[7] = ((y + 1) % widths[1]) * yrow + x        % yrow;
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
    int yrow = widths[0];
    int zrow = widths[0]*widths[1];
    
    int x = i % yrow + widths[0];
    int y = (i % zrow) / yrow + widths[1];
    int z = (i / zrow) % zrow + widths[2];
    
    vector<uint> v(27);
    uint n=0;
    for(int dz=-1; dz<=1; dz++)
    for(int dy=-1; dy<=1; dy++)
    for(int dx=-1; dx<=1; dx++){
        assert(n < v.size());
        v[n] = (
            ((z + dz) % widths[2]) * zrow + 
            ((y + dy) % widths[1]) * yrow + 
            ((x + dx) % widths[0])
            );
        n++;
    }
    return v;
    #endif
};

void Grid::optimize_widths(){
    if(minwidth <= 0) return;
    flt width_per_atom = pow(box->V() * goalwidth / atoms->size(), OVERNDIM);
    if (width_per_atom < minwidth) width_per_atom = minwidth;
    Vec bshape = box->boxshape();
    widths[0] = (uint) floorflt(bshape[0] / width_per_atom);
    widths[1] = (uint) floorflt(bshape[1] / width_per_atom);
    #ifdef VEC2D
    if(widths[0] < 3 or widths[1] < 3){
        widths[0] = 1;
        widths[1] = 1;
    }
    #else
    widths[2] = (uint) floorflt(bshape[2] / width_per_atom);
    if(widths[0] < 3 or widths[1] < 3 or widths[2] < 3){
        widths[0] = 1;
        widths[1] = 1;
        widths[2] = 1;
    }
    #endif
};

uint Grid::get_loc(Vec v, Vec bsize){
    v = vecmod(v - bsize/2., bsize) + bsize/2;
    uint x = (uint)floorflt(v[0] * widths[0] / bsize[0]);
    if(x == widths[0]) x = 0;
    uint y = (uint)floorflt(v[1] * widths[1] / bsize[1]);
    if(y == widths[1]) y = 0;
    #ifdef VEC2D
    return y * widths[0] + x;
    #else
    uint z = (uint)floorflt(v[2] * widths[2] / bsize[2]);
    if(z == widths[2]) z = 0;
    return (z * widths[1] + y) * widths[0] + x;
    #endif
};

void Grid::make_grid(){
    #ifdef VEC2D
    gridlocs = vector<set<atomid> >(widths[0] * widths[1]);
    #else
    gridlocs = vector<set<atomid> >(widths[0] * widths[1] * widths[2]);
    #endif
    Vec bsize = box->boxshape();
    atomgroup &g = *atoms;
    for(uint ai = 0; ai < g.size(); ai++){
        uint i = get_loc(g[ai].x, bsize);
        gridlocs[i].insert(g.get_id(ai));
    }
};

Grid::iterator Grid::begin(){
    return GridIterator(*this);
};

Grid::iterator Grid::end(){
    return GridIterator(*this, gridlocs.end());
};

Grid::pair_iter Grid::pairs(atomid a){
    return pair_iter(*this, a);
};

flt Grid::time_to_edge(atom &a){
    Vec bsize = box->boxshape();
    Vec v = vecmod(a.x - bsize/2., bsize) + bsize/2;
    
    flt t = 0;
    for(uint i=0; i<NDIM; ++i){
        if(a.v[i] == 0) continue;
        flt cellwidth = bsize[i] / widths[i];
        uint xd = (uint)floorflt(v[i] / cellwidth + 1.0);
        if(a.v[i] < 0) --xd;
        flt dist = abs(xd*cellwidth - v[i]);
        if(dist < 1e-14) dist += cellwidth;
        flt newt = dist / abs(a.v[i]);
        if(t == 0 or newt < t) t = newt;
    };
    
    return t;
};

vector<idpair> Grid::allpairs(){
    vector<idpair> v = vector<idpair>();
    for(iterator p=begin(); p!=end(); ++p){
        idpair pr = *p;
        cout << "n1: " << pr.first().n()
                << "  n2: "  << pr.last().n() << '\n';
        assert(pr.first().n() < atoms->size());
        assert(pr.last().n() < atoms->size());
        v.push_back(pr);
    }
    return v;
};

vector<atomid> Grid::allpairs(atomid a){
    vector<atomid> v = vector<atomid>();
    for(pair_iter p=pairs(a); p!=p.end(); ++p) v.push_back(*p);
    return v;
};

GridPairedIterator::GridPairedIterator(Grid & grid, atomid a) : 
        grid(grid), atom1(a){
    uint cellnum = grid.get_loc(a.x(), grid.box->boxshape());
    neighbor_cells = grid.neighbors(cellnum);
    cellnum2 = neighbor_cells.begin();
    
    assert(cellnum2 != neighbor_cells.end());
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while(atom2 == cell2->end()){
        if(!increment_cell2()) return;
        atom2 = cell2->begin();
    };
    
    while(*atom2 == atom1){
        if(!increment_atom2()) return;
    }
    if(cellnum2 != neighbor_cells.end()){
        assert(atom2->n() < grid.atoms->size());
    };
};

bool GridPairedIterator::increment_cell2(){
    if(cellnum2 == neighbor_cells.end()) return false;
    cellnum2++;
    if(cellnum2 == neighbor_cells.end()) return false;
    cell2 = &(grid.gridlocs[*cellnum2]);
    return true;
};

bool GridPairedIterator::increment_atom2(){
    atom2++;
    while(atom2 == cell2->end()){
        if (!increment_cell2()) return false;
        atom2 = cell2->begin();
    };
    if(*atom2 == atom1) return increment_atom2();
    return true;
};

GridPairedIterator& GridPairedIterator::operator++(){
    increment_atom2();
    return *this;
};

bool GridPairedIterator::operator==(const GridPairedIterator &other){
    if(atom1 != other.atom1) return false;
    if(cellnum2 != other.cellnum2) return false;
    if(cellnum2 == neighbor_cells.end()) return true; // both at end
    return (atom2 == other.atom2);
}

GridIterator::GridIterator(Grid & grid) : grid(grid), cell1(grid.gridlocs.begin()){
    if(cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    while(atom1 == cell1->end()){
        if(!increment_cell1()) return;
        atom1 = cell1->begin();
    };
    uint n = cell1 - grid.gridlocs.begin();
    neighbor_cells = grid.neighbors(n);
    assert(neighbor_cells.size() > 0);
    cellnum2 = neighbor_cells.begin();
    while(cellnum2 == neighbor_cells.end()){
        if (!increment_atom1()) return;
        assert(neighbor_cells.size() > 0);
        cellnum2 = neighbor_cells.begin();
    }
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while(atom2 == cell2->end()){
        if(!increment_cell2()) return;
        atom2 = cell2->begin();
    };
    
    while(*atom2 == *atom1){
        if(!increment_atom2()) return;
    }
    if(cell1 != grid.gridlocs.end()){
        assert(atom1->n() < grid.atoms->size());
        assert(atom2->n() < grid.atoms->size());
    };
};

GridIterator::GridIterator(Grid & grid, vector<set<atomid> >::iterator cell1) : 
                        grid(grid), cell1(cell1){
    if(cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    while(atom1 == cell1->end()){
        if(!increment_cell1()) return;
        atom1 = cell1->begin();
    };
    uint n = cell1 - grid.gridlocs.begin();
    neighbor_cells = grid.neighbors(n);
    assert(neighbor_cells.size() > 0);
    cellnum2 = neighbor_cells.begin();
    while(cellnum2 == neighbor_cells.end()){
        if (!increment_atom1()) return;
        assert(neighbor_cells.size() > 0);
        cellnum2 = neighbor_cells.begin();
    }
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
    while(atom2 == cell2->end()){
        if(!increment_cell2()) return;
        atom2 = cell2->begin();
    };
    while(*atom2 == *atom1){
        if(!increment_atom2()) return;
    };
    
    assert((cell1 == grid.gridlocs.end()) or 
        (atom1->n() < grid.atoms->size() and atom2->n() < grid.atoms->size()));
};

bool GridIterator::increment_cell1(){
    if(cell1 == grid.gridlocs.end()) return false;
    cell1++;
    if(cell1 == grid.gridlocs.end()) return false;
    uint n = cell1 - grid.gridlocs.begin();
    neighbor_cells = grid.neighbors(n);
    return true;
};

bool GridIterator::increment_atom1(){
    if(cell1 == grid.gridlocs.end()) return false;
    assert(atom1 != cell1->end());
    atom1++;
    while(atom1 == cell1->end()){
        if (!increment_cell1()) return false;
        atom1 = cell1->begin();
    };
    return true;
};

bool GridIterator::increment_cell2(){
    if(cell1 == grid.gridlocs.end()) return false;
    assert(cellnum2 != neighbor_cells.end());
    cellnum2++;
    while(cellnum2 == neighbor_cells.end()){
        if (!increment_atom1()) return false;
        cellnum2 = neighbor_cells.begin();
    }
    cell2 = &(grid.gridlocs[*cellnum2]);
    return true;
};

bool GridIterator::increment_atom2(){
    if(cell1 == grid.gridlocs.end()) return false;
    assert(atom2 != cell2->end());
    atom2++;
    while(atom2 == cell2->end()){
        if (!increment_cell2()) return false;
        atom2 = cell2->begin();
    }
    if(atom2 == atom1) return increment_atom2();
    return true;
};

GridIterator& GridIterator::operator++(){
    increment_atom2();
    assert((cell1 == grid.gridlocs.end()) or 
        (atom1->n() < grid.atoms->size() and atom2->n() < grid.atoms->size()));
    return *this;
};

bool GridIterator::operator==(const GridIterator &other){
    if(cell1 != other.cell1) return false;
    if(cell1 == grid.gridlocs.end()) return true; // both at end
    if(atom1 != other.atom1) return false;
    if(cellnum2 != other.cellnum2) return false;
    return (atom2 == other.atom2);
};

