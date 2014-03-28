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

vector<uint> Grid::neighbors(uint i){
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

uint get_width(flt L, flt minwidth, flt goalwidth){
    uint n = (uint)roundflt(L / goalwidth);
    uint nmin = (uint)ceilflt(L / minwidth);
    return n > nmin ? n : nmin;
}

void Grid::update_widths(){
    if(goalwidth <= 0) return;
    Vec bsize = box->boxshape();
    widths[0] = get_width(bsize[0], minwidth, goalwidth);
    widths[1] = get_width(bsize[1], minwidth, goalwidth);
    #ifdef VEC2D
    if(widths[0] < 3 or widths[1] < 3){
        widths[0] = 1;
        widths[1] = 1;
    }
    #else
    widths[2] = get_width(bsize[2], minwidth, goalwidth);
    
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
    uint y = (uint)floorflt(v[1] * widths[1] / bsize[1]);
    #ifdef VEC2D
    return y * widths[0] + x;
    #else
    uint z = (uint)floorflt(v[2] * widths[2] / bsize[2]);
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
    for(uint ai = 0; ai < atoms.size(); ai++){
        uint i = get_loc(atoms[ai].x, bsize);
        gridlocs[i].insert(atoms.get_id(ai));
    }
};

Grid::iterator Grid::begin(){
    return GridIterator(*this);
};

Grid::iterator Grid::end(){
    return GridIterator(*this, gridlocs.end());
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
    atom1++;
    while(atom1 == cell1->end()){
        if (!increment_cell1()) return false;
        atom1 = cell1->begin();
    };
    return true;
};

bool GridIterator::increment_cell2(){
    cellnum2++;
    while(cellnum2 == neighbor_cells.end()){
        if (!increment_atom1()) return false;
        cellnum2 = neighbor_cells.begin();
    }
    cell2 = &(grid.gridlocs[*cellnum2]);
    return true;
};

bool GridIterator::increment_atom2(){
    atom2++;
    while(atom2 == cell2->end() or atom2 == atom1){
        if (!increment_cell2()) return false;
        atom2 = cell2->begin();
    }
    if(atom2 == atom1) return increment_atom2();
    return true;
};

GridIterator::GridIterator(Grid & grid) : grid(grid), cell1(grid.gridlocs.begin()){
    if(cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    uint n = cell1 - grid.gridlocs.begin();
    neighbor_cells = grid.neighbors(n);
    cellnum2 = neighbor_cells.begin();
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
};

GridIterator::GridIterator(Grid & grid, vector<set<atomid> >::iterator cell1) : 
                        grid(grid), cell1(cell1){
    if(cell1 == grid.gridlocs.end()) return;
    atom1 = cell1->begin();
    uint n = cell1 - grid.gridlocs.begin();
    neighbor_cells = grid.neighbors(n);
    cellnum2 = neighbor_cells.begin();
    // neighbor_cells should be non-empty
    cell2 = &(grid.gridlocs[*cellnum2]);
    atom2 = cell2->begin();
};

GridIterator& GridIterator::operator++(){
    increment_atom2();
    return *this;
};

bool GridIterator::operator==(const GridIterator &other){
    if(cell1 != other.cell1) return false;
    if(cell1 == grid.gridlocs.end()) return true; // both at end
    if(atom1 != other.atom1) return false;
    if(cellnum2 != other.cellnum2) return false;
    if(atom2 != other.atom2) return false;
    
    return true;
};
