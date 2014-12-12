#include "constraints.hpp"

ContactTracker::ContactTracker(sptr<Box> box, sptr<atomgroup> atoms, vector<flt> dists) :
    atoms(atoms), dists(dists), contacts(), breaks(0), formations(0),
        incontact(0){
    //~ cout << "Making contact tracker." << endl;
    uint N = atoms->size();
    dists.resize(N);
    contacts.resize(N);
    for(uint i=0; i<N; ++i){
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
    for(uint i=0; i<N; ++i){
        Vec ri = atoms->get(i).x;
        for(uint j=0; j<i; ++j){
            Vec rj = atoms->get(j).x;
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
    for(uint i=0; i<Natoms; ++i){
        atom &curatom = atoms->get(i);
        curK += curatom.v.sq() * curatom.m / 2;
    }
    
    vector<sptr<interaction> >::iterator it;
    for(it = interactions.begin(); it != interactions.end(); ++it){
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
    for(it = interactions.begin(); it != interactions.end(); ++it){
        curU += (*it)->energy(box);
    }
    setU0(curU);
};


RsqTracker1::RsqTracker1(atomgroup& atoms, unsigned long skip, Vec com) :
            pastlocs(atoms.size(), Vec()), xyz2sums(atoms.size(), Vec()),
            xyz4sums(atoms.size(), Vec()), r4sums(atoms.size(), 0),
            skip(skip), count(0){
    for(uint i = 0; i<atoms.size(); ++i){
        pastlocs[i] = atoms[i].x - com;
    };
};
        
void RsqTracker1::reset(atomgroup& atoms, Vec com){
    pastlocs.resize(atoms.size());
    xyz2sums.assign(atoms.size(), Vec());
    xyz4sums.assign(atoms.size(), Vec());
    r4sums.assign(atoms.size(), 0);
    for(uint i = 0; i<atoms.size(); ++i){
        pastlocs[i] = atoms[i].x - com;
    };
    count = 0;
};

bool RsqTracker1::update(Box& box, atomgroup& atoms, unsigned long t, Vec com){
    if(t % skip != 0) return false;
    
    for(uint i = 0; i<atoms.size(); ++i){
        //flt dist = box.diff(atoms[i].x, pastlocs[i]).sq();
        Vec r = atoms[i].x - com;
        // We don't want the boxed distance - we want the actual distance moved!
        Vec distsq = r - pastlocs[i];
        Vec distsqsq = Vec();
        flt dist4 = 0;
        for(uint j=0; j<NDIM; ++j){
            distsq[j] *= distsq[j];
            distsqsq[j] = distsq[j]*distsq[j];
            dist4 += distsq[j];
        };
        
        dist4 *= dist4;
        
        xyz2sums[i] += distsq;
        xyz4sums[i] += distsqsq;
        r4sums[i] += dist4;
        pastlocs[i] = r;
    };
    count += 1;
    return true;
};

vector<Vec> RsqTracker1::xyz2(){
    vector<Vec> means(xyz2sums.size(), Vec());
    for(uint i=0; i<xyz2sums.size(); ++i){
        means[i] = xyz2sums[i] / ((flt) count);
    }
    return means;
};

vector<Vec> RsqTracker1::xyz4(){
    vector<Vec> means(xyz4sums.size(), Vec());
    for(uint i=0; i<xyz4sums.size(); ++i){
        means[i] = xyz4sums[i] / ((flt) count);
    }
    return means;
};

vector<flt> RsqTracker1::r4(){
    vector<flt> means(r4sums.size(), 0);
    for(uint i=0; i<r4sums.size(); ++i){
        means[i] = r4sums[i] / ((flt) count);
    }
    return means;
};
    

RsqTracker::RsqTracker(sptr<atomgroup> atoms, vector<unsigned long> ns, bool usecom) : 
            atoms(atoms), curt(0), usecom(usecom){
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<unsigned long>::iterator n=ns.begin(); n!=ns.end(); ++n){
        singles.push_back(RsqTracker1(*atoms, *n, com));
    }
};

void RsqTracker::update(Box &box){
    curt++;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        it->update(box, *atoms, curt, com);
    }
};

void RsqTracker::reset(){
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        it->reset(*atoms, com);
    }
};

vector<vector<Vec> > RsqTracker::xyz2(){
    vector<vector<Vec> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back(it->xyz2());
    }
    return vals;
};

vector<vector<flt> > RsqTracker::r2(){
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vector<flt> val;
        vector<Vec> xyz2 = it->xyz2();
        val.reserve(xyz2.size());
        
        for(vector<Vec>::iterator it2=xyz2.begin(); it2!=xyz2.end(); ++it2){
            Vec v = *it2;
            flt r2 = 0;
            for(uint i=0; i<NDIM; ++i){
                r2 += v[i];
            }
            val.push_back(r2);
        }
        vals.push_back(val);
    }
    return vals;
};

vector<vector<Vec> > RsqTracker::xyz4(){
    vector<vector<Vec> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back(it->xyz4());
    }
    return vals;
};

vector<vector<flt> > RsqTracker::r4(){
    vector<vector<flt> > vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back(it->r4());
    }
    return vals;
};
vector<flt> RsqTracker::counts(){
    vector<flt> vals;
    vals.reserve(singles.size());
    for(vector<RsqTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back((flt) it->get_count());
    }
    return vals;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
ISFTracker1::ISFTracker1(atomgroup& atoms, unsigned long skip, vector<flt> ks, Vec com) :
            pastlocs(atoms.size(), Vec()), ISFsums(ks.size(), vector<Array<cmplx, NDIM> >()),
            ks(ks), skip(skip), count(0){
    for(uint i = 0; i<atoms.size(); ++i){
        pastlocs[i] = atoms[i].x - com;
    };
    
    for(uint ki=0; ki<ks.size(); ++ki){
        ISFsums[ki] = vector<Array<cmplx, NDIM> >(atoms.size(), Array<cmplx, NDIM>());
    };
};
        
void ISFTracker1::reset(atomgroup& atoms, Vec com){
    pastlocs.resize(atoms.size());
    for(uint ki=0; ki<ks.size(); ki++){
        ISFsums[ki].assign(atoms.size(), Array<cmplx, NDIM>());
    }
    for(uint i = 0; i<atoms.size(); ++i){
        pastlocs[i] = atoms[i].x - com;
    };
    count = 0;
};

bool ISFTracker1::update(Box& box, atomgroup& atoms, unsigned long t, Vec com){
    if(t % skip != 0) return false;
    
    for(uint i = 0; i<atoms.size(); ++i){
        //flt dist = box.diff(atoms[i].x, pastlocs[i]).sq();
        Vec r = atoms[i].x - com;
        // We don't want the boxed distance - we want the actual distance moved!
        Vec dr = r - pastlocs[i];
        for(uint ki=0; ki<ks.size(); ++ki){
            for(uint j=0; j<NDIM; ++j){
                //~ cout << "i: " << i << "  ki: " << ki << "  j: " << j;
                //~ cout << "  dr[j]: " << dr[j];
                //~ cout << "  ks[ki]: " << ks[ki];
                //~ cout << endl; 
                
                ISFsums[ki][i][j] += exp(cmplx(0, ks[ki]*dr[j]));
            }
        }
        
        pastlocs[i] = r;
    };
    count += 1;
    return true;
};

vector<vector<cmplx> > ISFTracker1::ISFs() {
    vector<vector<cmplx> > means(ks.size(), vector<cmplx>());
    for(uint ki=0; ki<ks.size(); ++ki){
        means[ki].assign(ISFsums[ki].size(), 0);
        for(uint i = 0; i<ISFsums[ki].size(); ++i){
            for(uint j=0; j<NDIM; ++j){
                means[ki][i] += ISFsums[ki][i][j];
            }
            means[ki][i] /= NDIM;
            means[ki][i] /= (flt) count;
        }
    }
    
    return means;
};

vector<vector<Array<cmplx, NDIM> > > ISFTracker1::ISFxyz() {
    vector<vector<Array<cmplx, NDIM> > > means(ks.size(), vector<Array<cmplx, NDIM> >());
    for(uint ki=0; ki<ks.size(); ++ki){
        means[ki].assign(ISFsums[ki].size(), Array<cmplx, NDIM>());
        for(uint i = 0; i<ISFsums[ki].size(); ++i){
            for(uint j = 0; j<NDIM; ++j)
                means[ki][i][j] = ISFsums[ki][i][j] / (cmplx((flt) count, 0));
        }
    }
    
    return means;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

ISFTracker::ISFTracker(sptr<atomgroup> atoms, vector<flt> ks, 
                    vector<unsigned long> ns, bool usecom) : atoms(atoms), curt(0), usecom(usecom){
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<unsigned long>::iterator n=ns.begin(); n!=ns.end(); ++n){
        singles.push_back(ISFTracker1(*atoms, *n, ks, com));
    }
};

void ISFTracker::update(Box &box){
    curt++;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<ISFTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        it->update(box, *atoms, curt, com);
    }
};

void ISFTracker::reset(){
    curt = 0;
    Vec com = usecom ? atoms->com() : Vec();
    for(vector<ISFTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        it->reset(*atoms, com);
    }
};

vector<vector<vector<cmplx> > > ISFTracker::ISFs(){
    vector<vector<vector<cmplx> > > vals;
    vals.reserve(singles.size());
    for(vector<ISFTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back(it->ISFs());
    }
    return vals;
};

vector<vector<vector<Array<cmplx, NDIM> > > > ISFTracker::ISFxyz(){
    vector<vector<vector<Array<cmplx, NDIM> > > > vals;
    vals.reserve(singles.size());
    for(vector<ISFTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back(it->ISFxyz());
    }
    return vals;
};

vector<flt> ISFTracker::counts(){
    vector<flt> vals;
    vals.reserve(singles.size());
    for(vector<ISFTracker1>::iterator it=singles.begin(); it!=singles.end(); ++it){
        vals.push_back((flt) it->get_count());
    }
    return vals;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

SmoothLocs::SmoothLocs(sptr<atomgroup> atoms, Box &box, uint smoothn, uint skipn, bool usecom) : 
        atoms(atoms), smoothn(smoothn), skipn(skipn), curlocs(atoms->size()), numincur(0), locs(),
        curt(0), usecom(usecom){
    update(box);
};


void SmoothLocs::reset(){
    curlocs.assign(atoms->size(), Vec());
    numincur = 0;
    curt = 0;
    locs.clear();
};

void SmoothLocs::update(Box &box){
    if(curt % skipn != 0){
        curt++;
        return;
    };
    
    numincur++;
    bool smooth_time = (numincur >= smoothn);
    
    Vec com = usecom ? atoms->com() : Vec();
    for(uint i = 0; i<atoms->size(); ++i){
        Vec curloc = (*atoms)[i].x - com;
        curlocs[i] += curloc;
        if(smooth_time) {
            curlocs[i] /= numincur;
        }
    };
    
    if(smooth_time){
        locs.push_back(curlocs);
        curlocs.assign(atoms->size(), Vec());
        numincur = 0;
    };
    
    curt++;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
RDiffs::RDiffs(sptr<atomgroup> atoms, unsigned long skip, bool usecom) :
        atoms(atoms), pastlocs(atoms->size(), Vec()), dists(), skip(skip), curt(0), usecom(usecom){
    Vec com = usecom ? atoms->com() : Vec();
    for(uint i = 0; i<atoms->size(); ++i){
        pastlocs[i] = (*atoms)[i].x - com;
    }
}

void RDiffs::reset(){
    curt = 0;
    for(uint i = 0; i<atoms->size(); ++i){
        pastlocs[i] = (*atoms)[i].x;
    }
    dists.clear();
}

void RDiffs::update(Box &box){
    curt++;
    if(curt < skip) return;
    Vec com = usecom ? atoms->com() : Vec();
    vector<flt> rdiff = vector<flt>(atoms->size(), 0.0);
    for(uint i = 0; i<atoms->size(); ++i){
        Vec loc = (*atoms)[i].x - com;
        rdiff[i] = (loc - pastlocs[i]).mag();
        pastlocs[i] = loc;
    }
    
    dists.push_back(rdiff);
    
    curt = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
bool jamminglist::operator<(const jamminglist& other ){
    //return distsq < other.distsq;
    if(other.distsq  - distsq > 1e-8) return true;
    if(distsq  - other.distsq > 1e-8) return false;
    //~ cout << "\nWithin 1e-8\n";
    
    uint sz = size();
    uint osz = other.size();
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for(uint i=0; i<sz; ++i){
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing sizes...\n";
    //~ cout << "Equal!\n";
    return false; // consider them equal
};

#ifdef VEC2D
/* There are two ways of looking at the different arrangements.
 * In both cases, we leave A the same as it was, and rotate / flip / translate B.
 * Also in both cases, we wrap A, then subtract off its COM (in an infinite box).
 * 
 * Method 1:
   * "unwrap" the box into 9 boxes (2D), and then choose a box for each
   * of the particles.
   * (N-1)⁹, right? (×8×N!, with rotations / flips / permutations)
 * Method 2:
   * Pick a particle (in B), move the box so that's on the left.
   * Pick another particle (maybe the same), move the box so its on the bottom.
   * Calculate COM of B without PBC, subtract that off.
   * Compare A and B, but using the PBC to calculate distances.
   * N² (×8×N!)
 * Method 3:
   * Only try different rotations, but measure distances as Σ (\vec r_{ij} - \vec s_{ij})²
 
 * We'll go with method 3.
*/

bool jamminglistrot::operator<(const jamminglistrot& other ){
    //return distsq < other.distsq;
    if(other.distsq  - distsq > 1e-8) return true;
    if(distsq  - other.distsq > 1e-8) return false;
    
    uint sz = size();
    uint osz = other.size();
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << "\nWithin 1e-8. ";
    
    //~ cout << sz << ' ' << osz << ' ' << N << " Indices:";
    for(uint i=0; i<sz; ++i){
        if (assigned[i] < other.assigned[i]) return true;
        if (assigned[i] > other.assigned[i]) return false;
        //~ cout << " " << i;
    }
    //~ cout << " Done. Comparing rotations...";
    if(rotation < other.rotation) return true;
    if(rotation > other.rotation) return false;
    
    //~ cout << " Comparing sizes...";
    if(sz < osz) return true;
    if(sz > osz) return false;
    //~ cout << "Equal!\n";
    return false; // consider them equal
};

jammingtree2::jammingtree2(sptr<Box>box, vector<Vec>& A0, vector<Vec>& B0)
            : box(box), jlists(), A(A0), Bs(8, B0){
    for(uint rot=0; rot < 8; ++rot){
        for(uint i=0; i<B0.size(); ++i){
                        Bs[rot][i] = B0[i].rotate_flip(rot); }
        if(A0.size() <= B0.size()) jlists.push_back(jamminglistrot(rot));
        //~ cout << "Created, now size " << jlists.size() << endl;
    }
    
};

flt jammingtree2::distance(jamminglistrot& jlist){
    flt dist = 0;
    uint rot = jlist.rotation;
    for(uint i=1; i<jlist.size(); ++i){
        uint si = jlist.assigned[i];
        for(uint j=0; j<i; ++j){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            dist += box->diff(rij, sij).sq();
        }
    }
    return dist / ((flt) jlist.assigned.size());
};

list<jamminglistrot> jammingtree2::expand(jamminglistrot curjlist){
    vector<uint>& curlist = curjlist.assigned;
    list<jamminglistrot> newlists = list<jamminglistrot>();
    if(curlist.size() >= A.size()){
        return newlists;
    }
    
    uint N = (uint) Bs[curjlist.rotation].size();
    for(uint i=0; i < N; ++i){
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;
        
        jamminglistrot newjlist = jamminglistrot(curjlist, i, 0);
        newjlist.distsq = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool jammingtree2::expand(){
    jamminglistrot curjlist = jlists.front();
    list<jamminglistrot> newlists = expand(curjlist);
    
    if(newlists.empty()){
        //~ cout << "No lists made\n";
        return false;
    }
    //~ cout << "Have " << newlists.size() << "\n";
    newlists.sort();
    //~ cout << "Sorted.\n";
    jlists.pop_front();
    //~ cout << "Popped.\n";
    jlists.merge(newlists);
    //~ cout << "Merged to size " << jlists.size() << "best dist now " << jlists.front().distsq << "\n";
    return true;
};


jammingtreeBD::jammingtreeBD(sptr<Box> box, vector<Vec>& A, vector<Vec>& B, 
                    uint cutoffA, uint cutoffB) :
            jammingtree2(box, A, B), cutoff1(cutoffA), cutoff2(cutoffB){
    if(cutoffA > cutoffB){jlists.clear();}
    if(A.size() - cutoffA > B.size() - cutoffB){jlists.clear();}
};

list<jamminglistrot> jammingtreeBD::expand(jamminglistrot curjlist){
    vector<uint>& curlist = curjlist.assigned;
    list<jamminglistrot> newlists = list<jamminglistrot>();
    if(curlist.size() >= A.size()){
        return newlists;
    }
    
    uint N = (uint) Bs[curjlist.rotation].size();
    uint start = 0, end=cutoff2;
    if(curlist.size() >= cutoff1){
        start=cutoff2;
        end = N;
    }
    for(uint i=start; i < end; ++i){
        vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
        //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
        if (found != curlist.end()) continue;
        
        jamminglistrot newjlist = jamminglistrot(curjlist, i, 0);
        newjlist.distsq = distance(newjlist);
        newlists.push_back(newjlist);
    }
    return newlists;
};

bool jammingtreeBD::expand(){
    jamminglistrot curjlist = jlists.front();
    list<jamminglistrot> newlists = expand(curjlist);
    
    if(newlists.empty()) return false;
    newlists.sort();
    jlists.pop_front();
    jlists.merge(newlists);
    return true;
};


vector<Vec> jammingtree2::locationsB(jamminglistrot jlist){
    uint rot = jlist.rotation;
    vector<Vec> locs = vector<Vec>(jlist.size());
    
    uint N = jlist.size();
    for(uint i=0; i<N; ++i){
        uint si = jlist.assigned[i];
        locs[i] = A[i];
        for(uint j=0; j<N; ++j){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            locs[i] -= box->diff(rij, sij)/N;
        }
    }
    return locs;
};


vector<Vec> jammingtree2::locationsA(jamminglistrot jlist){
    uint rot = jlist.rotation;
    vector<Vec> locs = vector<Vec>(Bs[rot].size(), Vec(NAN,NAN));
    
    uint N = jlist.size();
    for(uint i=0; i<N; ++i){
        uint si = jlist.assigned[i];
        locs[si] = Bs[rot][si];
        for(uint j=0; j<N; ++j){
            uint sj = jlist.assigned[j];
            Vec rij = box->diff(A[i], A[j]);
            Vec sij = box->diff(Bs[rot][si], Bs[rot][sj]);
            locs[si] += box->diff(rij, sij)/N;
        }
        
        // this is an inverse rotateflip
        locs[si] = locs[si].rotate_flip_inv(rot);
    }
    return locs;
};

Vec jammingtree2::straight_diff(Box &bx, vector<Vec>& As, vector<Vec>& Bs){
    uint N = (uint) As.size();
    if(Bs.size() != N) return Vec(NAN,NAN);
    
    Vec loc = Vec();
    for(uint i=0; i<N; ++i){
        for(uint j=0; j<N; ++j){
            Vec rij = bx.diff(As[i], As[j]);
            Vec sij = bx.diff(Bs[i], Bs[j]);
            loc += bx.diff(rij, sij);
        }
    }
    return loc / N;
};

flt jammingtree2::straight_distsq(Box &bx, vector<Vec>& As, vector<Vec>& Bs){
    uint N = (uint) As.size();
    if(Bs.size() != N) return NAN;
    
    flt dist = 0;
    for(uint i=0; i<N; ++i){
        for(uint j=0; j<N; ++j){
            Vec rij = bx.diff(As[i], As[j]);
            Vec sij = bx.diff(Bs[i], Bs[j]);
            dist += bx.diff(rij, sij).sq();
        }
    }
    return dist / N;
};



#endif

////////////////////////////////////////////////////////////////////////////////////////////////////

void CNodePath::add(CNode node, OriginBox& box){
    if(size() > 0){
        distance += box.diff(node.x, nodes.back().x);
    }
    nodes.push_back(node);
};

void Connectivity::add_edge(CNode node1, CNode node2){
    if(node1.n < 0 || node2.n < 0){
        throw invalid_argument("Got node with negative n");
    }
    nodes.insert(node1);
    nodes.insert(node2);
    
    neighbors[node1.n].push_back(node2);
    std::sort(neighbors[node1.n].begin(), neighbors[node1.n].end());
    neighbors[node2.n].push_back(node1);
    std::sort(neighbors[node2.n].begin(), neighbors[node2.n].end());
};

array<bool, NDIM> Connectivity::nonzero(Vec diff_vec){
    array<bool, NDIM> nonzeros;
    Vec half_shape = box->boxshape() / 2;
    for(uint i=0; i<NDIM; i++){
        nonzeros[i] = (abs(diff_vec[i]) > half_shape[i]);
    }
    return nonzeros;
};

void Connectivity::add(vector<Vec> locs, vector<flt> diameters){
    vector<CNode> cnodes;
    
    uint n=0;
    for(vector<Vec>::const_iterator it=locs.begin(); it!=locs.end(); it++){
        CNode cn = CNode(n, *it);
        
        uint n2=0;
        for(vector<CNode>::const_iterator cit=cnodes.begin(); cit!=cnodes.end(); cit++){
            Vec dx = box->diff(cn.x, cit->x);
            if(dx.mag() <= (diameters[n] + diameters[n2])/2.0){
                neighbors[n].push_back(*cit);
                std::sort(neighbors[n].begin(), neighbors[n].end());
                neighbors[n2].push_back(cn);
                std::sort(neighbors[n2].begin(), neighbors[n2].end());
            }
            n2++;
        }
        nodes.insert(cn);
        cnodes.push_back(cn);
        
        n++;
    }
};

CNodePath Connectivity::make_cycle(CNodePath forward, CNodePath backward){
    //forward and backward should both start and end in the same place
    //and vice versa
    
    if(forward.size() == 0)
        throw invalid_argument("make_cycle got an empty forward path. This should not happen.");
    else if(backward.size() == 0)
        throw invalid_argument("make_cycle got an empty forward path. This should not happen.");
    
    if(forward.nodes.front() != backward.nodes.front())
        throw invalid_argument("make_cycle: forward.front() != backward.front()");
    if(backward.nodes.back() != forward.nodes.back())
        throw invalid_argument("make_cycle: backward.back() != forward.back()");
    
    CNodePath cycle = forward;
    
    for(vector<CNode>::reverse_iterator it=++(backward.nodes.rbegin()); it != backward.nodes.rend(); it++){
        cycle.add(*it, *box);
    }
    return cycle;
};

map<uint, CNodePath> Connectivity::circular_from(CNode node, set<uint>& visited, bool check_all){
    map<uint, CNodePath> full_paths; // complete roots around the box. key is DIMENSION
    map<uint, CNodePath> prev_paths; // other ways we've found to get to any node
    queue<CNodePath> paths;
    
    CNodePath path0 = CNodePath(node);
    paths.push(path0);
    prev_paths[node.n] = path0;
    visited.insert(node.n);
    
    while(!paths.empty()){
        CNodePath path = paths.front();
        paths.pop();
        assert(path.size() > 0);
        CNode lastnode = path.nodes.back();
        vector<CNode>& nextnodes = neighbors[lastnode.n];
        for(vector<CNode>::iterator it=nextnodes.begin(); it < nextnodes.end(); ++it){
            if(it->n < node.n) continue;
            CNode nextnode = *it;
            CNodePath newpath = CNodePath(path, nextnode, *box);
            map<uint, CNodePath>::iterator found_path_it = prev_paths.find(nextnode.n);
            if(found_path_it != prev_paths.end()){
                // We've been to this node before.
                // But have we found a circle (around the box), or just another route?
                bool found_nonzero = false;
                CNodePath& found_path = found_path_it->second;
                
                Vec pathdiff = found_path.distance - newpath.distance;
                array<bool, NDIM> nonzeros = nonzero(pathdiff);
                for(uint i=0; i<NDIM; i++){
                    if(nonzeros[i]){
                        found_nonzero = true;
                        // We've found a cycle around the whole box!
                        map<uint, CNodePath>::iterator found_cycle = full_paths.find(i);
                        CNodePath cycle_path = make_cycle(found_path, newpath);
                        if(found_cycle == full_paths.end()){
                            // And its along a dimension we've never found before
                            // Note that now both *found_path and newpath
                            // begin with CNode "node" and end with "nextnode"
                            
                            full_paths[i] = cycle_path;
                            // Have we found enough paths?
                            // If so, we're done. We're not trying to find the "best" cycle,
                            // just any cycle (or any NDIM cycles if check_all)
                            if(!check_all) return full_paths;
                            else if(full_paths.size() >= NDIM) return full_paths;
                        } else {
                            // Already found a cycle in this dimension, maybe replace it
                            if(full_paths[i].nodes.size() > cycle_path.nodes.size())
                                full_paths[i] = cycle_path;
                        }
                    }
                }
                
                if(!found_nonzero) {
                    // Different route to the same node, but through the same box.
                    // If this is a shorter route, might as well hold onto it...
                    if(newpath.size() < found_path.size()){
                        prev_paths[nextnode.n] = newpath;
                    }
                    continue;
                }
            } else {
                // we've never been to this node before
                visited.insert(nextnode.n);
                prev_paths[nextnode.n] = newpath;
                paths.push(newpath);
                continue;
            }
            // Now same path, new neighbor to add to it
        }
        // We're done, move on to the next incomplete path
    };
    
    return full_paths;
};

map<uint, CNodePath> Connectivity::find_percolation(bool check_all){
    map<uint, CNodePath> full_paths;
    set<uint> visited;
    for(set<CNode>::iterator it=nodes.begin(); it!=nodes.end(); ++it){
        if(visited.count(it->n) > 0) continue; 
        // we've already been to this node, no need to start there
        
        map<uint, CNodePath> new_paths = circular_from(*it, visited, check_all);
        for(map<uint, CNodePath>::iterator fit=new_paths.begin(); fit!=new_paths.end(); ++fit){
            uint dim = fit->first;
            CNodePath& new_path = fit->second;
            if(full_paths.find(dim) == full_paths.end()){
                // never been along this dimension before
                full_paths[dim] = new_path;
            } else {
                // found one along this dimension before. If the new one is shorter, might as
                // well take it instead
                CNodePath& old_path = full_paths[dim];
                if(new_path.nodes.size() < old_path.nodes.size()){
                    full_paths[dim] = new_path;
                }
            }
        }
        if(!check_all and !full_paths.empty()) return full_paths;
        else if(full_paths.size() >= NDIM) return full_paths;
    }
    return full_paths;
};
