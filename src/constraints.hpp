#include "interaction.hpp"

#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <complex>

using namespace std;

typedef std::complex<flt> cmplx; // need the std:: for SWIG complex.i, not sure why

class Constraint : public boost::enable_shared_from_this<Constraint> {
    public:
        virtual void apply(Box &box) = 0;
        virtual int ndof() = 0;
        virtual ~Constraint(){};
};

class CoordConstraint : public Constraint {
    private:
        Atom* a;
        bool fixed[3];
        Vec loc;
    public:
        CoordConstraint(Atom* atm, bool fixx, bool fixy, bool fixz, Vec loc) :
            a(atm), loc(loc) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        CoordConstraint(Atom* atm, bool fixx, bool fixy, bool fixz) :
            a(atm), loc(a->x) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        CoordConstraint(Atom* atm) :
            a(atm), loc(a->x) {fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            for(uint i=0; i<3; i++){
                if(not fixed[i]) continue;
                a->f[i] = 0;
                a->v[i] = 0;
                a->x[i] = loc[i];
            }
        }
};

class CoordCOMConstraint : public Constraint {
    private:
        sptr<AtomGroup> a;
        bool fixed[3];
        Vec loc;
    public:
        CoordCOMConstraint(sptr<AtomGroup> atm, bool fixx, bool fixy, bool fixz, Vec loc) :
            a(atm), loc(loc) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        CoordCOMConstraint(sptr<AtomGroup> atm, bool fixx, bool fixy, bool fixz) :
            a(atm), loc(a->com()) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        CoordCOMConstraint(sptr<AtomGroup> atm) :
            a(atm), loc(a->com()) {fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            Vec com = a->com() - loc;
            Vec comv = a->comv();
            Vec totf = Vec();
            for(uint i=0; i< a->size(); i++){
                totf += (*a)[i].f;
            }
            Vec tota = totf / a->mass();

            for(uint i=0; i< a->size(); i++){
                Atom &atm = (*a)[i];
                Vec df = (tota * (atm.m));
                for(uint j=0; j<3; j++){
                    if(not fixed[j]) continue;
                    atm.f[j] -= df[j];
                    atm.v[j] -= comv[j];
                    atm.x[j] -= com[j];
                }
            }
        }
};

class RelativeConstraint : public Constraint {
    private:
        Atom *a1, *a2;
        bool fixed[3];
        Vec loc;
    public:
        RelativeConstraint(Atom* atm1, Atom* atm2, bool fixx, bool fixy, bool fixz, Vec loc) :
            a1(atm1), a2(atm2), loc(loc) {
                fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        RelativeConstraint(Atom* atm1, Atom* atm2, bool fixx, bool fixy, bool fixz) :
            a1(atm1), a2(atm2), loc(a2->x - a1->x) {
                fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        RelativeConstraint(Atom* atm1, Atom* atm2) :
            a1(atm1), a2(atm2), loc(a2->x - a1->x) {
                fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            flt mratio1 = a1->m / (a1->m + a2->m);
            flt mratio2 = a2->m / (a1->m + a2->m);
            Vec totf = a2->f + a1->f;
            Vec dv = a2->v - a1->v;
            Vec dx = a2->x - a1->x;
            for(uint i=0; i<NDIM; i++){
                if(not fixed[i]) continue;
                a1->f[i] = totf[i]*mratio1;
                a2->f[i] = totf[i]*mratio2;
                a1->v[i] += dv[i]/2;
                a2->v[i] -= dv[i]/2;
                a1->x[i] += dx[i]/2;
                a2->x[i] -= dx[i]/2;
                assert(abs(a2->v[i] - a1->v[i]) < 1e-5);
                assert(abs(a2->x[i] - a1->x[i]) < 1e-5);
            }
        }
};

class DistConstraint : public Constraint {
    private:
        AtomID a1, a2;
        flt dist;
    public:
        DistConstraint(AtomID atm1, AtomID atm2, flt dist) :
            a1(atm1), a2(atm2), dist(dist) {};
        DistConstraint(AtomID atm1, AtomID atm2) :
            a1(atm1), a2(atm2), dist((a1->x - a2->x).mag()){};
        int ndof(){return 1;};
        void apply(Box &box){
            flt M = (a1->m + a2->m);
            flt mratio1 = a1->m / M;
            flt mratio2 = a2->m / M;

            Vec dx = a2->x - a1->x;
            flt dxmag = dx.mag();
            Vec dxnorm = dx / dxmag;

            a1->x += dx * ((1 - dist/dxmag)*mratio2);
            a2->x -= dx * ((1 - dist/dxmag)*mratio1);
            //~ dx = a2->x - a1->x;
            //~ dxmag = dx.mag();
            //~ dxnorm = dx / dxmag; dxnorm should still be the same

            Vec baddv = dxnorm * ((a2->v - a1->v).dot(dxnorm)/2);
            a1->v += baddv;
            a2->v -= baddv;

            // newv2 • u = v2 - |baddv|
            // newv1 • u = v1 + |baddv|
            // (newv2 - newv1) • u = (v2 - v1 - (2*baddv)) • u
            //                     = ((v2 - v1)•u - (2*baddv)•u)
            //                     = (|baddv|*2 - |baddv|*2) = 0
            // assert((a2->v - a1->v).dot(dxnorm) < 1e-8);
            if((a2->v - a1->v).dot(dxnorm) > 1e-8){
                throw std::overflow_error("Velocities are not minimal.");
            }

            // TODO: Fix mass ratio stuff
            Vec baddf = dxnorm * ((a2->f - a1->f).dot(dxnorm)/2);
            a1->f += baddf;
            a2->f -= baddf;
            // assert((a2->f - a1->f).dot(dxnorm) < 1e-8);
            if((a2->f - a1->f).dot(dxnorm) > 1e-8){
                throw std::overflow_error("Forces are not minimal.");
            }

        }
};


class LinearConstraint : public Constraint {
    private:
        sptr<AtomGroup> atms;
        flt dist;
        flt lincom, I, M;
    public:
        LinearConstraint(sptr<AtomGroup> atms, flt dist) :
            atms(atms), dist(dist), lincom(0), I(0), M(0) {
            for(uint i = 0; i < atms->size(); i++){
                M += (*atms)[i].m;
                lincom += (dist*i)*(*atms)[i].m;
            }
            lincom /= M;

            for(uint i = 0; i < atms->size(); i++){
                flt dx = (dist*i - lincom);
                I += (*atms)[i].m * dx * dx;
            }
        };
        int ndof(){return (int)atms->size()-1;};

        void apply(Box &box){
            Vec com = atms->com();
            Vec comv = atms->comv();
            Vec comf = Vec();

            uint sz = atms->size();
            Vec lvec = Vec();
            #ifdef VEC3D
            Vec L = Vec();
            Vec omega  = Vec();
            Vec tau = Vec();
            Vec alpha = Vec();
            #else
            flt L = 0;
            flt omega = 0;
            flt tau = 0;
            flt alpha = 0;
            #endif

            for(uint i = 0; i < sz; i++){
                flt chaindist = i * dist - lincom;
                Atom& ai = (*atms)[i];
                Vec dx = ai.x - com;
                comf += ai.f;
                lvec += dx.norm() * chaindist;
                L += dx.cross(ai.v) * ai.m;
                tau += dx.cross(ai.f);
            }

            lvec.normalize();
            omega = L / I;
            alpha = tau / I;

            for(uint i = 0; i < sz; i++){
                flt chaindist = i * dist - lincom;
                Vec dx = lvec*chaindist;
                Atom& ai = (*atms)[i];
                ai.x = com + dx;
                ai.v = comv + dx.cross(omega);
                ai.f = comf + (dx.cross(alpha)*ai.m);
            }
        }
};

class ContactTracker : public StateTracker{
    protected:
        sptr<AtomGroup> atoms;
        vector<flt> dists;
        vector<vector<bool> > contacts;

        unsigned long long breaks;
        unsigned long long formations;
        unsigned long long incontact;
    public:
        ContactTracker(sptr<Box> box, sptr<AtomGroup> atoms, vector<flt> dists);
        void update(Box &box);

        void reset(){
            breaks=0;
            formations=0;
            incontact=0;
            uint N = atoms->size();
            contacts.resize(N);
            for(uint i=0; i<N; ++i){
                contacts[i].assign(i, false);
            }
        };
        unsigned long long broken(){return breaks;};
        unsigned long long formed(){return formations;};
        unsigned long long number(){return incontact;};
};

inline ContactTracker* ContactTrackerD(sptr<Box> box, sptr<AtomGroup> atoms, vector<double> dists){
    vector<flt> newdists = vector<flt>();
    for(uint i=0; i<dists.size(); i++){
        newdists.push_back(dists[i]);
    }
    return new ContactTracker(box, atoms, newdists);
}

class EnergyTracker : public StateTracker{
    protected:
        sptr<AtomGroup> atoms;
        vector<sptr<Interaction> > interactions;

        uint N;
        uint nskip, nskipped;
        flt U0;
        flt Es, Us, Ks;
        flt Esq, Usq, Ksq;
    public:
        EnergyTracker(sptr<AtomGroup> atoms,
            vector<sptr<Interaction> > interactions, uint nskip=1)
             : atoms(atoms),
            interactions(interactions), N(0), nskip(max(nskip,1u)), nskipped(0),
            U0(0),Es(0),Us(0),Ks(0), Esq(0), Usq(0), Ksq(0){};
        void update(Box &box);
        void reset(){
            nskipped=0;
            N=0; Es=0; Us=0; Ks=0;
            Esq=0; Usq=0; Ksq=0;
        };
        void setU0(flt newU0){
            U0 = newU0;
            reset();
        };
        void setU0(Box &box);
        flt getU0(){return U0;};

        flt E(){return Es/((flt) N);};
        flt U(){return Us/((flt) N);};
        flt K(){return Ks/((flt) N);};
        flt Estd(){return sqrt(Esq/N -Es*Es/N/N);};
        flt Kstd(){return sqrt(Ksq/N -Ks*Ks/N/N);};
        flt Ustd(){return sqrt(Usq/N -Us*Us/N/N);};
        flt Esqmean(){return Esq/N;};
        flt Ksqmean(){return Ksq/N;};
        flt Usqmean(){return Usq/N;};
        //~ flt Ustd(){return sqrt((Usq -(U*U)) / ((flt) N));};
        //~ flt Kstd(){return sqrt((Ksq -(K*K)) / ((flt) N));};
        uint n(){return N;};
};

class RsqTracker1 {
    // Tracks only a single dt (skip)
    public:
        vector<Vec> pastlocs;
        vector<Vec> xyz2sums;
        vector<Vec> xyz4sums;
        vector<flt> r4sums;
        unsigned long skip, count;

    public:
        RsqTracker1(AtomGroup& atoms, unsigned long skip, Vec com);

        void reset(AtomGroup& atoms, Vec com);

        bool update(Box& box, AtomGroup& atoms, unsigned long t, Vec com); // updates if necessary.
        vector<Vec> xyz2();
        vector<Vec> xyz4();
        vector<flt> r4();

        unsigned long get_skip(){return skip;};
        unsigned long get_count(){return count;};
};

class RsqTracker : public StateTracker {
    public:
        sptr<AtomGroup> atoms;
        vector<RsqTracker1> singles;
        unsigned long curt;
        bool usecom;

    public:
        RsqTracker(sptr<AtomGroup> atoms, vector<unsigned long> ns, bool usecom=true);

        void reset();
        void update(Box &box);

        vector<vector<Vec> > xyz2();
        vector<vector<flt> > r2();
        vector<vector<Vec> > xyz4();
        vector<vector<flt> > r4();
        vector<flt> counts();
};

////////////////////////////////////////////////////////////////////////
// ISF tracking
// code is similar to Rsqtracker.
// It tracks ISF(k, Δt) with one ISFTracker1 per Δt.
// k is of type 'flt', representing a length; it will average over
// k(x hat), k(y hat), k(z hat).

class ISFTracker1 {
    // Tracks only a single dt (skip)
    public:
        vector<Vec> pastlocs;
        vector<vector<Array<cmplx, NDIM> > > ISFsums; // (number of ks x number of particles x number of dimensions)
        vector<flt> ks;
        unsigned long skip, count;

    public:
        ISFTracker1(AtomGroup& atoms, unsigned long skip, vector<flt> ks, Vec com);

        void reset(AtomGroup& atoms, Vec com);

        bool update(Box& box, AtomGroup& atoms, unsigned long t, Vec com); // updates if necessary.
        vector<vector<cmplx> > ISFs();
        vector<vector<Array<cmplx, NDIM> > > ISFxyz();

        unsigned long get_skip(){return skip;};
        unsigned long get_count(){return count;};
};

class ISFTracker : public StateTracker {
    public:
        sptr<AtomGroup> atoms;
        vector<ISFTracker1> singles;
        unsigned long curt;
        bool usecom;

    public:
        ISFTracker(sptr<AtomGroup> atoms, vector<flt> ks,
                    vector<unsigned long> ns, bool usecom=false);

        void reset();
        void update(Box &box);

        vector<vector<vector<Array<cmplx, NDIM> > > > ISFxyz();
        vector<vector<vector<cmplx> > > ISFs();
        vector<flt> counts();
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Get the "smoothed" version of each particles location, "smoothed" over windows of a given size

class SmoothLocs : public StateTracker {
    public:
        sptr<AtomGroup> atoms;
        uint smoothn; // number of skipns to "smooth" over
        uint skipn; // number of dts to skip
        vector<Vec> curlocs;
        uint numincur;
        vector<vector<Vec> > locs;
        unsigned long curt;
        bool usecom;

    public:
        SmoothLocs(sptr<AtomGroup> atoms, Box &box, uint smoothn, uint skipn=1, bool usecom=false);

        void reset();
        void update(Box &box);

        vector<vector<Vec> > smooth_locs(){return locs;};
};

////////////////////////////////////////////////////////////////////////////////////////////////////
class RDiffs : public StateTracker {
    // Tracks only a single dt (skip)
    public:
        sptr<AtomGroup> atoms;
        vector<Vec> pastlocs;
        vector<vector<flt> > dists;
        unsigned long skip;
        unsigned long curt;
        bool usecom;
    public:
        RDiffs(sptr<AtomGroup> atoms, unsigned long skip, bool usecom=false);

        void reset();

        void update(Box& box); // updates if necessary.
        vector<vector<flt> > rdiffs(){return dists;};

        unsigned long get_skip(){return skip;};
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// For comparing two jammed structures

/* We have two packings, A and B, and want to know the sequence {A1, A2, A3...}
 * such that particle A1 of packing 1 matches particle 1 of packing B.
 * A JammingList is a partial list; it has a list {A1 .. An}, with n / N
 * particles assigned, with a total distance² of distsq.
*/
class JammingList {
    public:
        vector<uint> assigned;
        flt distsq;

        JammingList() : assigned(), distsq(0){};
        JammingList(const JammingList& other)
            : assigned(other.assigned), distsq(other.distsq){};
        JammingList(const JammingList& other, uint expand, flt addeddist)
            : assigned(other.size() + 1, 0), distsq(other.distsq + addeddist){
            for(uint i=0; i < other.size(); i++){
                assigned[i] = other.assigned[i];
            }
            assigned[assigned.size()-1] = expand;
        }
        inline uint size() const {return (uint) assigned.size();};

        bool operator<(const JammingList& other) const;
};

class JammingTree {
    private:
        sptr<Box> box;
        list<JammingList> jlists;
        vector<Vec> A;
        vector<Vec> B;
    public:
        JammingTree(sptr<Box> box, vector<Vec>& A, vector<Vec>& B)
            : box(box), jlists(), A(A), B(B) {
            jlists.push_back(JammingList());
            assert(A.size() <= B.size());
        };

        bool expand(){
            JammingList curjlist = jlists.front();
            vector<uint>& curlist = curjlist.assigned;
            if(curlist.size() >= A.size()){
                //~ cout << "List already too big\n";
                return false;
            }

            list<JammingList> newlists = list<JammingList>();
            for(uint i=0; i < B.size(); i++){
                vector<uint>::iterator found = find(curlist.begin(), curlist.end(), i);
                //if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
                if (found != curlist.end()){
                    //~ cout << "Found " << i << "\n";
                    //cout << found << '\n';
                    continue;
                }
                flt newdist = box->diff(A[curlist.size()], B[i]).sq();
                JammingList newjlist = JammingList(curjlist, i, newdist);
                newlists.push_back(newjlist);
                //~ cout << "Made " << i << "\n";
            }

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
        }
        bool expand(uint n){
            bool retval=false;
            for(uint i=0; i<n; i++){
                retval = expand();
            }
            return retval;
        }
        list<JammingList> &mylist(){return jlists;};
        list<JammingList> copylist(){return jlists;};

        JammingList curbest(){
            JammingList j = JammingList(jlists.front());
            //~ cout << "Best size: " << j.size() << " dist: " << j.distsq;
            //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " << j.assigned[j.size()-1] << "]";
            //~ cout << '\n';
            return j;
            //return JammingList(jlists.front());
            };
        uint size(){return (uint) jlists.size();};
};

#ifdef VEC2D

class JammingListRot : public JammingList {
    public:
        uint rotation;

        JammingListRot() : JammingList(), rotation(0){};
        JammingListRot(uint rot) : JammingList(), rotation(rot){};
        JammingListRot(const JammingListRot& other)
            : JammingList(other), rotation(other.rotation){};
        JammingListRot(const JammingListRot& other, uint expand, flt addeddist)
            : JammingList(other, expand, addeddist), rotation(other.rotation){};

        bool operator<(const JammingListRot& other) const;
};

// Includes rotations, flips, and translations.
class JammingTree2 {
    protected:
        sptr<Box> box;
        list<JammingListRot> jlists;
        vector<Vec> A;
        vector<vector<Vec> > Bs;
    public:
        // make all 8 possible rotations / flips
        // then subtract off all possible COMVs
        JammingTree2(sptr<Box>box, vector<Vec>& A, vector<Vec>& B);
        flt distance(JammingListRot& jlist);
        list<JammingListRot> expand(JammingListRot curjlist);

        virtual bool expand();

        bool expand(uint n){
            bool retval=false;
            for(uint i=0; i<n; i++){
                retval = expand();
                if(!retval) break;
            }
            return retval;
        }
        bool expandto(flt maxdistsq){
            bool retval = true;
            while((maxdistsq <= 0 or jlists.front().distsq < maxdistsq) and retval){
                retval = expand();
            };
            return retval;
        }
        static Vec straight_diff(Box &bx, vector<Vec>& A, vector<Vec>& B);
        static flt straight_distsq(Box &bx, vector<Vec>& A, vector<Vec>& B);

        list<JammingListRot> &mylist(){return jlists;};
        list<JammingListRot> copylist(){return jlists;};
        list<JammingListRot> copylist(uint n){
            list<JammingListRot>::iterator last = jlists.begin();
            advance(last, n);
            return list<JammingListRot>(jlists.begin(), last);
        };


        JammingListRot curbest(){
            if(jlists.empty()){
                JammingListRot bad_list = JammingListRot();
                bad_list.distsq = -1;
                return bad_list;
                }
            JammingListRot j = JammingListRot(jlists.front());
            //~ cout << "Best size: " << j.size() << " dist: " << j.distsq;
            //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " << j.assigned[j.size()-1] << "]";
            //~ cout << '\n';
            return j;
            //return JammingList(jlists.front());
            };

        //JammingListRot operator[](uint i){
        //    assert(i < jlists.size());
        //    return JammingListRot(jlists[i]);
        //};

        uint size(){return (uint) jlists.size();};

        vector<Vec> locationsB(JammingListRot jlist);
        vector<Vec> locationsB(){return locationsB(curbest());};
        vector<Vec> locationsA(JammingListRot jlist);
        vector<Vec> locationsA(){return locationsA(curbest());};
        virtual ~JammingTree2(){};
};


class JammingTreeBD : public JammingTree2 {
    /* For a bi-disperse packing.
     * 'cutoff' is the number of particles of the first kind; i.e., the
     * A vector should have A[0]..A[cutoff-1] be of particle type 1,
     * and A[cutoff]..A[N-1] of particle type 2.
     * This does much the same as JammingTree2, but doesn't check any
     * reordering in which particles of one type are relabeled as another.
     * For exampe, with 2+2 particles (cutoff 2), we check
     * [0123],[1023],[0132],[1032]
     * But not
     * [0213],[0231],[0312],[0321],[1203],[1230],[1302],[1320],...
     * This means at most (cutoff! (N-cutoff)!) combinations are checked,
     * and not all N!, which can save a lot of time (as well as
     *  rejecting false combinations).
     */
    protected:
        uint cutoff1,cutoff2;
    public:
        JammingTreeBD(sptr<Box>box, vector<Vec>& A, vector<Vec>& B, uint cutoff) :
            JammingTree2(box, A, B), cutoff1(cutoff), cutoff2(cutoff){};
        JammingTreeBD(sptr<Box>box, vector<Vec>& A, vector<Vec>& B,
                    uint cutoffA, uint cutoffB);// :
            //JammingTree2(box, A, B), cutoff1(cutoffA), cutoff2(cutoffB){};

        list<JammingListRot> expand(JammingListRot curjlist);
        bool expand();
        bool expand(uint n){return JammingTree2::expand(n);};
};
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Finding Percolation
 *
 * Plan:
 * struct Node {
 *      uint n;
 *      Vec x;
 * }
 *
 * struct Connectivity {
 *      sptr<OriginBox>;
 *      vector<Node> nodes;
 *      set<Node, vector<Node> > neighbors; # with both directions in
 * }
 *
 * make_connectivity(sptr<OriginBox>, NListed<A,P>) {
 *      # for each pair, if energy != 0, add it to the Connectivity
 * }
 *
 * or
 *
 * make_connectivity(sptr<OriginBox>, NeighborList, ... sigmas) {
 *      # for each pair, if box.diff(a1.x, a2.x) < (sigma1 < sigma2)/2, add it to the list
 * }
 *
 * struct Path {
 *      Vec<uint> distance          # Euclidean distance from first node to last. Note this is not
 *                                  #   the same as box.diff(last - first), because it might go
 *                                  #   a longer way around the box
 *      vector<uint> nodes          # nodes visited
 * }
 *
 * Maybe implement <, ==, > by using length of nodes first and what the nodes are second
 *
 * Implement ... circular_from(uint n, check_all=true):
 *  -   Do a breadth-first search starting from n, including only nodes > n.
 *  -   While searching, build map<node, path>, which is the Euclidean distance from the start node
 *      to the current node via the path traveled
 *  -   If you get to an already visited node, compare the current Vec path to the previous one. If
 *      the same:
 *      -   either replace the old path (if new_path < old_path) or don't
 *      -   Do not follow connections
 *      If different, we've found a percolation:
 *      -   Mark percolations with "map<uint, path> cycles"; the key is the dimension, and
 *          the value (vector<uint>) is the path.
 *      -   If(check_all && cycles.size() >= 1) return cycles
 *      -   else if cycles.size() >= NDIM return cycles
 *  -   If you finish without finding any (or enough), return cycles
 *
 * Implement ... find_percolation(...):
 *  -   for each n, run circular_from(uint n, false)
 *  -   if any have a cycle, return it!
 *
 * Implement ... find_percolations(...):
 *  -   for each n, run circular_from(uint n, true)
 *  -   keep merging results
 *  -   if cycles.size() >= NDIM return cycles
 *  -   if you get to the end, return cycles
 */

class CNode {
    public:
        int n;
        Vec x;

        CNode() : n(-1){};
        CNode(int n, Vec x) : n(n), x(x){};

        bool operator!=(const CNode& other)  const { return n != other.n;};
        bool operator==(const CNode& other) const { return n == other.n;};
        bool operator<(const CNode& other) const { return n < other.n;};
        bool operator>(const CNode& other) const { return n > other.n;};
};

class CNodePath {
    public:
        Vec distance;
        vector<CNode> nodes;

    public:
        CNodePath(){};
        CNodePath(CNode node){nodes.push_back(node);};
        CNodePath(CNodePath other, CNode node, OriginBox& box) :
            distance(other.distance), nodes(other.nodes){add(node, box);};
        void add(CNode node, OriginBox& box);
        uint size(){return (uint) nodes.size();};
};

// The class that actually finds percolation.

// TODO:
/*
 * add a way to work with a NeighborList for initial construction
 */
class Connectivity {
    public:
        sptr<OriginBox> box;
        set<CNode> nodes;
        map<int, vector<CNode> > neighbors;

        array<bool, NDIM> nonzero(Vec diff_vec); // is it nonzero. Specifically, is any dimension >  L/2?
        CNodePath make_cycle(CNodePath forward, CNodePath backward);

        map<uint, CNodePath> circular_from(CNode node, set<uint>& visited, bool check_all);

    public:
        Connectivity(sptr<OriginBox> box) : box(box){};
        void add_edge(CNode node1, CNode node2);

        // assumes diameters are additive
        // Note that "diameter" should generally be the attractive diameter
        void add(vector<Vec> locs, vector<flt> diameters);

        map<uint, CNodePath> find_percolation(bool check_all_dims=true);
};

#endif
