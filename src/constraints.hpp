#include "interaction.hpp"

#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <complex>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <vector>

using namespace std;

class Constraint {
   public:
    virtual void apply_positions(Box& box) = 0;
    virtual void apply_velocities(Box& box) = 0;
    virtual void apply_forces(Box& box) = 0;
    virtual int constrained_dof() = 0;
    virtual ~Constraint(){};
};

class CoordConstraint : public Constraint {
   private:
    Atom* a;
    bool fixed[3];
    Vec loc;

   public:
    CoordConstraint(Atom* atm, bool fixx, bool fixy, bool fixz, Vec loc)
        : a(atm), loc(loc) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    CoordConstraint(Atom* atm, bool fixx, bool fixy, bool fixz)
        : a(atm), loc(a->x) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    CoordConstraint(Atom* atm) : a(atm), loc(a->x) {
        fixed[0] = fixed[1] = fixed[2] = true;
    };
    int constrained_dof() {
        return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];
    };
    void apply_positions(Box& box) {
        for (uint i = 0; i < 3; i++) {
            if (not fixed[i]) continue;
            a->x[i] = loc[i];
        }
    }
    void apply_velocities(Box& box) {
        for (uint i = 0; i < 3; i++) {
            if (not fixed[i]) continue;
            a->v[i] = 0;
        }
    }
    void apply_forces(Box& box) {
        for (uint i = 0; i < 3; i++) {
            if (not fixed[i]) continue;
            a->f[i] = 0;
        }
    }
};

class CoordCOMConstraint : public Constraint {
   private:
    sptr<AtomGroup> a;
    bool fixed[3];
    Vec loc;

   public:
    CoordCOMConstraint(sptr<AtomGroup> atm, bool fixx, bool fixy, bool fixz,
                       Vec loc)
        : a(atm), loc(loc) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    CoordCOMConstraint(sptr<AtomGroup> atm, bool fixx, bool fixy, bool fixz)
        : a(atm), loc(a->com()) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    CoordCOMConstraint(sptr<AtomGroup> atm) : a(atm), loc(a->com()) {
        fixed[0] = fixed[1] = fixed[2] = true;
    };
    int constrained_dof() {
        return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];
    };
    void apply_positions(Box& box);
    void apply_velocities(Box& box);
    void apply_forces(Box& box);
};

class RelativeConstraint : public Constraint {
   private:
    Atom *a1, *a2;
    bool fixed[3];
    Vec loc;

   public:
    RelativeConstraint(Atom* atm1, Atom* atm2, bool fixx, bool fixy, bool fixz,
                       Vec loc)
        : a1(atm1), a2(atm2), loc(loc) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    RelativeConstraint(Atom* atm1, Atom* atm2, bool fixx, bool fixy, bool fixz)
        : a1(atm1), a2(atm2), loc(a2->x - a1->x) {
        fixed[0] = fixx;
        fixed[1] = fixy;
        fixed[2] = fixz;
    };
    RelativeConstraint(Atom* atm1, Atom* atm2)
        : a1(atm1), a2(atm2), loc(a2->x - a1->x) {
        fixed[0] = fixed[1] = fixed[2] = true;
    };
    int constrained_dof() {
        return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];
    };

    void apply_positions(Box& box);
    void apply_velocities(Box& box);
    void apply_forces(Box& box);
};

class DistConstraint : public Constraint {
   private:
    AtomID a1, a2;
    flt dist;

   public:
    DistConstraint(AtomID atm1, AtomID atm2, flt dist)
        : a1(atm1), a2(atm2), dist(dist){};
    DistConstraint(AtomID atm1, AtomID atm2)
        : a1(atm1), a2(atm2), dist((a1->x - a2->x).norm()){};
    int constrained_dof() { return 1; };
    void apply_positions(Box& box);
    void apply_velocities(Box& box);
    void apply_forces(Box& box);
};

class LinearConstraint : public Constraint {
   private:
    sptr<AtomGroup> atms;
    flt dist;
    flt lincom, I, M;
    Vec lvec, com;

    void set_lvec_com();

   public:
    LinearConstraint(sptr<AtomGroup> atms, flt dist)
        : atms(atms),
          dist(dist),
          lincom(0),
          I(0),
          M(0),
          lvec(Vec::Zero()),
          com(Vec::Zero()) {
        for (uint i = 0; i < atms->size(); i++) {
            M += (*atms)[i].m;
            lincom += (dist * i) * (*atms)[i].m;
        }
        lincom /= M;

        for (uint i = 0; i < atms->size(); i++) {
            flt dx = (dist * i - lincom);
            I += (*atms)[i].m * dx * dx;
        }

        set_lvec_com();
    };
    int constrained_dof() { return (int)atms->size() - 1; };
    void apply_positions(Box& box);
    void apply_velocities(Box& box);
    void apply_forces(Box& box);
};

#ifdef VEC3D
//! A class that enforces rigid-body dynamics
class RigidConstraint : public Constraint {
   private:
    sptr<AtomGroup> atms;
    flt M;

    // moment of inertia matrix
    Matrix MoI;
    // moment of inertia matrix inverse
    Eigen::JacobiSVD<Matrix> MoI_solver;
    // locations relative to center of mass, no rotation
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> expected;

    // updated at every apply_positions()
    // Curent rotation matrix
    Matrix rot;
    // center of mass
    Vec com;
    // Angular velocity, acceleration
    Vec omega, alpha;

   public:
    RigidConstraint(sptr<Box> box, sptr<AtomGroup> atms);

    // TODO: should probably handle cases with atms.size() < 3
    int constrained_dof() { return (int)atms->size() * NDIM - 6; };

    void apply_positions(Box& box);
    void apply_velocities(Box& box);
    void apply_forces(Box& box);
    Matrix get_rotation();
    Matrix get_MoI() { return MoI; };
};
#endif

class ContactTracker : public StateTracker {
   protected:
    sptr<AtomGroup> atoms;
    vector<flt> dists;
    vector<vector<bool> > contacts;

    unsigned long long breaks;
    unsigned long long formations;
    unsigned long long incontact;

   public:
    ContactTracker(sptr<Box> box, sptr<AtomGroup> atoms, vector<flt> dists);
    void update(Box& box);

    void reset() {
        breaks = 0;
        formations = 0;
        incontact = 0;
        uint N = atoms->size();
        contacts.resize(N);
        for (uint i = 0; i < N; ++i) {
            contacts[i].assign(i, false);
        }
    };
    unsigned long long broken() { return breaks; };
    unsigned long long formed() { return formations; };
    unsigned long long number() { return incontact; };
};

inline ContactTracker* ContactTrackerD(sptr<Box> box, sptr<AtomGroup> atoms,
                                       vector<double> dists) {
    vector<flt> newdists = vector<flt>();
    for (uint i = 0; i < dists.size(); i++) {
        newdists.push_back(dists[i]);
    }
    return new ContactTracker(box, atoms, newdists);
}

class EnergyTracker : public StateTracker {
   protected:
    sptr<AtomGroup> atoms;
    vector<sptr<Interaction> > interactions;

    uint N;
    uint n_skip, n_skipped;
    flt U0;
    flt Es, Us, Ks;
    flt Esq, Usq, Ksq;

   public:
    EnergyTracker(sptr<AtomGroup> atoms,
                  vector<sptr<Interaction> > interactions, uint n_skip = 1)
        : atoms(atoms),
          interactions(interactions),
          N(0),
          n_skip(max(n_skip, 1u)),
          n_skipped(0),
          U0(0),
          Es(0),
          Us(0),
          Ks(0),
          Esq(0),
          Usq(0),
          Ksq(0){};
    void update(Box& box);
    void reset() {
        n_skipped = 0;
        N = 0;
        Es = 0;
        Us = 0;
        Ks = 0;
        Esq = 0;
        Usq = 0;
        Ksq = 0;
    };
    void set_U0(flt newU0) {
        U0 = newU0;
        reset();
    };
    void set_U0(Box& box);
    flt get_U0() { return U0; };

    flt E() { return Es / ((flt)N); };
    flt U() { return Us / ((flt)N); };
    flt K() { return Ks / ((flt)N); };
    flt E_std() { return sqrt(Esq / N - Es * Es / N / N); };
    flt K_std() { return sqrt(Ksq / N - Ks * Ks / N / N); };
    flt U_std() { return sqrt(Usq / N - Us * Us / N / N); };
    flt E_squared_mean() { return Esq / N; };
    flt K_squared_mean() { return Ksq / N; };
    flt U_squared_mean() { return Usq / N; };
    //~ flt U_std(){return sqrt((Usq -(U*U)) / ((flt) N));};
    //~ flt K_std(){return sqrt((Ksq -(K*K)) / ((flt) N));};
    uint n() { return N; };
};

class RsqTracker1 {
    // Tracks only a single dt (skip)
   public:
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> pastlocs;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> xyz2sums;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> xyz4sums;
    vector<flt> r4sums;
    unsigned long skip, count;

   public:
    RsqTracker1(AtomGroup& atoms, unsigned long skip, Vec com);

    void reset(AtomGroup& atoms, Vec com);

    bool update(Box& box, AtomGroup& atoms, unsigned long t,
                Vec com);  // updates if necessary.
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> xyz2();
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> xyz4();
    vector<flt> r4();

    unsigned long get_skip() { return skip; };
    unsigned long get_count() { return count; };
};

class RsqTracker : public StateTracker {
   public:
    sptr<AtomGroup> atoms;
    vector<RsqTracker1> singles;
    unsigned long curt;
    bool usecom;

   public:
    RsqTracker(sptr<AtomGroup> atoms, vector<unsigned long> ns,
               bool usecom = true);

    void reset();
    void update(Box& box);

    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > xyz2();
    vector<vector<flt> > r2();
    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > xyz4();
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
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> pastlocs;
    vector<vector<barray<cmplx, NDIM> > >
        ISFsums;  // (number of ks x number of particles x number of dimensions)
    vector<flt> ks;
    unsigned long skip, count;

   public:
    ISFTracker1(AtomGroup& atoms, unsigned long skip, vector<flt> ks, Vec com);

    void reset(AtomGroup& atoms, Vec com);

    bool update(Box& box, AtomGroup& atoms, unsigned long t,
                Vec com);  // updates if necessary.
    vector<vector<cmplx> > ISFs();
    vector<vector<barray<cmplx, NDIM> > > ISFxyz();

    unsigned long get_skip() { return skip; };
    unsigned long get_count() { return count; };
};

class ISFTracker : public StateTracker {
   public:
    sptr<AtomGroup> atoms;
    vector<ISFTracker1> singles;
    unsigned long curt;
    bool usecom;

   public:
    ISFTracker(sptr<AtomGroup> atoms, vector<flt> ks, vector<unsigned long> ns,
               bool usecom = false);

    void reset();
    void update(Box& box);

    vector<vector<vector<barray<cmplx, NDIM> > > > ISFxyz();
    vector<vector<vector<cmplx> > > ISFs();
    vector<flt> counts();
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Get the "smoothed" version of each particles location, "smoothed" over
// windows of a given size

class SmoothLocs : public StateTracker {
   public:
    sptr<AtomGroup> atoms;
    uint smoothn;  // number of skipns to "smooth" over
    uint skipn;    // number of dts to skip
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> curlocs;
    uint numincur;
    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > locs;
    unsigned long curt;
    bool usecom;

   public:
    SmoothLocs(sptr<AtomGroup> atoms, Box& box, uint smoothn, uint skipn = 1,
               bool usecom = false);

    void reset();
    void update(Box& box);

    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > smooth_locs() {
        return locs;
    };
};

////////////////////////////////////////////////////////////////////////////////////////////////////
class RDiffs : public StateTracker {
    // Tracks only a single dt (skip)
   public:
    sptr<AtomGroup> atoms;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> pastlocs;
    vector<vector<flt> > dists;
    unsigned long skip;
    unsigned long curt;
    bool usecom;

   public:
    RDiffs(sptr<AtomGroup> atoms, unsigned long skip, bool usecom = false);

    void reset();

    void update(Box& box);  // updates if necessary.
    vector<vector<flt> > rdiffs() { return dists; };

    unsigned long get_skip() { return skip; };
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// For comparing two jammed structures

/* We have two packings, A and B, and want to know the sequence {A1, A2, A3...}
 * such that particle A1 of packing 1 matches particle 1 of packing B.
 * A JammingList is a partial list; it has a list {A1 .. An}, with n / N
 * particles assigned, with a total distance² of distance_squared.
*/
class JammingList {
   public:
    vector<uint> assigned;
    flt distance_squared;

    JammingList() : assigned(), distance_squared(0){};
    JammingList(const JammingList& other)
        : assigned(other.assigned), distance_squared(other.distance_squared){};
    JammingList(const JammingList& other, uint expand, flt addeddist)
        : assigned(other.size() + 1, 0),
          distance_squared(other.distance_squared + addeddist) {
        for (uint i = 0; i < other.size(); i++) {
            assigned[i] = other.assigned[i];
        }
        assigned[assigned.size() - 1] = expand;
    }
    inline uint size() const { return (uint)assigned.size(); };

    bool operator<(const JammingList& other) const;
};

//! A class for determining if two packings are the same.
//! Given a list of particle positions A and a second list of particle positions
//! B, after repeated calls to expand, this will (after repeated calls to
//! expand) return a JammingList that
class JammingTree {
   private:
    sptr<Box> box;
    list<JammingList> jlists;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> A;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> B;

   public:
    JammingTree(sptr<Box> box, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B)
        : box(box), jlists(), A(A), B(B) {
        jlists.push_back(JammingList());
        assert(A.size() <= B.size());
    };

    bool expand() {
        JammingList curjlist = jlists.front();
        vector<uint>& curlist = curjlist.assigned;
        if (curlist.size() >= (uint)A.rows()) {
            //~ cout << "List already too big\n";
            return false;
        }

        list<JammingList> newlists = list<JammingList>();
        for (uint i = 0; i < B.size(); i++) {
            vector<uint>::iterator found =
                find(curlist.begin(), curlist.end(), i);
            // if (find(curlist.begin(), curlist.end(), i) != curlist.end()){
            if (found != curlist.end()) {
                //~ cout << "Found " << i << "\n";
                // cout << found << '\n';
                continue;
            }
            flt newdist =
                box->diff(A.row(curlist.size()), B.row(i)).squaredNorm();
            JammingList newjlist = JammingList(curjlist, i, newdist);
            newlists.push_back(newjlist);
            //~ cout << "Made " << i << "\n";
        }

        if (newlists.empty()) {
            //~ cout << "No lists made\n";
            return false;
        }
        //~ cout << "Have " << newlists.size() << "\n";
        newlists.sort();
        //~ cout << "Sorted.\n";
        jlists.pop_front();
        //~ cout << "Popped.\n";
        jlists.merge(newlists);
        //~ cout << "Merged to size " << jlists.size() << "best dist now " <<
        // jlists.front().distance_squared << "\n";
        return true;
    }
    bool expand(uint n) {
        bool retval = false;
        for (uint i = 0; i < n; i++) {
            retval = expand();
        }
        return retval;
    }
    list<JammingList>& my_list() { return jlists; };
    list<JammingList> copy_list() { return jlists; };

    JammingList current_best() {
        JammingList j = JammingList(jlists.front());
        //~ cout << "Best size: " << j.size() << " dist: " <<
        // j.distance_squared;
        //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " <<
        // j.assigned[j.size()-1] << "]";
        //~ cout << '\n';
        return j;
        // return JammingList(jlists.front());
    };
    uint size() { return (uint)jlists.size(); };
};

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
class JammingTreeRot {
   protected:
    sptr<Box> box;
    list<JammingListRot> jlists;
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> A;
    vector<Eigen::Matrix<flt, Eigen::Dynamic, NDIM> > Bs;

   public:
    // make all 8 possible rotations / flips
    // then subtract off all possible COMVs
    JammingTreeRot(sptr<Box> box, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                   Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B,
                   bool use_rotations = true, bool use_inversions = true);
    flt distance(JammingListRot& jlist);
    list<JammingListRot> expand(JammingListRot curjlist);

    virtual bool expand();

    bool expand(uint n) {
        bool retval = false;
        for (uint i = 0; i < n; i++) {
            retval = expand();
            if (!retval) break;
        }
        return retval;
    }
    bool expand_to(flt maxdistsq) {
        bool retval = true;
        while (
            (maxdistsq <= 0 or jlists.front().distance_squared < maxdistsq) and
            retval) {
            retval = expand();
        };
        return retval;
    }
    static Vec straight_diff(Box& bx,
                             Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                             Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B);
    static flt straight_distsq(Box& bx,
                               Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                               Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B);

    list<JammingListRot>& my_list() { return jlists; };
    list<JammingListRot> copy_list() { return jlists; };
    list<JammingListRot> copy_list(uint n) {
        list<JammingListRot>::iterator last = jlists.begin();
        advance(last, n);
        return list<JammingListRot>(jlists.begin(), last);
    };

    JammingListRot current_best() {
        if (jlists.empty()) {
            JammingListRot bad_list = JammingListRot();
            bad_list.distance_squared = -1;
            return bad_list;
        }
        JammingListRot j = JammingListRot(jlists.front());
        //~ cout << "Best size: " << j.size() << " dist: " <<
        // j.distance_squared;
        //~ if(j.size() > 0) cout << " Elements: [" << j.assigned[0] << ", " <<
        // j.assigned[j.size()-1] << "]";
        //~ cout << '\n';
        return j;
        // return JammingList(jlists.front());
    };

    // JammingListRot operator[](uint i){
    //    assert(i < jlists.size());
    //    return JammingListRot(jlists[i]);
    //};

    uint size() { return (uint)jlists.size(); };

    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locations_B(JammingListRot jlist);
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locations_B() {
        return locations_B(current_best());
    };
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locations_A(JammingListRot jlist);
    Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locations_A() {
        return locations_A(current_best());
    };
    virtual ~JammingTreeRot(){};
};

class JammingTreeBD : public JammingTreeRot {
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
    uint cutoff1, cutoff2;

   public:
    JammingTreeBD(sptr<Box> box, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                  Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B, uint cutoff,
                  bool use_rotations = true, bool use_inversions = true)
        : JammingTreeRot(box, A, B, use_rotations, use_inversions),
          cutoff1(cutoff),
          cutoff2(cutoff){};
    JammingTreeBD(sptr<Box> box, Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& A,
                  Eigen::Matrix<flt, Eigen::Dynamic, NDIM>& B, uint cutoffA,
                  uint cutoffB, bool use_rotations = true,
                  bool use_inversions = true);  // :
    // JammingTree2(box, A, B), cutoff1(cutoffA), cutoff2(cutoffB){};

    list<JammingListRot> expand(JammingListRot curjlist);
    bool expand();
    bool expand(uint n) { return JammingTreeRot::expand(n); };
};

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
 *      # for each pair, if box.diff(a1.x, a2.x) < (sigma1 < sigma2)/2, add it
 * to the list
 * }
 *
 * struct Path {
 *      Vec<uint> distance          # Euclidean distance from first node to
 * last. Note this is not
 *                                  #   the same as box.diff(last - first),
 * because it might go
 *                                  #   a longer way around the box
 *      vector<uint> nodes          # nodes visited
 * }
 *
 * Maybe implement <, ==, > by using length of nodes first and what the nodes
 * are second
 *
 * Implement ... circular_from(uint n, check_all=true):
 *  -   Do a breadth-first search starting from n, including only nodes > n.
 *  -   While searching, build map<node, path>, which is the Euclidean distance
 * from the start node
 *      to the current node via the path traveled
 *  -   If you get to an already visited node, compare the current Vec path to
 * the previous one. If
 *      the same:
 *      -   either replace the old path (if new_path < old_path) or don't
 *      -   Do not follow connections
 *      If different, we've found a percolation:
 *      -   Mark percolations with "map<uint, path> cycles"; the key is the
 * dimension, and
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

    bool operator!=(const CNode& other) const { return n != other.n; };
    bool operator==(const CNode& other) const { return n == other.n; };
    bool operator<(const CNode& other) const { return n < other.n; };
    bool operator>(const CNode& other) const { return n > other.n; };
};

class CNodePath {
   public:
    Vec distance;
    vector<CNode> nodes;

   public:
    CNodePath(){};
    CNodePath(CNode node) { nodes.push_back(node); };
    CNodePath(CNodePath other, CNode node, OriginBox& box)
        : distance(other.distance), nodes(other.nodes) {
        add(node, box);
    };
    void add(CNode node, OriginBox& box);
    uint size() { return (uint)nodes.size(); };
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

    barray<bool, NDIM> nonzero(
        Vec diff_vec);  // is it nonzero. Specifically, is any dimension >  L/2?
    CNodePath make_cycle(CNodePath forward, CNodePath backward);

    map<uint, CNodePath> circular_from(CNode node, set<uint>& visited,
                                       bool check_all);

   public:
    Connectivity(sptr<OriginBox> box) : box(box){};
    void add_edge(CNode node1, CNode node2);

    // assumes diameters are additive
    // Note that "diameter" should generally be the attractive diameter
    void add(Eigen::Matrix<flt, Eigen::Dynamic, NDIM> locs,
             vector<flt> diameters);

    map<uint, CNodePath> find_percolation(bool check_all_dims = true);
};

#endif
