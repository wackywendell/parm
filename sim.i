#ifdef VEC2D
#ifdef LONGFLOAT
%module sim2dlong
#else
%module sim2d
#endif
#else
#define VEC3D
#ifdef LONGFLOAT
%module sim3dlong
#else
%module sim
#endif
#endif
#pragma SWIG nowarn=302,321,389
%feature("autodoc", "1");

%include pyabc.i
%include typemaps.i
%include exception.i
%include carrays.i
%include std_vector.i
%include std_list.i
%include std_set.i
%include boost_shared_ptr.i
%include "vec.hpp"

%apply double { long double } 

%include <pybuffer.i>
%pybuffer_mutable_binary(double *buffer, size_t sizet)

%shared_ptr(Box)
%shared_ptr(OriginBox)
%shared_ptr(InfiniteBox)
%shared_ptr(LeesEdwardsBox)
%shared_ptr(SCbox)
%ignore AtomIter;
%shared_ptr(atomgroup)
%shared_ptr(subgroup)
%shared_ptr(statetracker)
%shared_ptr(interaction)
%shared_ptr(constraint)
%shared_ptr(atomvec)
%shared_ptr(metagroup)
%shared_ptr(neighborlist)
%shared_ptr(ContactTracker)
%shared_ptr(EnergyTracker)
%shared_ptr(RsqTracker)
%shared_ptr(SCatomvec)
%shared_ptr(coordConstraint)
%shared_ptr(fixedForce)
%shared_ptr(fixedSpring)
%shared_ptr(SoftWall)
%shared_ptr(WalledBox2D)
%shared_ptr(COMSpring)
%shared_ptr(bondpairs)
%shared_ptr(angletriples)
%shared_ptr(dihedrals)
%shared_ptr(interactionpairsx)
%shared_ptr(LJsimple)
%shared_ptr(Charges)
%shared_ptr(SCSpringList)
%shared_ptr(coordCOMConstraint)
%shared_ptr(relativeConstraint)
%shared_ptr(distConstraint)
%shared_ptr(linearConstraint)
%shared_ptr(NPHGaussianConstraint)
%shared_ptr(SimpleListed< HertzianAtom,HertzianPair >)
%shared_ptr(NListed< LJatom,LJpair >)
%shared_ptr(NListed< LJatom,LJpair >)
%shared_ptr(NListed< LJatomcut,LJAttractPair >)
%shared_ptr(NListed< HydroAtom,HydroPair >)
%shared_ptr(NListed< LJAtomIndexed,LJFullPair >)
%shared_ptr(NListed< LJAttractRepulseAtom,LJAttractRepulsePair >)
%shared_ptr(NListed< LJAttractFixedRepulseAtom,LJAttractFixedRepulsePair >)
%shared_ptr(NListed< LJDoubleAtom,LJDoublePair >)
%shared_ptr(NListed< EisMclachlanAtom,EisMclachlanPair >)
%shared_ptr(NListed< LJishAtom,LJishPair >)
%shared_ptr(NListed< LoisOhernAtom,LoisOhernPair >)
%shared_ptr(NListed< HertzianAtom,HertzianPair >)
%shared_ptr(SCboxed< HertzianAtom,HertzianPair >)
%shared_ptr(NListedVirial< HertzianAtom,HertzianPair >)

%{
#include "vec.hpp"
#include "vecrand.hpp"
#include "vecrand.cpp"
#include "box.hpp"
#include "box.cpp"
#include "trackers.hpp"
#include "trackers.cpp"
#include "interaction.hpp"
#include "interaction.cpp"
#include "constraints.hpp"
#include "constraints.cpp"
#include "collection.hpp"
#include "collection.cpp"
static int myErr = 0;
%}

%typemap(in) bool value[3] (bool temp[3]) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != 3) {
    PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected 3 elements");
    return NULL;
  }
  for (i = 0; i < 3; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = (bool) PyInt_AsLong(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");      
      return NULL;
    }
  }
  $1 = temp;
};

%extend Vector3 {
    char* __str__() {
        static char temp[256];
        sprintf(temp,"(%g, %g, %g)", $self->getxd(),$self->getyd(),$self->getzd());
        return &temp[0];
    };
    char* __repr__() {
        static char temp[256];
        sprintf(temp,"Vec(%g, %g, %g)", $self->getxd(),$self->getyd(),$self->getzd());
        return &temp[0];
    };
    
    Vector3 __truediv__(const double n) const{
        return Vector3<T>($self->operator/(n));
    };
    
    Vector3 __mul__(const double n) const{
        return Vector3<T>($self->operator*(n));
    };
    
    %insert("python") %{
    @property
    def X(self):
        return self.getxd()
    
    @property
    def Y(self):
        return self.getyd()
    
    @property
    def Z(self):
        return self.getzd()
    %}
    #~ %insert("python") %{
        #~ def __iter__(self):
            #~ return iter((self.getxd(), self.getyd(), self.getzd()))
    #~ %};
};

%extend Vector2 {
    char* __str__() {
        static char temp[256];
        sprintf(temp,"(%g, %g)", $self->getxd(),$self->getyd());
        return &temp[0];
    };
    char* __repr__() {
        static char temp[256];
        sprintf(temp,"Vec(%g, %g)", $self->getxd(),$self->getyd());
        return &temp[0];
    };
    
    Vector2 __truediv__(const double n) const{
        return Vector2<T>($self->operator/(n));
    };
    
    Vector2 __mul__(const double n) const{
        return Vector2<T>($self->operator*(n));
    };
    
    %insert("python") %{
    @property
    def X(self):
        return self.getxd()
    
    @property
    def Y(self):
        return self.getyd()
    %}
};

%extend Nvector {
    %template() operator*<double>;
    %template() operator/<double>;
    
    Nvector __truediv__(const double n) const{
        return $self->operator/(n);
    };
    
    T __getitem__(const unsigned int n) const{
        return $self->get(n);
    };
    
    void __setitem__(const unsigned int n, const T val){
        $self->set(n, val);
    };
    
    unsigned int __len__(){ return N;};
    
        
    %insert("python") %{
        #~ def __setitem__(self, n, val):
            #~ return self.set(n, val)
        
        def __iter__(self):
            for i in range(len(self)):
                yield self.get(i)
        
        #~ def __len__(self):
            #~ return self.len()
    %};
};

%template(_Nvector3) Nvector<double, 3>;
%template(_Nvector2) Nvector<double, 2>;
%template(_Numvector3) Numvector<double, 3>;
%template(_Numvector2) Numvector<double, 2>;
%template(_Nvector3L) Nvector<long double, 3>;
%template(_Nvector2L) Nvector<long double, 2>;
%template(_Numvector3L) Numvector<long double, 3>;
%template(_Numvector2L) Numvector<long double, 2>;

#ifdef VEC2D
%template(Vec) Vector2<double>;
%template(VecL) Vector2<long double>;
namespace std {
    %template(vecptrvector) vector<Vector2<double>*>;
    %template(vecvector) vector<Vector2<double> >;
    %template(_jamminglist) list<jamminglist>;
    %template(_jamminglistrot) list<jamminglistrot>;
    
    %template(vecptrvectorL) vector<Vector2<long double>*>;
    %template(vecvectorL) vector<Vector2<long double> >;
}
#else
%template(Vec) Vector3<double>;
%template(VecL) Vector3<long double>;
%template(_MatrixBase) Nvector<Vector3<double>, 3>;
%template(Matr) Matrix<double>;
namespace std {
    %template(vecptrvector) vector<Vector3<double>*>;
    %template(vecvector) vector<Vector3<double> >;
    %template(vecptrvectorL) vector<Vector3<long double>*>;
    %template(vecvectorL) vector<Vector3<long double> >;
}
#endif
%template(Pair) Numvector<double, 2>;
%template(VecPair) Nvector<Vec, 2>;
%template(_atomarray2) Array<atom*, 2>;
%template(_atompair2) Array<atom, 2>;
%template(_idarray2) Array<atomid, 2>;
%template(_atomarray3) Array<atom*, 3>;
%template(_atomarray4) Array<atom*, 4>;


namespace std {
    %template(fvector) vector<float>;
    %template(_ffvector) vector<vector<float> >;
    %template(dvector) vector<double>;
    %template(_ddvector) vector<vector<double> >;
    %template(ldvector) vector<long double>;
    //%template(avector) vector<shared_ptr<atomgroup> >;
    //%template(aptrvector) vector<shared_ptr<atom> >;
    %template(ivector) vector<shared_ptr<interaction> >;
    %template(ifxvector) vector<shared_ptr<interactionpairsx> >;
    %template(tvector) vector<shared_ptr<statetracker> >;
    %template(constraintvector) vector<shared_ptr<constraint> >;
    #ifdef VEC2D
    %template(wallvector) vector<shared_ptr<SoftWall> >;
    #endif
    %template(idvector) vector<atomid>;
    %template(idpairvector) vector<idpair>;
    %template(_uintvector) vector<unsigned int>;
    %template(_eventset) set<event>;
}

%extend atom {
    %insert("python") %{
        def __getstate__(self):
            return (tuple(self.x),tuple(self.v),tuple(self.f), tuple(self.a))
        
        def __setstate__(self, lst):
            self.x, self.v, self.f, self.a = [Vec(*r) for r in lst]
        
        def __str__(self):
            if hasattr(self, 'name'):
                return "<atom %s>" % self.name
            return "<atom>"
        
        def __repr__(self):
            #ifdef VEC2D
            x,y = tuple(self.x)
            if hasattr(self, 'name'):
                return "<atom %s at (%.2f,%.2f)>" % (self.name,x,y)
            return "<atom at (%.2f,%.2f)>" % (x,y)
            #else
            x,y,z = tuple(self.x)
            if hasattr(self, 'name'):
                return "<atom %s at (%.2f,%.2f,%.2f)>" % (self.name,x,y,z)
            return "<atom at (%.2f,%.2f,%.2f)>" % (x,y,z)
            #endif
    %};
}

%extend atomgroup {
    %insert("python") %{
        def __iter__(self):
            for i in range(self.size()):
                yield self.get(i)
    %};
};


%extend atomvec {
    %insert("python") %{
        def __iter__(self):
            for i in range(self.size()):
                yield self[i]
        
        def __len__(self):
            return self.size()
        
        def __getitem__(self, obj):
            return self.get_id(obj)
        
        #def __setitem__(self, obj, val):
        #    return self.set(obj, val)
        
        def __getstate__(self):
            return ([self.getmass(i) for i in range(self.N())],
                        [a.__getstate__() for a in self])
        
        def __setstate__(self, lst):
            masses, atomstates = lst
            self.__init__(fvector(masses))
            for i, atomstate in enumerate(atomstates):
                #~ print "i:",i
                a = self.get(i)
                a.__setstate__(atomstate)
    %};
};

%exception neighborlist::__getitem__ {
  assert(!myErr);
  $action
  if (myErr) {
    myErr = 0; // clear flag for next time
    // You could also check the value in $result, but it's a PyObject here
    SWIG_exception(SWIG_IndexError, "Index out of bounds");
  }
};

%include "vecrand.hpp"
%include "box.hpp"
%include "trackers.hpp"
%include "interaction.hpp"
%include "constraints.hpp"
%template(LJgroup) NListed<LJatom, LJpair>; // Pure repulsive
%template(LJfull) NListed<LJatomcut, LJAttractPair>;  // Both repulsive and attractive
%template(Hydrophobicity) NListed<HydroAtom, HydroPair>;  // Pure attractive
%template(LJdetailed) NListed<LJAtomIndexed, LJFullPair>;  // Both repulsive and attractive, with indices for sigma and epsilon

// Both repulsive and attractive, with indices for epsilon, sigma fixed per atom.
// Positive epsilon -> Attractive + repulsive, negative for pure repulsive.
%template(LJAttractRepulse) NListed<LJAttractRepulseAtom, LJAttractRepulsePair>;
%template(LJAttractFixedRepulse) NListed<LJAttractFixedRepulseAtom, LJAttractFixedRepulsePair>;
%template(LJDouble) NListed<LJDoubleAtom, LJDoublePair>;
%template(EisMclachlan) NListed<EisMclachlanAtom, EisMclachlanPair>;
%template(LJish) NListed<LJishAtom, LJishPair>;
%template(HertzianSimple) SimpleListed<HertzianAtom, HertzianPair>;
%template(Hertzian) NListed<HertzianAtom, HertzianPair>;
%template(HertzianSC) SCboxed<HertzianAtom, HertzianPair>;
%template(HertzianVirial) NListedVirial<HertzianAtom, HertzianPair>;
%template(LoisOhern) NListed<LoisOhernAtom, LoisOhernPair>;
//%rename(__lt__) jamminglist::operator<;

//%{
//    shared_ptr<NListedVirial<HertzianAtom, HertzianPair> > Hertzian(neighborlist *neighbors){
//        NListedVirial<HertzianAtom, HertzianPair> *h = new NListedVirial<HertzianAtom, HertzianPair>(neighbors);
//        return shared_ptr<NListedVirial<HertzianAtom, HertzianPair> >(h);
//    };
//%}

%extend neighborlist {
  idpair __getitem__(size_t i) {
    if (i >= $self->numpairs()) {
      myErr = 1;
      return idpair(atomid(), atomid());
    }
    return $self->get((uint) i);
  }
  
  %insert("python") %{
    def __iter__(self):
        for i in range(self.numpairs()):
            yield self[i]
  %}
};

// %extend collection {
//     %insert("python") %{
//         def __init__(self, box, groups, interactions=None, trackers=None, constraints=None):
//             self.box = box
//             self.groups  = groups
//             self.interactions = interactions
//             self.trackers = trackers
//             self.constraints = constraints
//     %};
// }

%extend NListed<LJatom, LJpair> {
    %insert("python") %{
        def add_atom(self, epsilon, sigma, a):
            self.add(LJatom(epsilon, sigma, a))
    %};
}

%extend NListed<LJatomcut, LJAttractPair> {
    %insert("python") %{
        def add_atom(self, epsilon, sigma, a, cut):
            self.add(LJatomcut(epsilon, sigma, a, cut))
    %};
}

%include "collection.hpp"
%include "collection.cpp"
