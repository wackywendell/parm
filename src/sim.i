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
%module sim3d
#endif
#endif
#pragma SWIG nowarn=302,321,389
%feature("autodoc", "1");

%include pyabc.i
%include typemaps.i
%include exception.i
%include carrays.i
%include complex.i
%include std_vector.i
%include std_list.i
%include std_set.i
%include boost_shared_ptr.i
%include array.i
%include std_map.i

%apply double { long double } 

%include <pybuffer.i>
%pybuffer_mutable_binary(double *buffer, size_t sizet)

%shared_ptr(Box)
%shared_ptr(OriginBox)
%shared_ptr(InfiniteBox)
%shared_ptr(LeesEdwardsBox)
%shared_ptr(SCBox)
%ignore AtomIter;
%shared_ptr(AtomGroup)
%shared_ptr(SubGroup)
%shared_ptr(StateTracker)
%shared_ptr(Interaction)
%shared_ptr(Constraint)
%shared_ptr(AtomVec)
%shared_ptr(MetaGroup)
%shared_ptr(NeighborList)
%shared_ptr(ContactTracker)
%shared_ptr(EnergyTracker)
%shared_ptr(RsqTracker)
%shared_ptr(ISFTracker)
%shared_ptr(SmoothLocs)
%shared_ptr(RDiffs)
%shared_ptr(SCAtomVec)
%shared_ptr(CoordConstraint)
%shared_ptr(RigidConstraint)
%shared_ptr(RandomForce)
%shared_ptr(FixedForce)
%shared_ptr(FixedForceRegion)
%shared_ptr(FixedSpring)
%shared_ptr(SoftWall)
%shared_ptr(SoftWallCylinder)
%shared_ptr(WalledBox2D)
%shared_ptr(COMSpring)
%shared_ptr(BondPairs)
%shared_ptr(AngleTriples)
%shared_ptr(Dihedrals)
%shared_ptr(InteractionPairsX)
%shared_ptr(LJsimple)
%shared_ptr(Charges)
%shared_ptr(SCSpringList)
%shared_ptr(CoordCOMConstraint)
%shared_ptr(RelativeConstraint)
%shared_ptr(DistConstraint)
%shared_ptr(LinearConstraint)
%shared_ptr(NPHGaussianConstraint)
%shared_ptr(SimpleListed< EpsSigExpAtom,RepulsionPair >)
%shared_ptr(NListed< EpsSigAtom,LJRepulsePair >)
%shared_ptr(NListed< EpsSigAtom,LJRepulsePair >)
%shared_ptr(NListed< EpsSigCutAtom,LJAttractCutPair >)
%shared_ptr(NListed< IEpsSigCutAtom,LJAttractCutPair >)
%shared_ptr(NListed< IEpsISigCutAtom,LJAttractCutPair >)
%shared_ptr(NListed< IEpsISigCutAtom,LennardJonesCutPair >)
%shared_ptr(NListed< IEpsSigCutAtom,LJAttractRepulsePair >)
%shared_ptr(NListed< IEpsRepsSigCutAtom,LJAttractFixedRepulsePair >)
%shared_ptr(NListed< LJDoubleAtom,LJDoublePair >)
%shared_ptr(NListed< EisMclachlanAtom,EisMclachlanPair >)
%shared_ptr(NListed< IEpsRepsSigExpCutAtom,LJishPair >)
%shared_ptr(NListed< LoisOhernAtom,LoisOhernPair >)
%shared_ptr(NListed< LoisOhernAtom,LoisOhernPairMinCLs >)
%shared_ptr(NListed< LoisLinAtom,LoisLinPair >)
%shared_ptr(NListed< LoisLinAtom,LoisLinPairMin >)
%shared_ptr(NListed< EpsSigExpAtom,RepulsionPair >)
%shared_ptr(NListed< IEpsISigExpAtom,RepulsionPair >)
%shared_ptr(NListed< EpsSigExpDragAtom,RepulsionDragPair >)
%shared_ptr(SCBoxed< EpsSigExpAtom,RepulsionPair >)
%shared_ptr(NListedVirial< EpsSigExpAtom,RepulsionPair >)

%{
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
#include <Python.h>

// Part of the numpy initialization sequence.
// See http://wiki.scipy.org/Cookbook/SWIG_NumPy_examples#head-6c11cb03512f8fd5bfa20f8d8c6f69e7cc1ce494
#define SWIG_FILE_WITH_INIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static int myErr = 0;
%}

// Part of the numpy initialization sequence.
%init %{
    import_array()
%}

%exception {
    try {
        $action
    } catch(std::range_error &e) {
        SWIG_exception(SWIG_ValueError, e.what());
    //~ } catch(DivisionByZero) {
    //~     SWIG_exception(SWIG_DivisionByZero, "Division by zero");
    //~ } catch(OutOfMemory) {
    //~     SWIG_exception(SWIG_MemoryError, "Out of memory");
    } catch(std::invalid_argument &e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch(std::overflow_error &e) {
        SWIG_exception(SWIG_OverflowError, e.what());
    } catch(std::exception &e) {
        SWIG_exception(SWIG_RuntimeError,e.what());
    }catch(...) {
        SWIG_exception(SWIG_RuntimeError,"Unknown exception");
    }
}

// Convert Vec2, Vec3 output into numpy arrays
// Taken from
// http://stackoverflow.com/questions/24375198/error-wrapping-eigen-c-with-python-using-swig
%typemap(out) Vec2 { 
    npy_intp dims[1] = {2}; 
    PyObject* array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double* data = ((double *)PyArray_DATA((PyArrayObject *) array)); 
    data[0] = (double) $1(0);
    data[1] = (double) $1(1);
    $result = array; 
};

// Written myself, using http://www.swig.org/Doc3.0/Typemaps.html#Typemaps
%typemap(out) Vec3 { 
    npy_intp dims[1] = {3}; 
    PyObject* array = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    double* data = ((double *)PyArray_DATA((PyArrayObject *) array)); 
    data[0] = (double) $1(0);
    data[1] = (double) $1(1);
    data[2] = (double) $1(2);
    $result = array; 
};

%typemap(out) Matrix { 
    npy_intp dims[2] = {NDIM,NDIM}; 
    PyObject* array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double* data = ((double *)PyArray_DATA((PyArrayObject *) array));
    for(uint i=0; i<NDIM; i++){
        for(uint j=0; j<NDIM; j++){
            data[i*NDIM+j] = (double) $1(i, j);
        }
    }
    $result = array; 
};

// from http://stackoverflow.com/questions/24375198/error-wrapping-eigen-c-with-python-using-swig
%typemap(in) Eigen::Matrix<flt, Eigen::Dynamic, 3> & (Eigen::Matrix<flt, Eigen::Dynamic, 3> temp) {
    if (!PySequence_Check($input)) {
        PyErr_SetString(PyExc_TypeError,"expected a sequence.");
        return NULL;
    }
    int len = PySequence_Length($input);
    if(len == 0){
        PyErr_SetString(PyExc_TypeError,"expected a sequence of length greather than 0.");
        return NULL;
    } else if (len == -1){
        PyErr_SetString(PyExc_TypeError, "Failure converting.");
        return NULL;
    }
    
    temp = Eigen::Matrix<flt, Eigen::Dynamic, 3>(len, 3);
    for(int i=0; i < len; i++){
        PyObject *obj_i = PySequence_GetItem($input, i);
        if (obj_i == NULL) {
            PyErr_SetString(PyExc_TypeError, "Failure parsing sequence.");
            return NULL;
        }
        if (!PySequence_Check(obj_i)) {
            PyErr_SetString(PyExc_TypeError, "expected a sequence of sequences.");
            return NULL;
        }
        
        PyObject* tup = PySequence_Tuple(obj_i);
        if (!PyArg_ParseTuple(tup,"ddd", 
                &temp(i, 0), &temp(i, 1), &temp(i, 2))){
            PyErr_SetString(PyExc_TypeError,"inner sequences must have 3 doubles.");
            return NULL;
        }
    }
    
    $1 = &temp;
};

%typemap(out) Eigen::Matrix<flt, Eigen::Dynamic, 3> { 
    unsigned int rows = ($1.rows());
    npy_intp dims[2] = {rows,NDIM};
    PyObject* array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    double* data = ((double *)PyArray_DATA((PyArrayObject *) array));
    for(uint i=0; i<rows; i++){
        for(uint j=0; j<NDIM; j++){
            data[i*NDIM+j] = (double) $1(i, j);
        }
    }
    $result = array; 
};

%typemap(out) vector<Eigen::Matrix<flt, Eigen::Dynamic, 3> > { 
    unsigned int vlength = $1.size();
    PyObject* pylist = PyList_New(vlength);
    if(pylist == NULL){
        PyErr_SetString(PyExc_TypeError, "Failure creating list.");
        return NULL;
    }
    for(unsigned int i=0; i<vlength; i++){
        Eigen::Matrix<flt, Eigen::Dynamic, 3> &m = $1.at(i);
        unsigned int rows = m.rows();
        npy_intp dims[2] = {rows,NDIM};
        PyObject* array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        double* data = ((double *)PyArray_DATA((PyArrayObject *) array));
        for(uint j=0; j<rows; j++){
            for(uint k=0; k<NDIM; k++){
                data[j*NDIM+k] = (double) m(j, k);
            }
        }
        if(PyList_SetItem(pylist, i, array) != 0){
            PyErr_SetString(PyExc_TypeError, "Failure setting item.");
            return NULL;
        }
    }
    $result = pylist; 
};

// Take any sequence as input for a Vec2 (or Vec3, below)
%typemap(in) Vec2 (Vec2 temp) {
  if (PySequence_Check($input)) {
    PyObject* tup = PySequence_Tuple($input);
    if (!PyArg_ParseTuple(tup,"dd", &temp(0), &temp(1))) {
      PyErr_SetString(PyExc_TypeError,"sequence must have 2 doubles.");
      return NULL;
    }
    $1 = temp;
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a sequence.");
    return NULL;
  }
};

%typemap(in) Vec2& (Vec2 temp) {
  if (PySequence_Check($input)) {
    PyObject* tup = PySequence_Tuple($input);
    if (!PyArg_ParseTuple(tup,"dd", &temp(0), &temp(1))) {
      PyErr_SetString(PyExc_TypeError,"sequence must have 2 doubles.");
      return NULL;
    }
    $1 = &temp;

  } else {
    PyErr_SetString(PyExc_TypeError,"expected a sequence.");
    return NULL;
  }
};

%typemap(in) Vec3 (Vec3 temp) {
  if (PySequence_Check($input)) {
    PyObject* tup = PySequence_Tuple($input);
    if (!PyArg_ParseTuple(tup,"ddd", &temp(0), &temp(1), &temp(2))) {
      PyErr_SetString(PyExc_TypeError,"sequence must have 3 doubles.");
      return NULL;
    }
    $1 = temp;
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a sequence.");
    return NULL;
  }
};

%typemap(in) Vec3& (Vec3 temp) {
  if (PySequence_Check($input)) {
    PyObject* tup = PySequence_Tuple($input);
    if (!PyArg_ParseTuple(tup,"ddd", &temp(0), &temp(1), &temp(2))) {
      PyErr_SetString(PyExc_TypeError,"sequence must have 3 doubles.");
      return NULL;
    }
    $1 = &temp;

  } else {
    PyErr_SetString(PyExc_TypeError,"expected a sequence.");
    return NULL;
  }
};

// Note that Vec2 has higher precedence than Vec3
%typecheck(SWIG_TYPECHECK_FLOAT_ARRAY) Vec2, Vec2&, Vec2* {
    $1 = (PySequence_Check($input) && (PySequence_Size($input) == 2)) ? 1 : 0;
};

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) Vec3, Vec3&, Vec3* {
    $1 = (PySequence_Check($input) && (PySequence_Size($input) == 3)) ? 1 : 0;
};

#ifdef VEC2D
%typemap(in) Vec = Vec2;
%typemap(out) Vec = Vec2;
%typemap(typecheck) Vec = Vec2;
%typemap(in) Vec& = Vec2&;
%typemap(typecheck) Vec& = Vec2&;
#else 
%typemap(in) Vec = Vec3;
%typemap(out) Vec = Vec3;
%typemap(typecheck) Vec = Vec3;
%typemap(in) Vec& = Vec3&;
%typemap(typecheck) Vec& = Vec3&;
#endif 

#ifdef VEC2D
%template(_JammingList) std::list<JammingList>;
%template(_JammingListRot) std::list<JammingListRot>;
#endif

%template(fvector) std::vector<float>;
%template(_ffvector) std::vector<std::vector<float> >;
%template(dvector) std::vector<double>;
%template(_ddvector) std::vector<std::vector<double> >;
%template(_dddvector) std::vector<std::vector<std::vector<double> > >;
%template(cvector) std::vector<std::complex<double> >;
%template(_ccvector) std::vector<std::vector<std::complex<double> > >;
%template(_cccvector) std::vector<std::vector<std::vector<std::complex<double> > > >;
%template(ldvector) std::vector<long double>;
// %template(_vvector) vector<Vec>;
// %template(_vvvector) vector<vector<Vec> >;
// %template(_vvvvector) vector<vector<vector<Vec> > >;
// %template(_vvector2) vector<Vec2>;
// %template(_vvvector2) vector<vector<Vec2> >;
// %template(_vvvvector2) vector<vector<vector<Vec2> > >;
// %template(_vvector3) std::vector<Vec3>;
// %template(_vvvector3) vector<vector<Vec3> >;
// %template(_vvvvector3) vector<vector<vector<Vec3> > >;
//%template(avector) vector<shared_ptr<Atomgroup> >;
//%template(aptrvector) vector<shared_ptr<Atom> >;
%template(ivector) std::vector<shared_ptr<Interaction> >;
%template(ifxvector) std::vector<shared_ptr<InteractionPairsX> >;
%template(tvector) std::vector<shared_ptr<StateTracker> >;
%template(constraintvector) std::vector<shared_ptr<Constraint> >;
#ifdef VEC2D
%template(wallvector) std::vector<shared_ptr<SoftWall> >;
#endif
%template(idvector) std::vector<AtomID>;
%template(idpairvector) std::vector<IDPair>;
%template(intvector) std::vector<int>;
%template(uintvector) std::vector<unsigned int>;
%template(ulongvector) std::vector<unsigned long>;
%template(_eventset) std::set<Event>;
%template(pair_uint_CNodePath) std::pair<unsigned int, CNodePath>;
%template(map_uint_CNodePath) std::map<unsigned int, CNodePath>;
%template(vector_CNode) std::vector<CNode>;
%template(pair_int_CNode) std::pair<int, vector<CNode> >;
%template(map_int_CNode) std::map<int, vector<CNode> >;

%template(_carray2) boost::array<std::complex<double>, 2>;
%template(_cavector2) std::vector<boost::array<std::complex<double>, 2> >;
%template(_ccavector2) std::vector<std::vector<boost::array<std::complex<double>, 2> > >;
%template(_cccavector2) std::vector<std::vector<std::vector<boost::array<std::complex<double>, 2> > > >;
%template(_carray3) boost::array<std::complex<double>, 3>;
%template(_cavector3) std::vector<boost::array<std::complex<double>, 3> >;
%template(_ccavector3) std::vector<std::vector<boost::array<std::complex<double>, 3> > >;
%template(_cccavector3) std::vector<std::vector<std::vector<boost::array<std::complex<double>, 3> > > >;

%extend Atom {
    %insert("python") %{
        def __getstate__(self):
            return (tuple(self.x),tuple(self.v),tuple(self.f), tuple(self.a))
        
        def __setstate__(self, lst):
            self.x, self.v, self.f, self.a = [Vec(*r) for r in lst]
        
        def __str__(self):
            if hasattr(self, 'name'):
                return "<Atom %s>" % self.name
            return "<Atom>"
        
        def __repr__(self):
            #ifdef VEC2D
            x,y = tuple(self.x)
            if hasattr(self, 'name'):
                return "<Atom %s at (%.2f,%.2f)>" % (self.name,x,y)
            return "<Atom at (%.2f,%.2f)>" % (x,y)
            #else
            x,y,z = tuple(self.x)
            if hasattr(self, 'name'):
                return "<Atom %s at (%.2f,%.2f,%.2f)>" % (self.name,x,y,z)
            return "<Atom at (%.2f,%.2f,%.2f)>" % (x,y,z)
            #endif
    %};
}

%extend AtomID {
    %insert("python") %{
        def __str__(self):
            return "<AtomID %s>" % self.n()
        
        def __repr__(self):
            #ifdef VEC2D
            x,y = tuple(self.x)
            return "<AtomID %s at (%.2f,%.2f)>" % (self.n(),x,y)
            #else
            x,y,z = tuple(self.x)
            return "<AtomID %s at (%.2f,%.2f,%.2f)>" % (self.n(),x,y,z)
            #endif
    %};
}

%extend AtomGroup {
    %insert("python") %{
        def __iter__(self):
            for i in range(self.size()):
                yield self.get_id(i)
    %};
};


%extend AtomVec {
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

%extend SCAtomVec {
    %insert("python") %{
        def __iter__(self):
            for i in range(self.size()):
                yield self[i]
        
        def all_pairs(self):
            for i in range(self.pairs()):
                yield self.pair(i)
        
        def __len__(self):
            return self.size()
        
        def __getitem__(self, obj):
            return self.get_id(obj)
    %};
};

%extend IDPair {
    %insert("python") %{
        def __iter__(self):
            return iter((self.first(), self.last()))
    %};
};

%exception NeighborList::__getitem__ {
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
%template(LJRepulse) NListed<EpsSigAtom, LJRepulsePair>; // Pure repulsive
%template(LJAttractCut) NListed<EpsSigCutAtom, LJAttractCutPair>;  // Both repulsive and attractive
%template(LJAttractICut) NListed<IEpsSigCutAtom, LJAttractCutPair>;  // Pure attractive
%template(LJAttractIICut) NListed<IEpsISigCutAtom, LJAttractCutPair>;  // Pure attractive, with indices for sigma and epsilon
%template(LJIICut) NListed<IEpsISigCutAtom, LennardJonesCutPair>;  // Both repulsive and attractive, with indices for sigma and epsilon

// Both repulsive and attractive, with indices for epsilon, sigma fixed per Atom.
// Positive epsilon -> Attractive + repulsive, negative for pure repulsive.
%template(LJAttractRepulse) NListed<IEpsSigCutAtom, LJAttractRepulsePair>;
%template(LJAttractFixedRepulse) NListed<IEpsRepsSigCutAtom, LJAttractFixedRepulsePair>;
%template(EisMclachlan) NListed<EisMclachlanAtom, EisMclachlanPair>;
%template(LJish) NListed<IEpsRepsSigExpCutAtom, LJishPair>;
%template(RepulsionSimple) SimpleListed<EpsSigExpAtom, RepulsionPair>;
%template(Repulsion) NListed<EpsSigExpAtom, RepulsionPair>;
%template(RepulsionII) NListed<IEpsISigExpAtom, RepulsionPair>;
%template(HertzianDrag) NListed<EpsSigExpDragAtom, RepulsionDragPair>;
%template(HertzianSC) SCBoxed<EpsSigExpAtom, RepulsionPair>;
%template(HertzianVirial) NListedVirial<EpsSigExpAtom, RepulsionPair>;
%template(LoisOhern) NListed<LoisOhernAtom, LoisOhernPair>;
%template(LoisLin) NListed<LoisLinAtom, LoisLinPair>;
%template(LoisLinMin) NListed<LoisLinAtom, LoisLinPairMin>;
%template(LoisOhernMin) NListed<LoisOhernAtom, LoisOhernPairMinCLs>;
//%rename(__lt__) JammingList::operator<;

//%{
//    shared_ptr<NListedVirial<EpsSigExpAtom, RepulsionPair> > Repulsion(NeighborList *neighbors){
//        NListedVirial<EpsSigExpAtom, RepulsionPair> *h = new NListedVirial<EpsSigExpAtom, RepulsionPair>(neighbors);
//        return shared_ptr<NListedVirial<EpsSigExpAtom, RepulsionPair> >(h);
//    };
//%}

%extend NeighborList {
  IDPair __getitem__(size_t i) {
    if (i >= $self->numpairs()) {
      myErr = 1;
      return IDPair(AtomID(), AtomID());
    }
    return $self->get((unsigned int) i);
  }
  
  %insert("python") %{
    def __iter__(self):
        for i in range(self.numpairs()):
            yield self[i]
  %}
};

%include "collection.hpp"
%include "collection.cpp"
