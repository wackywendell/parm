#include "collection.hpp"
//#include <boost/python.hpp>
//#include <boost/python/detail/wrap_python.hpp>
//#include <boost/python/module.hpp>

#include <boost/python/operators.hpp>
#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost/python/pure_virtual.hpp>
#include <boost/python/list.hpp>
//~ #include <boost/python/args.hpp>
//~ #include <boost/operators.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <iostream>

using namespace boost::python;

template <class T> 
vector<T> tovec(list l){
    ssize_t n = len(l);
    vector<T> v = vector<T>(n);
  for(ssize_t i=0;i<n;i++) 
    v[i] = extract<T>(l[i]);
  return v;
}

template <class T, class P> 
vector<P> tovec(list l){
    ssize_t n = len(l);
    vector<P> v = vector<P>(n);
    for(ssize_t i=0;i<n;i++){
        v[i] = (P)ptr<P>(extract<T>(l[i]));
    }
    return v;
}

struct atomwrap {
    private:
        atomvec vec;
    public:
        atomwrap(const flt L, cuint N, vector<flt> ms):vec(L,N,ms){};
        atomwrap(const flt L, cuint N, list ms):vec(L,N,tovec<flt>(ms)){};
        //~ atomwrap(const flt L, cuint N, list ms):vec(L,N,vector<flt>(N)){
            //~ for(uint i=0; i<N; i++) vec.setmass(i,extract<flt>(ms[i]));
        //~ };
        
        operator atomvec*(){return &vec;};
        operator atomgroup*(){return &vec;};
        atom get(cuint n){return vec[n];};
        void set(cuint n, atom a){vec[n]=a;};
        uint N() const{return vec.N();};
        Vec com() const{return vec.com();};//center of mass
        //~ Vec comvel() const; //center of mass velocity
        //~ virtual flt getmass(const unsigned int n) const = 0;
        //~ virtual Vec diff(cuint n, cuint m) const = 0;
        //~ virtual Vec diff(cuint n, const Vec r) const = 0;
        //~ flt mass() const;
        //~ flt kinetic(const Vec &originvelocity=Vec(0,0,0)) const;
        //~ Vec momentum() const;
        //~ Vec angmomentum(const Vec &loc) const;
        //~ flt mominertia(const Vec &loc, const Vec &axis) const;
        //~ void resetForces();
        //~ void setAccel();
};

class bondLength : public intraMolNNPair {
    private:
        spring s;
    public:
        bondLength(list gvec, const flt k, const flt x0)
        :s(k,x0),
        intraMolNNPair(tovec<atomwrap, atomgroup*>(gvec), &s){};
};

BOOST_PYTHON_MODULE(simulation){   
    class_<Vec>("Vec")
        .def(init<flt, flt, flt>())
        .add_property("x", &Vec::getx, &Vec::setx)
        .add_property("y", &Vec::gety, &Vec::sety)
        .add_property("z", &Vec::getz, &Vec::setz)
        .def("cross", &Vec::cross)
        .def("dot", &Vec::dot)
        .def("sq", &Vec::sq)
        .def("mag", &Vec::mag)
        .def("norm", &Vec::norm)
        .def("normalize", &Vec::normalize)
        .def(self + self)
        .def(self - self)
        .def(-self)
        .def(self += self)
        .def(self -= self)
        .def(self * flt())
        .def(self / flt())
        .def(self *= flt())
        .def(self /= flt())
        
        //~ .def(flt() * self) // not defined
        //~ .def(flt() / self)
        
        .def(str(self))
        .def(repr(self))
    ;
    
    class_<atom>("atom")
        .def_readwrite("x", &atom::x)
        .def_readwrite("v", &atom::v)
        .def_readwrite("a", &atom::a)
        .def_readwrite("f", &atom::f)
    ;
    
    //~ struct atomgroupWrap : atomgroup, wrapper<atomgroup>{
        //~ atom& operator[](cuint n){return this->get_override("[]")();}
        //~ atom& operator[](cuint n) const{return this->get_override("[]")();}
        //~ flt getmass(cuint n) const{return this->get_override("getmass")();}
        //~ Vec diff(cuint n, cuint m) const{return this->get_override("diff")();}
        //~ Vec diff(cuint n, const Vec r) const{return this->get_override("diff")();}
    //~ };
    //~ 
    //~ atom& (atomgroupWrap::*op1)(cuint) = &atomgroupWrap::operator[];
    //~ atom& (atomgroupWrap::*op2)(cuint) const = &atomgroupWrap::operator[];
    //~ Vec (atomgroupWrap::*d1)(cuint, cuint) const = &atomgroupWrap::diff;
    //~ Vec (atomgroupWrap::*d2)(cuint, Vec) const = &atomgroupWrap::diff;
    

    //~ class_<atomgroupWrap, boost::noncopyable>("atomgroup",init<const flt, cuint, flt*>())
        //~ .def("[]", pure_virtual(op1));
        //~ .def("getmass", pure_virtual(getmass))
        //~ .def("diff", pure_virtual(&atomgroup::diff);
    ;
    
    class_<vector<flt> >("fvector")
        .def(vector_indexing_suite<std::vector<flt> >())
    ;
    
    class_<atomwrap>("atomvec", init<const flt, cuint, vector<flt> >())
        .def(init<const flt, cuint, list >())
        .def("com", &atomwrap::com)
        .def("N", &atomwrap::N)
        .def("__getitem__", &atomwrap::get)
        .def("__setitem__", &atomwrap::set)
        .def("__len__", &atomwrap::N)
    ;
    
    class_<LJforce>("LJforce", init<const flt, const flt>())
        //~ .def("energy",&LJforce::energy)
        //~ .def("forces",&LJforce::forces)
    ;
    
    class_<LJcutoff>("LJcutoff", init<const flt, const flt, const flt>());
    
    class_<spring>("spring", init<const flt, const flt>());
    class_<bondangle>("bondangle", init<const flt, const flt, bool>())
        .def(init<const flt,const flt>());
    class_<dihedral>("dihedral", init<const vector<flt> >());

    class_<bondLength>("bondLength", init<list, const flt, const flt>())
        .def("energy",&bondLength::energy)
        .def("forces",&bondLength::setForces)
    ;
}
