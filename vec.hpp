/* An implementation of a number of vector classes.
 * 
 * Nvector is an n-dimensional vector, allowing addition, subtraction,
 * multiplication by a scalar. You can have Nvectors of Nvectors.
 * Numvector requires numbers as elements, and allows for calculation of
 * a dot product, magnitude, and normalizing.
 * Vector is a 3D vector, and adds a cross-product, as well as accessors
 * for x, y, and z components by name.
 */

#ifndef VEC_H
#define VEC_H

#include <ostream>
#include <cmath>
using namespace std;

template <class T, unsigned int N>
class Nvector {
    protected:
        T vals[N];
    public:
        typedef T* iterator;
        Nvector();
        Nvector(const Nvector &rhs);
        template <class U> Nvector(const Nvector<U, N> &rhs);
        Nvector(const T locs[N]);
        const T& get(const unsigned int n) const {return vals[n];}
        void set(const unsigned int n, const T a){vals[n]=a;}
        unsigned int len() const {return N;};
        
        Nvector& operator+=(const Nvector &rhs);
        Nvector& operator-=(const Nvector &rhs);
        template <class U> Nvector& operator*=(const U rhs);
        template <class U> Nvector& operator/=(const U rhs);
        Nvector operator-() const;
        Nvector operator+(const Nvector &rhs) const {return Nvector(*this) += rhs;}
        Nvector operator-(const Nvector &rhs) const {return Nvector(*this) -= rhs;}
        T& operator[](const unsigned int i){return vals[i];};
        template <class U> Nvector operator*(const U rhs) const {return Nvector(*this) *= rhs;}
        template <class U> Nvector operator/(const U rhs) const {return Nvector(*this) /= rhs;}
        T* begin(){return vals;};
        T* end(){return vals + N;};
        ~Nvector(){};
        
        template <class U, unsigned int M>
        friend ostream& operator<<(ostream& out, const Nvector<U, M> v);
};

//~ typedef typename std::iterator_traits<>::value_type cont;
//~ typedef typename cont::const_iterator const_iterator;

template <class T, unsigned int N>
class Numvector : public Nvector<T, N> {
    public:
        Numvector() {for(unsigned int i=0; i<N; i++) this->vals[i]=0;}
        Numvector(const Nvector<T, N> &rhs) {
                    for(unsigned int i=0; i<N; i++) this->vals[i]=rhs.get(i);}
        Numvector(const T rhs[N]) {
                    for(unsigned int i=0; i<N; i++) this->vals[i]=rhs[i];}
        T dot (const Numvector &other) const;
        inline T sq() const {return dot(*this);};
        inline T mag() const {return sqrt(sq());};
        inline T distance(const Numvector &rhs) const;
        Numvector perp(const Numvector &other) const;
        //returns the component of this perpendicular to other
        
        void normalize();
        Numvector norm() const;
        ~Numvector(){};
        
        template <class U, unsigned int M>
        friend ostream& operator<<(ostream& out, const Numvector<U, M> v);
};

template <class T>
class Vector : public Numvector<T, 3> {
    public:
        Vector() {setx(0); sety(0); setz(0);}
        Vector(const T a, const T b, const T c) {setx(a); sety(b); setz(c);}
        Vector(const Numvector<T, 3> rhs) {setx(rhs.get(0)); sety(rhs.get(1)); setz(rhs.get(2));}
        Vector(const Nvector<T, 3> rhs) {setx(rhs.get(0)); sety(rhs.get(1)); setz(rhs.get(2));}
        inline const T getx() const {return this->get(0);}
        inline const T gety() const {return this->get(1);}
        inline const T getz() const {return this->get(2);}
        inline void setx(const T a){this->vals[0]=a;}
        inline void sety(const T b){this->vals[1]=b;}
        inline void setz(const T c){this->vals[2]=c;}
        inline Vector operator-() const{
            return Vector(-getx(),-gety(),-getz());}
        inline Vector operator+(const Vector &rhs) const {
            return Vector(getx()+rhs.getx(),gety()+rhs.gety(),getz()+rhs.getz());}
        inline Vector operator-(const Vector &rhs) const {
            return Vector(getx()-rhs.getx(),gety()-rhs.gety(),getz()-rhs.getz());}
        template <class U> inline Vector operator*(const U rhs) const {
            return Vector(getx()*rhs,gety()*rhs,getz()*rhs);}
        template <class U> inline Vector operator/(const U rhs) const {
            return Vector(getx()/rhs,gety()/rhs,getz()/rhs);}
        Vector cross (const Vector &rhs) const;
        inline Vector norm() const {return Vector(Numvector<T,3>::norm());};
        template <class U> Vector& operator*=(const U rhs);
        template <class U> Vector& operator/=(const U rhs);
        ~Vector(){};
        
        template <class U>
        friend ostream& operator<<(ostream& out, const Vector<U> v);

};

template <class T, unsigned int N>
Nvector<T, N>::Nvector(const Nvector<T, N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = rhs.get(i);
    }
}

template <class T, unsigned int N>
Nvector<T, N>::Nvector() {}

template <class T, unsigned int N>
Nvector<T, N>::Nvector(const T locs[N]) {
    for(unsigned int i=0; i < N; i++) vals[i] = locs[i];
}

template <class T, unsigned int N> template<class U>
Nvector<T, N>::Nvector(const Nvector<U,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = T(rhs.get(i));
    }
}

template <class T, unsigned int N>
Nvector<T, N> Nvector<T, N>::operator-() const {
    Nvector<T, N> newvec = Nvector(*this);
    newvec *= -1;
    return newvec;
}

template <class T, unsigned int N>
Nvector<T, N>& Nvector<T, N>::operator-=(const Nvector<T,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] -= rhs.get(i);
    }
    return *this;
}

template <class T, unsigned int N>
Nvector<T, N>& Nvector<T, N>::operator+=(const Nvector<T,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] += rhs.get(i);
    }
    return *this;
}

template <class T, unsigned int N> template <class U>
Nvector<T,N>& Nvector<T,N>::operator*=(const U rhs) {
    for(unsigned int i=0; i<N; i++){
        this->vals[i] = this->vals[i] * rhs;
        //set(i, T(this->get(i) * rhs));
    }
    return *this;
}

template <class T, unsigned int N> template <class U>
Nvector<T,N>& Nvector<T,N>::operator/=(const U rhs) {
    for(unsigned int i=0; i<N; i++){
        vals[i] = T(this->vals[i] / rhs);
    }
    return *this;
}

template <class U, unsigned int M>
ostream& operator<< (ostream& out, const Nvector<U,M> v){
    out << "[" << v.get(0);
    for(unsigned int i = 1; i < M; i++)
        out << ',' << v.get(i);
    return out << ']';
}

template <class T, unsigned int N>
T Numvector<T,N>::dot(const Numvector<T, N> &other) const{
    T m = 0;
    for(unsigned int i=0; i<N; i++){
        m += this->get(i) * other.get(i);
    }
    return m;
}

template <class T, unsigned int N>
Numvector<T,N> Numvector<T,N>::perp(const Numvector<T,N>& other) const {
    //~ Numvector<T,N> parallel = other * (this->dot(other)) / other.dot(other);
    return (*this) - other * ((this->dot(other)) / other.dot(other));
    //~ return other - parallel;
}

template <class T, unsigned int N>
void Numvector<T,N>::normalize(){
    T m = mag();
    *this /= m;
}

template <class T, unsigned int N>
Numvector<T,N> Numvector<T,N>::norm() const{
    T m = mag();
    return *this / m;
}


template <class T, unsigned int N>
T Numvector<T,N>::distance(const Numvector<T,N> &rhs) const{
    T sum = 0;
    for(unsigned int i=0; i<N; i++){
        sum += pow(this->get(i) - rhs.get(i), 2);
    }
    return sqrt(sum);
    // return Numvector<T,N>(*this - rhs).mag();
}

template <class U, unsigned int M>
ostream& operator<< (ostream& out, const Numvector<U,M> v){
    out << "[" << v.get(0);
    for(int i = 1; i < M; i++)
        out << ',' << v.get(i);
    return out << ']';
}

template <class T>
Vector<T> Vector<T>::cross(const Vector<T> &rhs) const {
    T newx = gety()*rhs.getz() - rhs.gety()*getz();
    T newy = getz()*rhs.getx() - rhs.getz()*getx();
    T newz = getx()*rhs.gety() - rhs.getx()*gety();
    return Vector<T>(newx, newy, newz);
}

template <class T> template <class U>
Vector<T>& Vector<T>::operator*=(const U rhs){
    this->vals[0] *= rhs;
    this->vals[1] *= rhs;
    this->vals[2] *= rhs;
    return *this;
}

template <class T> template <class U>
Vector<T>& Vector<T>::operator/=(const U rhs){
    this->vals[0] /= rhs;
    this->vals[1] /= rhs;
    this->vals[2] /= rhs;
    return *this;
}

template <class U>
ostream& operator<< (ostream& out, const Vector<U> v){
    out << "(" << v.get(0);
    for(int i = 1; i < 3; i++)
        out << ',' << v.get(i);
    return out << ')';
}

#endif
