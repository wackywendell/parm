/** An implementation of a number of vector classes.
 *
 * NVector is an n-dimensional vector, allowing addition, subtraction,
 * multiplication by a scalar. You can have NVectors of NVectors.
 *
 * NumVector requires numbers as elements, and allows for calculation of
 * a dot product, magnitude, and normalizing.
 *
 * Vector is a 3D vector, and adds a cross-product, as well as accessors
 * for x, y, and z components by name.
 */

#ifndef VEC_H
#define VEC_H

#include <ostream>
#include <cmath>
using namespace std;

typedef unsigned int uint;

/**
A fixed size array. This is just a wrapper with some convenience methods,
and the following is (roughly) equivalent:

    Array<flt, 4> arr;
    flt[4] arr;

@tparam T a type that can be added and subtracted.
@tparam N the number of dimensions.
*/
template <class T, unsigned int N>
class Array {
    protected:
        T vals[N];
    public:
        Array();
        Array(const Array &rhs);
        template <class U> Array(const Array<U, N> &rhs);
        Array(const T locs[N]);
        const T& get(const unsigned int n) const {return vals[n];}
        void set(const unsigned int n, const T a){vals[n]=a;}
        unsigned int len() const {return N;};

        T& operator[](const unsigned int i){return vals[i];};
        const T& operator[](const unsigned int i) const {return vals[i];};

        T* begin(){return vals;};
        T* end(){return vals + N;};
        ~Array(){};
};

/**
An N-dimensional vector, extending addition and subtraction from the type T to the
NVector class.

This is extended by NumVector for things like multiplication and division, but you can use NVector
for something like an array of NumVectors.

@tparam T a type that can be added and subtracted.
@tparam N the number of dimensions.
*/
template <class T, unsigned int N>
class NVector {
    protected:
        T vals[N];
    public:
        typedef T* iterator;
        NVector();
        NVector(const NVector &rhs);
        template <class U> NVector(const NVector<U, N> &rhs);
        NVector(const T locs[N]);
        const T& get(const unsigned int n) const {return vals[n];}
        void set(const unsigned int n, const T a){vals[n]=a;}
        unsigned int len() const {return N;};

        NVector& operator+=(const NVector &rhs);
        NVector& operator-=(const NVector &rhs);
        template <class U> NVector& operator*=(const U rhs);
        template <class U> NVector& operator/=(const U rhs);
        NVector operator-() const;
        NVector operator+(const NVector &rhs) const {return NVector(*this) += rhs;}
        NVector operator-(const NVector &rhs) const {return NVector(*this) -= rhs;}
        T& operator[](const unsigned int i){return vals[i];};
        const T& operator[](const unsigned int i) const{return vals[i];};

        //! Multiplication by a scalar
        template <class U> NVector operator*(const U rhs) const {return NVector(*this) *= rhs;}
        //! Division by a scalar
        template <class U> NVector operator/(const U rhs) const {return NVector(*this) /= rhs;}
        T* begin(){return vals;};
        T* end(){return vals + N;};
        ~NVector(){};

        template <class U, unsigned int M>
        friend ostream& operator<<(ostream& out, const NVector<U, M> v);
};

/**
An N-dimensional physics vector, extending NVector. This extends addition, subtraction, and
multiplication from the type T to the NVector class.

@tparam T a numerical type, as as float or double.
@tparam N the number of dimensions.
*/
template <class T, unsigned int N>
class NumVector : public NVector<T, N> {
    public:
        inline NumVector() {for(unsigned int i=0; i<N; i++) NVector<T,N>::vals[i]=0;}
        inline NumVector(const NVector<T, N> &rhs) {
                    for(unsigned int i=0; i<N; i++) NVector<T,N>::vals[i]=rhs.get(i);}
        inline NumVector(const T rhs[N]) {
                    for(unsigned int i=0; i<N; i++) NVector<T,N>::vals[i]=rhs[i];}

        T dot (const NumVector &other) const; //!< Inner product.
        inline T sq() const {return dot(*this);}; //!< The square of the vector, \f$\vec r^2 \f$.
        inline T mag() const {return sqrt(sq());}; //!< The magnitude of the vector, \f$ \left| \vec r \right| \f$.
        /*! The magnitude of the vector distance, \f$ \left| \vec r - \vec s \right| \f$.

        @param rhs The other vector, \f$\vec s \f$
        */
        inline T distance(const NumVector &rhs) const;
        /*! returns the component of this perpendicular to other.

        \f$ \vec r - \frac{\vec r \cdot \vec s}{\vec s^2} \vec s  \f$

        where \f$ \vec r \f$ is this vector, and \f$\vec s \f$ is `rhs`.

        @param other The other vector, \f$\vec s\f$
        */
        NumVector perpto(const NumVector &other) const;


        void normalize(); //!< Normalize in place
        NumVector norm() const; //!< Return the normalized version
        ~NumVector(){};

        template <class U, unsigned int M>
        friend ostream& operator<<(ostream& out, const NumVector<U, M> v);
};

/**
A 3D physics vector, with methods for adding, subtracting, dot product, etc.

This is aliased as Vec when compiled with NDIM=3.

@ingroup basics

@tparam T a numerical type, as as float or double.
*/
template <class T>
class Vector3 : public NumVector<T, 3> {
    public:
        Vector3() {setx(0); sety(0); setz(0);}
        Vector3(const T a, const T b, const T c) {setx(a); sety(b); setz(c);}
        Vector3(const NumVector<T, 3> rhs) {setx(rhs.get(0)); sety(rhs.get(1)); setz(rhs.get(2));}
        Vector3(const NVector<T, 3> rhs) {setx(rhs.get(0)); sety(rhs.get(1)); setz(rhs.get(2));}
        inline const T getx() const {return NVector<T,3>::get(0);}
        inline const T gety() const {return NVector<T,3>::get(1);}
        inline const T getz() const {return NVector<T,3>::get(2);}
        /// Return x as a double. Useful with some versions of Python and long doubles.
        inline double getxd() const {return double(NVector<T,3>::get(0));}
        inline double getyd() const {return double(NVector<T,3>::get(1));}
        inline double getzd() const {return double(NVector<T,3>::get(2));}
        inline void setx(const T a){NVector<T,3>::vals[0]=a;}
        inline void sety(const T b){NVector<T,3>::vals[1]=b;}
        inline void setz(const T c){NVector<T,3>::vals[2]=c;}
        /// Set x with a double. Useful with some versions of Python and long doubles.
        inline void setxd(const double a){NVector<T,3>::vals[0]=a;}
        inline void setyd(const double b){NVector<T,3>::vals[1]=b;}
        inline void setzd(const double c){NVector<T,3>::vals[2]=c;}
        inline void set(const T a, const T b, const T c){
            NVector<T,3>::vals[0]=a; NVector<T,3>::vals[1]=b; NVector<T,3>::vals[2]=c;}
        inline Vector3 operator-() const{
            return Vector3(-getx(),-gety(),-getz());}
        inline Vector3 operator+(const Vector3 &rhs) const {
            return Vector3(getx()+rhs.getx(),gety()+rhs.gety(),getz()+rhs.getz());}
        inline Vector3 operator-(const Vector3 &rhs) const {
            return Vector3(getx()-rhs.getx(),gety()-rhs.gety(),getz()-rhs.getz());}
        inline T operator*(const Vector3 &rhs) const {
            return (getx()*rhs.getx() + gety()*rhs.gety() + getz()*rhs.getz());}
        template <class U> inline Vector3 operator*(const U rhs) const {
            return Vector3(getx()*rhs,gety()*rhs,getz()*rhs);}
        template <class U> inline Vector3 operator/(const U rhs) const {
            return Vector3(getx()/rhs,gety()/rhs,getz()/rhs);}
        Vector3 cross (const Vector3 &rhs) const;///< Cross product, \f$\vec r \times \vec s \f$
        inline Vector3 norm() const {return Vector3(NumVector<T,3>::norm());};
        inline Vector3& operator-=(const Vector3 &rhs){NVector<T,3>::operator-=(rhs); return *this;};
        inline Vector3& operator+=(const Vector3 &rhs){NVector<T,3>::operator+=(rhs); return *this;};
        template <class U> Vector3& operator*=(const U rhs);
        template <class U> Vector3& operator/=(const U rhs);

        /** The angle between two vectors, assuming they start at the same point (i.e.\ the origin).

        Equal to \f$\theta = \arccos \left(\frac{\vec x_1 \cdot \vec x_2}{\left| \vec x_1 \right| \left|  \vec x_2 \right|}\right) \f$

        @param dx1 \f$\vec x_1 \f$
        @param dx2 \f$\vec x_2 \f$
        */
        static T angle(const Vector3 &dx1, const Vector3 &dx2){
            return acos(dx1.dot(dx2) / dx1.mag() / dx2.mag());
        }

        /** The angle between three points. Equivalent to `angle(x1 - x2, x3 - x2)`.

        @param x1 The first point.
        @param x2 The middle point, around which we are finding the angle.
        @param x3 The third point.
        */
        static T angle(const Vector3 &x1, const Vector3 &x2, const Vector3 &x3){
            Vector3 dx1 = x1 - x2, dx2 = x3 - x2;
            return acos(dx1.dot(dx2) / dx1.mag() / dx2.mag());
        }

        /**
        The dihedral angle between three vectors:

        \f$ \phi=\operatorname{arctan2}\left(
            \vec{r}_{1}\cdot\left(\vec{r}_{2}\times
                \vec{r}_{3}\right)\left|\vec{r}_{2}\right|,
            \left(\vec{r}_{1}\times\vec{r}_{2}\right)\cdot
                \left(\vec{r}_{2}\times\vec{r}_{3}\right)
            \right) \f$

        @param dx1 \f$\vec r_1\f$
        @param dx2 \f$\vec r_2\f$
        @param dx3 \f$\vec r_3\f$
        */
        static T dihedral(const Vector3 &dx1, const Vector3 &dx2, const Vector3 &dx3){
            return atan2(dx1.dot(dx2.cross(dx3))*dx2.mag(),
                                (dx1.cross(dx2).dot(dx2.cross(dx3))));
        }

        /** The dihedral angle between four points.

        Equivalent to `dihedral(x2 - x1, x3 - x2, x4 - x3)`. */
        static T dihedral(const Vector3 &x1, const Vector3 &x2,
                        const Vector3 &x3, const Vector3 &x4){
            Vector3 dx1 = x2 - x1, dx2 = x3 - x2, dx3 = x4 - x3;
            return atan2(dx1.dot(dx2.cross(dx3))*dx2.mag(),
                                (dx1.cross(dx2).dot(dx2.cross(dx3))));
        }
        ~Vector3(){};

        template <class U>
        friend ostream& operator<<(ostream& out, const Vector3<U> v);

};

/**
A 2D physics vector, with methods for adding, subtracting, dot product, etc.

This is aliased as Vec when compiled with NDIM=2.

@ingroup basics

@tparam T a numerical type, as as float or double.
*/
template <class T>
class Vector2 : public NumVector<T, 2> {
    public:
        Vector2() {setx(0); sety(0);}
        Vector2(const T a, const T b) {setx(a); sety(b);}
        Vector2(const NumVector<T, 2> rhs) {setx(rhs.get(0)); sety(rhs.get(1));}
        Vector2(const NVector<T, 2> rhs) {setx(rhs.get(0)); sety(rhs.get(1));}
        inline const T getx() const {return NVector<T,2>::get(0);}
        inline const T gety() const {return NVector<T,2>::get(1);}
        inline double getxd() const {return (double) NVector<T,2>::get(0);}
        inline double getyd() const {return (double) NVector<T,2>::get(1);}
        inline void setx(const T a){NVector<T,2>::vals[0]=a;}
        inline void sety(const T b){NVector<T,2>::vals[1]=b;}
        inline void setxd(const double a){NVector<T,2>::vals[0]= (T) a;}
        inline void setyd(const double b){NVector<T,2>::vals[1]= (T) b;}
        inline void set(const T a, const T b){
            NVector<T,2>::vals[0]=a; NVector<T,2>::vals[1]=b;}
        inline Vector2 operator-() const{
            return Vector2(-getx(),-gety());}
        inline Vector2 operator+(const Vector2 &rhs) const {
            return Vector2(getx()+rhs.getx(),gety()+rhs.gety());}
        inline Vector2 operator-(const Vector2 &rhs) const {
            return Vector2(getx()-rhs.getx(),gety()-rhs.gety());}
        inline T operator*(const Vector2 &rhs) const {
            return (getx()*rhs.getx() + gety()*rhs.gety());}
        template <class U> inline Vector2 operator*(const U rhs) const {
            return Vector2(getx()*rhs,gety()*rhs);}
        template <class U> inline Vector2 operator/(const U rhs) const {
            return Vector2(getx()/rhs,gety()/rhs);}
        /// The 2D cross product, returning a scalar.
        T cross (const Vector2 &rhs) const{return getx()*rhs.gety() - rhs.getx()*gety();};
        /// The 2D cross product with the "missing" third dimension, returning a vector.
        Vector2 cross (const T v) const{return Vector2(gety()*v, -getx()*v);};
        /// The vector perpendicular to this one, in the clockwise direction.
        Vector2 perp() const{return Vector2(gety(),-getx());};
        /// The normalized version of this vector.
        inline Vector2 norm() const {return Vector2(NumVector<T,2>::norm());};
        inline Vector2& operator-=(const Vector2 &rhs){NVector<T,2>::operator-=(rhs); return *this;};
        inline Vector2& operator+=(const Vector2 &rhs){NVector<T,2>::operator+=(rhs); return *this;};
        template <class U> Vector2& operator*=(const U rhs);
        template <class U> Vector2& operator/=(const U rhs);

        /// Rotate by 90 degrees counter-clockwise.
        /// \todo TODO: this direction is inconsistent with `cross` and `perp`.
        Vector2 rotate(uint i);
        inline Vector2 flip(){return Vector2(gety(), getx());};
        /** Rotate and flip, for \f$0 \le i < 8 \f$.

        For \f$0 \le i < 4 \f$, equivalent to `rotate(i)`.
        For \f$4 \le i < 8 \f$, equivalent to `flip().rotate(i % 4)`.
        For \f$i \ge 8 \f$, equivalent to `rotate_flip(i % 8)`.
        */
        inline Vector2 rotate_flip(uint i){
            if((i / 4) % 2 == 1) return flip().rotate(i%4);
            return rotate(i%4);
        };

        /**
        The inverse of `rotate_flip(i)`.
        */
        inline Vector2 rotate_flip_inv(uint i){
            Vector2 inv = rotate(4-(i%4));
            if((i / 4) % 2 == 0) return inv;
            return inv.flip();
        };
        /** The angle between two vectors, assuming they start at the same point (i.e.\ the origin).

        Equal to \f$\theta = \arccos \left(\frac{\vec x_1 \cdot \vec x_2}{\left| \vec x_1 \right| \left|  \vec x_2 \right|}\right) \f$

        @param dx1 \f$\vec x_1 \f$
        @param dx2 \f$\vec x_2 \f$
        */
        static T angle(const Vector2 &dx1, const Vector2 &dx2){
            return acos(dx1.dot(dx2) / dx1.mag() / dx2.mag());
        }

        /** The angle between three points. Equivalent to `angle(x1 - x2, x3 - x2)`.

        @param x1 The first point.
        @param x2 The middle point, around which we are finding the angle.
        @param x3 The third point.
        */
        static T angle(const Vector2 &x1, const Vector2 &x2, const Vector2 &x3){
            Vector2 dx1 = x1 - x2, dx2 = x3 - x2;
            return acos(dx1.dot(dx2) / dx1.mag() / dx2.mag());
        }
        ~Vector2(){};

        template <class U>
        friend ostream& operator<<(ostream& out, const Vector2<U> v);

};

/**
A 3x3 matrix, with methods for adding, subtracting, dot product, etc. Useful for rotations.

@ingroup basics

@tparam T a numerical type, as as float or double.
*/
template<class C>
class Matrix : public NVector<Vector3<C>,3> {
    public:
        Vector3<C> dot(Vector3<C> v) const;
        inline Vector3<C> operator *(Vector3<C> v) const{return dot(v);};
        Matrix<C> SymmetricInverse() const; /**<The inverse of a symmetric inverse*/
        C det() const; /**< The determinant */
};

template<class C>
C Matrix<C>::det() const{
    const Matrix<C> &M = *this;
    return 2*M[0][1]*M[0][2]*M[1][2] + M[0][0]*M[1][1]*M[2][2]
              - M[0][2]*M[0][2]*M[1][1] - M[0][0]*M[1][2]*M[1][2]
               - M[0][1]*M[0][1]*M[2][2];
};

template<class C>
Matrix<C> Matrix<C>::SymmetricInverse() const{
    const Matrix<C> &M = *this;
    C d = det();
    Matrix<C> I;
    I[0][0] = (M[1][1]*M[2][2]-M[1][2]*M[1][2]) / d;
    I[0][1] = (M[0][2]*M[1][2]-M[0][1]*M[2][2]) / d;
    I[1][0] = I[0][1];
    I[0][2] = (M[0][1]*M[1][2]-M[0][2]*M[1][1]) / d;
    I[2][0] = I[0][2];
    I[1][1] = (M[0][0]*M[2][2]-M[0][2]*M[0][2]) / d;
    I[1][2] = (M[0][1]*M[0][2]-M[0][0]*M[1][2]) / d;
    I[2][1] = I[1][2];
    I[2][2] = (M[0][0]*M[1][1]-M[0][1]*M[0][1]) / d;
    return I;
};

template<class C>
Vector3<C> Matrix<C>::dot(Vector3<C> vec) const{
    return Vector3<C>(this->get(0).dot(vec),this->get(1).dot(vec),this->get(2).dot(vec));
}

template <class T, unsigned int N>
Array<T, N>::Array(const Array<T, N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = rhs.get(i);
    }
}

template <class T, unsigned int N>
Array<T, N>::Array() {}

template <class T, unsigned int N>
Array<T, N>::Array(const T locs[N]) {
    for(unsigned int i=0; i < N; i++) vals[i] = locs[i];
}

template <class T, unsigned int N> template<class U>
Array<T, N>::Array(const Array<U,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = T(rhs.get(i));
    }
}

template <class T, unsigned int N>
NVector<T, N>::NVector(const NVector<T, N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = rhs.get(i);
    }
}

template <class T, unsigned int N>
NVector<T, N>::NVector() {}

template <class T, unsigned int N>
NVector<T, N>::NVector(const T locs[N]) {
    for(unsigned int i=0; i < N; i++) vals[i] = locs[i];
}

template <class T, unsigned int N> template<class U>
NVector<T, N>::NVector(const NVector<U,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] = T(rhs.get(i));
    }
}

template <class T, unsigned int N>
NVector<T, N> NVector<T, N>::operator-() const {
    NVector<T, N> newvec = NVector(*this);
    newvec *= -1;
    return newvec;
}

template <class T, unsigned int N>
NVector<T, N>& NVector<T, N>::operator-=(const NVector<T,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] -= rhs.get(i);
    }
    return *this;
}

template <class T, unsigned int N>
NVector<T, N>& NVector<T, N>::operator+=(const NVector<T,N> &rhs) {
    for(unsigned int i=0; i < N; i++){
        vals[i] += rhs.get(i);
    }
    return *this;
}

template <class T, unsigned int N> template <class U>
NVector<T,N>& NVector<T,N>::operator*=(const U rhs) {
    for(unsigned int i=0; i<N; i++){
        vals[i] = vals[i] * rhs;
        //set(i, T(this->get(i) * rhs));
    }
    return *this;
}

template <class T, unsigned int N> template <class U>
NVector<T,N>& NVector<T,N>::operator/=(const U rhs) {
    for(unsigned int i=0; i<N; i++){
        vals[i] = T(vals[i] / rhs);
    }
    return *this;
}

template <class U, unsigned int M>
ostream& operator<< (ostream& out, const NVector<U,M> &v){
    out << "[" << v.get(0);
    for(unsigned int i = 1; i < M; i++)
        out << ',' << v.get(i);
    return out << ']';
}

template <class T, unsigned int N>
T NumVector<T,N>::dot(const NumVector<T, N> &other) const{
    T m = 0;
    for(unsigned int i=0; i<N; i++){
        m += NVector<T,N>::get(i) * other.get(i);
    }
    return m;
}

template <class T, unsigned int N>
NumVector<T,N> NumVector<T,N>::perpto(const NumVector<T,N>& other) const {
    //~ NumVector<T,N> parallel = other * (this->dot(other)) / other.dot(other);
    return (*this) - other * ((dot(other)) / other.sq());
    //~ return other - parallel;
}

template <class T, unsigned int N>
void NumVector<T,N>::normalize(){
    T m = mag();
    *this /= m;
}

template <class T, unsigned int N>
NumVector<T,N> NumVector<T,N>::norm() const{
    T m = mag();
    return *this / m;
}


template <class T, unsigned int N>
inline T NumVector<T,N>::distance(const NumVector<T,N> &rhs) const{
    T sum = 0;
    for(unsigned int i=0; i<N; i++){
        sum += pow(NVector<T,N>::get(i) - rhs.get(i), 2);
    }
    return sqrt(sum);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline double NumVector<double,2>::distance(const NumVector<double,2> &rhs) const{
    double x=get(0) - rhs.get(0), y=get(1) - rhs.get(1);
    return hypot(x,y);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline long double NumVector<long double,2>::distance(const NumVector<long double,2> &rhs) const{
    long double x=get(0) - rhs.get(0), y=get(1) - rhs.get(1);
    return hypotl(x,y);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline double NumVector<double,3>::distance(const NumVector<double,3> &rhs) const{
    double x=get(0)-rhs.get(0), y=get(1)-rhs.get(1), z=get(2)-rhs.get(2);
    return sqrt(x*x + y*y + z*z);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline long double NumVector<long double,3>::distance(const NumVector<long double,3> &rhs) const{
    long double x=get(0)-rhs.get(0), y=get(1)-rhs.get(1), z=get(2)-rhs.get(2);
    return sqrtl(x*x + y*y + z*z);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline double NumVector<double,2>::mag() const{
    double x=get(0), y=get(1);
    return hypot(x,y);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline long double NumVector<long double,2>::mag() const{
    long double x=get(0), y=get(1);
    return hypotl(x,y);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline double NumVector<double,3>::mag() const{
    double x=get(0), y=get(1), z=get(2);
    return sqrt(x*x + y*y + z*z);
    // return NumVector<T,N>(*this - rhs).mag();
}

template<>
inline long double NumVector<long double,3>::mag() const{
    long double x=get(0), y=get(1), z=get(2);
    return sqrtl(x*x + y*y + z*z);
    // return NumVector<T,N>(*this - rhs).mag();
}

template <class U, unsigned int M>
ostream& operator<< (ostream& out, const NumVector<U,M> &v){
    out << "[" << v.get(0);
    for(int i = 1; i < M; i++)
        out << ',' << v.get(i);
    return out << ']';
}

template <class T>
Vector3<T> Vector3<T>::cross(const Vector3<T> &rhs) const {
    T newx = gety()*rhs.getz() - rhs.gety()*getz();
    T newy = getz()*rhs.getx() - rhs.getz()*getx();
    T newz = getx()*rhs.gety() - rhs.getx()*gety();
    return Vector3<T>(newx, newy, newz);
}

template <class T> template <class U>
Vector3<T>& Vector3<T>::operator*=(const U rhs){
    NVector<T,3>::vals[0] *= rhs;
    NVector<T,3>::vals[1] *= rhs;
    NVector<T,3>::vals[2] *= rhs;
    return *this;
}

template <class T> template <class U>
Vector3<T>& Vector3<T>::operator/=(const U rhs){
    NVector<T,3>::vals[0] /= rhs;
    NVector<T,3>::vals[1] /= rhs;
    NVector<T,3>::vals[2] /= rhs;
    return *this;
}

template <class U>
ostream& operator<< (ostream& out, const Vector3<U> v){
    out << "{" << v.get(0);
    for(int i = 1; i < 3; i++)
        out << ',' << v.get(i);
    return out << '}';
}

template <class T> template <class U>
Vector2<T>& Vector2<T>::operator*=(const U rhs){
    NVector<T,2>::vals[0] *= rhs;
    NVector<T,2>::vals[1] *= rhs;
    return *this;
}

template <class T> template <class U>
Vector2<T>& Vector2<T>::operator/=(const U rhs){
    NVector<T,2>::vals[0] /= rhs;
    NVector<T,2>::vals[1] /= rhs;
    return *this;
}

template <class T>
Vector2<T> Vector2<T>::rotate(uint i){
    if(i % 4 == 0) return (*this);
    else if(i % 4 == 3) return Vector2(gety(), -getx());
    else if(i % 4 == 2) return Vector2(-getx(), -gety());
    else return Vector2(-gety(), getx());
}

template <class U>
ostream& operator<< (ostream& out, const Vector2<U> v){
    out << "{" << v.get(0);
    for(int i = 1; i < 2; i++)
        out << ',' << v.get(i);
    return out << '}';
}

#endif
