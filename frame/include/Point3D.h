#ifndef POINT3D_H
#define POINT3D_H

/// 3D Point class.
/// publicly inherited from uBLAS vectors

#include <boost/numeric/bindings/traits/ublas_vector.hpp>

template <class T>
class Point3D : public boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,3> >
{
  typedef boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,3> > Base_vector;
 public:
  Point3D () : Base_vector(3) {}
  template <class R> Point3D (R x0, R x1, R x2) : Base_vector(3) {
    Base_vector::operator()(0) = x0;
    Base_vector::operator()(1) = x1;
    Base_vector::operator()(2) = x2;
  }
  template <class R> Point3D (boost::numeric::ublas::vector_expression<R>& r) : Base_vector(r) { }
  template <class R> void operator=(boost::numeric::ublas::vector_expression<R>& r) {
    Base_vector::operator=(r);
  }
  template <class R> void operator=( Base_vector& r) {
    Base_vector::operator=(r);
  }
};

#endif
