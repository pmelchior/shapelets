#ifndef POINT2D_H
#define POINT2D_H

#include <boost/numeric/bindings/traits/ublas_vector.hpp>

namespace shapelens {

/// Templated 2D Point class.
/// publicly inherited from uBLAS vectors

template <class T>
class Point2D : public boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,2> >
{
  typedef boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,2> > Base_vector;
 public:
  Point2D () : Base_vector(2) {}
  template <class R> Point2D (R x0, R x1) : Base_vector(2) {
    Base_vector::operator()(0) = x0;
    Base_vector::operator()(1) = x1;
  }
  template <class R> Point2D (boost::numeric::ublas::vector_expression<R>& r) : Base_vector(r) { }
  template <class R> void operator=(boost::numeric::ublas::vector_expression<R>& r) {
    Base_vector::operator=(r);
  }
  template <class R> void operator=( Base_vector& r) {
    Base_vector::operator=(r);
  }
  // comparison operator for sorted containers of the STL
  // 2nd dimension first ensures efficient lookup for images
  template <class R>
    bool operator<(const Point2D<R>& b) const {
    if (Base_vector::operator()(1) < b(1) || (Base_vector::operator()(1) == b(1) && Base_vector::operator()(0) < b(0)))
      return true;
    else
      return false;
  }
};
} // end namespace
#endif
