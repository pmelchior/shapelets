#ifndef SHAPELENS_POINT2D_H
#define SHAPELENS_POINT2D_H

#include <numla/boost/numeric/bindings/traits/ublas_vector.hpp>

namespace shapelens {

/// Templated 2D Point class.
/// publicly inherited from uBLAS vectors

template <class T>
class Point : public boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,2> >
{
  typedef boost::numeric::ublas::vector<T, boost::numeric::ublas::bounded_array<T,2> > Base_vector;
 public:
  /// Constructor.
  Point () : Base_vector(2) {}
  /// Constructor with given coordinates.
  template <class R> Point (R x0, R x1) : Base_vector(2) {
    Base_vector::operator()(0) = x0;
    Base_vector::operator()(1) = x1;
  }
  /// Copy constructor from base class.
  template <class R> Point (const boost::numeric::ublas::vector_expression<R>& r) : Base_vector(r) { }
  /// Copy operator from base class.
  template <class R> void operator=(const boost::numeric::ublas::vector_expression<R>& r) {
    Base_vector::operator=(r);
  }
  /// Copy operator from base-class.
  template <class R> void operator=( Base_vector& r) {
    Base_vector::operator=(r);
  }
  /// Comparison operator.
  /// This is important for sorted containers of the STL.
  /// To ensure efficient lookups for image-type data, the points are
  /// ordered according to their 2nd dimension
  template <class R>
    bool operator<(const Point<R>& b) const {
    if (Base_vector::operator()(1) < b(1) || (Base_vector::operator()(1) == b(1) && Base_vector::operator()(0) < b(0)))
      return true;
    else
      return false;
  }
  /// Equality operator.
  template <class R>
    bool operator==(const Point<R>& b) const {
    if (Base_vector::operator()(0) == b(0) && Base_vector::operator()(1) == b(1))
      return true;
    else
      return false;
  }
  /// Return pointer to data storage.
  T* c_array() {
    return boost::numeric::bindings::traits::vector_storage(*this);
  }
  /// Return pointer to data storage.
  const T* c_array() const {
    return boost::numeric::bindings::traits::vector_storage(*this);
  }
};
} // end namespace
#endif
