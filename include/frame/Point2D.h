#ifndef POINT2D_H
#define POINT2D_H

/// 2D Point class.
/// publicly inherited from uBLAS vectors

#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <Typedef.h>

class Point2D : public boost::numeric::ublas::vector<data_t, boost::numeric::ublas::bounded_array<data_t,2> >
{
  typedef boost::numeric::ublas::vector<data_t, boost::numeric::ublas::bounded_array<data_t,2> > Base_vector;
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
};

#endif
