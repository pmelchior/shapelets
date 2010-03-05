#ifndef SHAPELENS_MOMENTS_H
#define SHAPELENS_MOMENTS_H

#include "../Typedef.h"
#include "Object.h"

namespace shapelens {

/// Container class for monopole \f$l=0\f$ moments of a two-dimensional quantity.
class Moment0 {
 public:
  /// Constructor.
  Moment0();
  /// Constructor from on Object.
  Moment0(const Object& obj);
  /// Access operator.
  data_t& operator()(bool i = 0);
  /// Access operator.
  const data_t& operator()(bool i=0) const;
  //private:
  data_t M;
};


/// Container class for dipole \f$l=1\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment1 {
 public:
  /// Constructor.
  Moment1();
  /// Constructor from an Object.
  Moment1(const Object& obj);
  /// Access operator.
  data_t& operator()(bool i);
  /// Access operator.
  const data_t& operator()(bool i) const;
  //private:
  data_t M[2];
};


/// Container class for quadrupole \f$l=2\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment2 {
 public:
  /// Constructor.
  Moment2();
  /// Constructor from an Object.
  Moment2(const Object& obj);
  /// Access operator with \f$Q_{ij} = Q_{ji}\f$.
  data_t& operator()(bool i, bool j);
  /// Access operator with \f$Q_{ij} = Q_{ji}\f$.
  const data_t& operator()(bool i, bool j) const;
  //private:
  data_t M[3];
};

/// Container class for octupole \f$l=3\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment3 {
  public:
  /// Constructor.
  Moment3();
  /// Constructor from an Object.
  Moment3(const Object& obj);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k) const;
  //private:
  data_t M[4];
};

/// Container class for hexadecupole \f$l=4\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment4 {
  public:
  /// Constructor.
  Moment4();
  /// Constructor from an Object.
  Moment4(const Object& obj);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k, bool l);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k, bool l) const;
  //private:
  data_t M[5];
};

/// Container class for \f$l=5\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment5 {
  public:
  /// Constructor.
  Moment5();
  /// Constructor from an Object.
  Moment5(const Object& obj);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k, bool l, bool m);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k, bool l, bool m) const;
  //private:
  data_t M[6];
};

/// Container class for \f$l=6\f$ moments of a two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Moment6 {
  public:
  /// Constructor.
  Moment6();
  /// Constructor from an Object.
  Moment6(const Object& obj);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k, bool l, bool m, bool n);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k, bool l, bool m, bool n) const;
  //private:
  data_t M[7];
};

/// Container class for arbitrary moments of a two-dimensional quantity.
/// Access is provided by an index or by the moment power in x and y.\n
/// All moments are computed in one sweep over the data.
class MomentsOrdered : public NumVector<data_t> {
public:
  /// Default constructor
  MomentsOrdered();
  /// Constructor for moments up to order \p N.
  MomentsOrdered(int N);
  /// Constructor for moments up to order \p N.
  /// The moments are populated from \p obj.
  MomentsOrdered(const Object& obj, int N);
  /// Access operator for \f$\langle x^{p_x}\, y^{p_y}\rangle\f$.
  data_t& operator()(unsigned int px, unsigned int py);
  /// Access operator for \f$\langle x^{p_x}\, y^{p_y}\rangle\f$.
  const data_t& operator()(unsigned int px, unsigned int py) const;
  /// Get maximum moment order.
  int getOrder() const;
  /// Get vector index of moment \f$\langle x^{p_x}\, y^{p_y}\rangle\f$ from 
  /// the powers.
  int getIndex(unsigned int px, unsigned int py) const;
 private:
  unsigned int pyramid_num(int n) const;
  int N;
};



} // end namespace
#endif
