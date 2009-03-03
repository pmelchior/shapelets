#ifndef MOMENTS_H
#define MOMENTS_H

#include <Typedef.h>

namespace shapelens {
//
/// Container class for Quadrupole moments of two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.

class Quadrupole {
 public:
  /// Constructor.
  Quadrupole();
  /// Access operator with \f$Q_{ij} = Q_{ji}\f$.
  data_t& operator()(bool i, bool j);
  /// Access operator with \f$Q_{ij} = Q_{ji}\f$.
  const data_t& operator()(bool i, bool j) const;
 private:
  data_t m[3];
};

/// Container class for Octupole moments of two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Octupole {
  public:
  /// Constructor.
  Octupole();
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k) const;
 private:
  data_t m[4];
};

/// Container class for Hexadecupole moments of two-dimensional quantity.
/// The class exploits the fact, that moments are invariant under 
/// index permutation.
class Hexadecupole {
  public:
  Hexadecupole();
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  data_t& operator()(bool i, bool j, bool k, bool l);
  /// Access operator.
  /// The access will result in the same result, when indices are only permuted.
  const data_t& operator()(bool i, bool j, bool k, bool l) const;
 private:
  data_t m[5];
};
} // end namespace
#endif
