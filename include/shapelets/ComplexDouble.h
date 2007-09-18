#ifndef COMPLEXDOUBLE_H
#define COMPLEXDOUBLE_H

#include <complex>
#include <Typedef.h>

/// Enhanced <tt>complex<data_t></tt> class.
/// This class provides an typecast operator from <tt>complex<data_t></tt> to <tt>data_t</tt> which is missing
/// in the base class.

class ComplexDouble : public complex<data_t> {
 public:
  /// Default constructor.
  ComplexDouble() : complex<data_t>() {
  }
  /// Argumented constructor to form a complex number from two data_t numbers.
  ComplexDouble(data_t real, data_t imag=0) : complex<data_t>(real,imag) {
  }
  /// Copy constructor.
  ComplexDouble (const complex<data_t>& c) {
    complex<data_t>::operator=(c);
  }
  /// Typecast operator.
  /// Returns the real part of the complex number.
  operator data_t() const {
    return std::real(*this);
  }
};

#endif
