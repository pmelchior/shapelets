#ifndef COMPLEXDOUBLE_H
#define COMPLEXDOUBLE_H

#include <complex>

/// Enhanced <tt>complex<double></tt> class.
/// This class provides an typecast operator from <tt>complex<double></tt> to <tt>double</tt> which is missing
/// in the base class.

class ComplexDouble : public complex<double> {
 public:
  /// Default constructor.
  ComplexDouble() : complex<double>() {
  }
  /// Argumented constructor to form a complex number from two double numbers.
  ComplexDouble(double real, double imag=0) : complex<double>(real,imag) {
  }
  /// Copy constructor.
  ComplexDouble (const complex<double>& c) {
    complex<double>::operator=(c);
  }
  /// Typecast operator.
  /// Returns the real part of the complex number.
  operator double() const {
    return std::real(*this);
  }
};

#endif
