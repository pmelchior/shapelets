#ifndef SHAPELENS_HERMITE_H
#define SHAPELENS_HERMITE_H

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <boost/numeric/ublas/triangular.hpp>
#include "../Typedef.h"

namespace shapelens {

/// Hermite polynomial class.
/// Provides calculation of values of Hermite polynomials and its coefficients
/// \f$H_n(x) = \sum_{2i+j=n} (-1)^i(2x)^j\frac{n!}{i!j!}\f$. \n\n
/// The code makes use of the following recurrance relation for Hermite polynomials:
/// \f$H_{n+1}(x) = 2xH_n(x) - 2nH_{n-1}(x)\f$ with \f$H_0 = 1\f$.

class Hermite {
 public:
  /// Default constructor. 
  /// Does no computation of the coefficients.
  Hermite ();
  /// Constructor with arbitrary order \f$n\f$
  Hermite (unsigned int n);
  
  /// Return highest order of H.
  int getOrder () const;
  /// Set the highest polynomial order \f$n\f$.
  void setOrder (unsigned int n);
  /// Return the coefficient of the \f$i\f$th power of the \f$n\f$th Hermite polynomial.
  /// Double values neccessary because these coefficients become very large at high orders.
  data_t getCoefficient(unsigned int n, unsigned int i) const;
  /// Evaluate \f$H_{n}(x)\f$.
  data_t eval (unsigned int n, data_t x);

private:
  boost::numeric::ublas::triangular_matrix<data_t,boost::numeric::ublas::lower> HermiteCoeffs;
  int computed;
  void computeHermiteCoeffs(unsigned int n);
};
} // end namespace
#endif


