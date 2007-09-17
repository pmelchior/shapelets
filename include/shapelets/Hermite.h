#ifndef HERMITE_H
#define HERMITE_H

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <boost/numeric/ublas/triangular.hpp>

/// Hermite polynomial class.
/// Provides calculation of values of Hermite polynomials and its coefficients
/// \f$H_n(x) = \sum_{2i+j=n} (-1)^i(2x)^j\frac{n!}{i!j!}\f$. \n\n
/// The code uses following recurrance relation for Hermite polynomials:
/// \f$H_{n+1}(x) = 2xH_n(x) - 2nH_{n-1}(x)\f$, where \f$H_0 = 1\f$.

class Hermite {
 public:
  /// Default constructor. 
  /// Does no computation of the coefficients.
  Hermite ();
  /// Constructor with arbitrary order \f$n\f$
  Hermite (unsigned int n);
  
  /// Return highest order of H.
  int getOrder ();
  /// Set the highest polynomial order \f$n\f$.
  void setOrder (unsigned int n);
  /// Return the coefficient of the \f$i\f$th power of the \f$n\f$th Hermite polynomial.
  /// Double values neccessary because these coefficients become very large at high orders.
  double getCoefficient(unsigned int n, unsigned int i);
  /// Evaluate \f$H_{n}(x)\f$.
  double eval (unsigned int n, double x);
  /// Print coeffs of \f$H_{n}\f$ to stdout.
  void printCoeffs (unsigned int n);

private:
  boost::numeric::ublas::triangular_matrix<double,boost::numeric::ublas::lower> HermiteCoeffs;
  int computed;
  void computeHermiteCoeffs(unsigned int n);
};

#endif

