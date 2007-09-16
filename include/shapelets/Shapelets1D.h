#ifndef SHAPELETS1D_H
#define SHAPELETS1D_H

#include <shapelets/Hermite.h>

/// 1D Shapelet class.
/// Provides calculation of values of 1D Shapelets basis functions 
/// \f$B_n(x;\beta) = \beta^{-\frac{1}{2}}\phi_n(\beta^{-1}x)\f$ where
/// \f$\phi_n(x) = [2^n\pi^{\frac{1}{2}}n!]^{-\frac{1}{2}} H_n(x) \exp(-\frac{x^2}{2})\f$ 
/// and \f$H_n(x)\f$ is the nth order Hermite polynomial. 

class Shapelets1D {
 public:
  /// Default constructor.
  Shapelets1D ();
  /// Constructor with scale size \f$\beta\f$ and maximal order
  Shapelets1D (int order, double beta);
  
  /// Return highest order of B.
  int getOrder ();
  /// Set the highest order of B.
  void setOrder (int order);
  /// Return \f$\beta\f$.
  double getBeta();
  /// Set \f$\beta\f$ to arbitrary value.
  void setBeta(double beta);
  /// Get smallest reproducible object size \f$\theta_{min}\f$.
  double getThetaMin(int order);
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  double getThetaMax(int order);
  /// Get integral over basis function \f$B_{order}\f$.
  double integrate(int order);
  /// Get integral over basis fuction  \f$B_{order}\f$ within interval xmin .. xmax.
  double integrate(int order, double xmin, double xmax);
  /// Evaluate \f$B_{order}(x;\beta)\f$.
  double eval (int order, double x);

private:
  Hermite H;
  double beta,sqrt_beta;
};

#endif
