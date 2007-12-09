#ifndef SHAPELETS1D_H
#define SHAPELETS1D_H

#include <Typedef.h>
#include <shapelets/Hermite.h>

/// 1D Shapelet class.
/// Provides calculation of values of 1D Shapelets basis functions 
/// \f$B_n(x;\beta) = \beta^{-\frac{1}{2}}\phi_n(\beta^{-1}x)\f$ where
/// \f$\phi_n(x) = [2^n\pi^{\frac{1}{2}}n!]^{-\frac{1}{2}} H_n(x) \exp(-\frac{x^2}{2})\f$ 
/// and \f$H_n(x)\f$ is the \f$n\f$th order Hermite polynomial. 

class Shapelets1D {
 public:
  /// Default constructor.
  Shapelets1D ();
  /// Constructor with scale size \f$\beta\f$ and maximal order
  Shapelets1D (int order, data_t beta);
  
  /// Return highest order of B.
  int getOrder ();
  /// Set the highest order of B.
  void setOrder (int order);
  /// Return \f$\beta\f$.
  data_t getBeta();
  /// Set \f$\beta\f$ to arbitrary value.
  void setBeta(data_t beta);
  /// Get smallest reproducible object size \f$\theta_{min}\f$.
  data_t getThetaMin(int order);
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  data_t getThetaMax(int order);
  /// Get integral over basis function \f$B_{order}\f$.
  data_t integrate(int order);
  /// Get integral over basis fuction  \f$B_{order}\f$ within interval xmin .. xmax.
  data_t integrate(int order, data_t xmin, data_t xmax);
  /// Evaluate \f$B_{order}(x;\beta)\f$.
  data_t eval (int order, data_t x);

private:
  Hermite H;
  data_t beta,sqrt_beta;
};

#endif
