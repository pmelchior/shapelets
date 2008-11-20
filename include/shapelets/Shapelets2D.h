#ifndef SHAPELETS2D_H
#define SHAPELETS2D_H

/// 2D Shapelet class.
/// Provides calculation of values of 2D Shapelets basis functions 
/// \f$B_{n_0,n_1}(x_0,x_1;\beta) \equiv B_{n_0}(x_0;\beta)B_{n_1}(x_1;\beta)\f$, where
/// \f$B_{n_i}(x_i;\beta)\f$ is the Shapelets1D basis function.

#include <Typedef.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <shapelets/Shapelets1D.h>

class Shapelets2D {
 public:
  /// Default constructor.
  Shapelets2D ();
  /// Constructor with scale size \f$\beta\f$.
  Shapelets2D (data_t beta);
  /// Return \f$\beta\f$.
  data_t getBeta() const;
  /// Set \f$\beta\f$ to arbitrary value.
  void setBeta(data_t beta);
  /// Get smallest reproducible object size \f$\theta_{min}\f$.
  data_t getThetaMin(int order0, int order1) const;
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  data_t getThetaMax(int order0, int order1) const;
  /// Get integral over basis function  \f$B_{n_0,n_1}\f$.
  data_t integrate(int order0, int order1) const;
  /// Get integral over basis function  \f$B_{n_0,n_1}\f$ within the area enclosed by
  /// \f$(x_0^{min}/x_1^{min})\ ..\ (x_0^{max}/x_1^{max})\f$.
  /// see Paper III. eq. (82)
  data_t integrate(int n0, int n1, data_t x0min, data_t x0max, data_t x1min,data_t x1max);
  /// Evaluate \f$B_{n_0,n_1}((x_0,x_1);\beta)\f$.
  data_t eval(int n0, int n1, Point2D<data_t>& x);

private:
  Shapelets1D S1D;
};

#endif
