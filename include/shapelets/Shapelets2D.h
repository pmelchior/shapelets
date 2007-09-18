#ifndef SHAPELETS2D_H
#define SHAPELETS2D_H

/// 2D Shapelet class.
/// Provides calculation of values of 2D Shapelets basis functions 
/// \f$B_{order0,order1}(x0,x1;\beta)\f$.

#include <Typedef.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <shapelets/Shapelets1D.h>

class Shapelets2D {
 public:
  /// Default constructor.
  Shapelets2D ();
  /// Constructor with scale size \f$\beta\f$ and orders.
  Shapelets2D (int order0, int order1, data_t beta);
  
  /// Return highest order of \f$B\f$ in direction (0/1).
  int getOrder (bool direction);
  /// Set the highest shapelets orders.
  void setOrders (int order0, int order1);
  /// Return \f$\beta\f$.
  data_t getBeta();
  /// Set \f$\beta\f$ to arbitrary value.
  void setBeta(data_t beta);
  /// Get smallest reproducible object size \f$\theta_{min}\f$.
  data_t getThetaMin(int order0, int order1);
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  data_t getThetaMax(int order0, int order1);
  /// Get integral over basis function  \f$B_{order0,order1}\f$.
  data_t integrate(int order0, int order1);
  /// Get integral over basis function  \f$B_{order}\f$ within range (x0min/x1min)..(x0max/x1max).
  /// see Paper III. eq. (82)
  data_t integrate(int order0, int order1, data_t x0min, data_t x0max, data_t x1min,data_t x1max);
  /// Evaluate \f$B_{order}((x0,x1);\beta)\f$.
  data_t eval(int order0, int order1, Point2D& x);

private:
  Shapelets1D S1D;
  data_t beta;
  int order0, order1;
};

#endif
