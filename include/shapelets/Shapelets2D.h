#ifndef SHAPELETS2D_H
#define SHAPELETS2D_H

/// 2D Shapelet class.
/// Provides calculation of values of 2D Shapelets basis functions 
/// \f$B_{order0,order1}(x0,x1;\beta)\f$.

#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <shapelets/Shapelets1D.h>

class Shapelets2D {
 public:
  /// Default constructor.
  Shapelets2D ();
  /// Constructor with scale size \f$\beta\f$ and orders.
  Shapelets2D (int order0, int order1, double beta);
  
  /// Return highest order of \f$B\f$ in direction (0/1).
  int getOrder (bool direction);
  /// Set the highest shapelets orders.
  void setOrders (int order0, int order1);
  /// Return \f$\beta\f$.
  double getBeta();
  /// Set \f$\beta\f$ to arbitrary value.
  void setBeta(double beta);
  /// Get smallest reproducible object size \f$\theta_{min}\f$.
  double getThetaMin(int order0, int order1);
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  double getThetaMax(int order0, int order1);
  /// Get integral over basis function  \f$B_{order0,order1}\f$.
  double integrate(int order0, int order1);
  /// Get integral over basis function  \f$B_{order}\f$ within range (x0min/x1min)..(x0max/x1max).
  /// see Paper III. eq. (82)
  double integrate(int order0, int order1, double x0min, double x0max, double x1min,double x1max);
  /// Evaluate \f$B_{order}((x0,x1);\beta)\f$.
  double eval(int order0, int order1, Point2D& x);

private:
  Shapelets1D S1D;
  double beta;
  int order0, order1;
};

#endif
