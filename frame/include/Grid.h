#ifndef GRID_H
#define GRID_H

/// Grid class.
/// Provides a grid of points as a matrix.
/// By now grids in 1D and 2D are implemented.

#include <math.h>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <Point2D.h>

class Grid {
 public:
  /// Default constructor.
  /// Don't use it!
  Grid();
  /// Build 1D Grid.
  /// It's assumed that start < stop holds.
  Grid(double start, double stop, double stepsize);
  /// Build 2D Grid.
  /// It's assumed that start0 < stop0 and start1 < stop1 holds.
  Grid(double start0, double stop0, double stepsize0, double start1, double stop1, double stepsize1);
  /// Copy operator.
  Grid& operator= (const Grid& h);
  /// Assignment operator.
  double& operator() (const unsigned int i, const unsigned int j);
  /// Index operator for const Grid.
  const double& operator() (const unsigned int i, const unsigned int j) const;
  /// Overloaded assignment operator.
  /// Return the ith point as a Point2D.
  const Point2D& operator() (const unsigned int i);
  /// Return stepsize between grid points in given direction.
  double getStepsize(bool direction) const;
  /// Return starting position in given direction.
  double getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  double getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  int getSize(bool direction) const;
  /// Return the number of grid points.
  int size() const;
 private:
  boost::numeric::ublas::matrix<double> grid;
  boost::numeric::ublas::vector<double> start,stop,stepsize;
  Point2D point;
};

#endif
