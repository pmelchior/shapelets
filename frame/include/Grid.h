#ifndef GRID_H
#define GRID_H

/// Grid class.

#include <math.h>
#include <Point2D.h>

class Grid {
 public:
  /// Default constructor.
  Grid();
  /// Argumented constructor for 1D grid.
  Grid(double start, double stop, double stepsize0);
  /// Argumented constructor 2D Grid.
  /// It's assumed that start0 < stop0 and start1 < stop1 holds.
  Grid(double start0, double stop0, double stepsize0, double start1, double stop1, double stepsize1);
  /// Index operator for const Grid.
  double operator() (unsigned int index, bool direction) const;
  /// Return the ith point as a Point2D.
  Point2D operator() (unsigned int index) const;
  /// Return stepsize between grid points in given direction.
  double getStepsize(bool direction) const;
  /// Return starting position in given direction.
  double getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  double getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned int size() const;
 private:
  unsigned int computeSize(bool direction) const;
  unsigned int axsize0, axsize1;
  Point2D start,stop,stepsize;
};

#endif
