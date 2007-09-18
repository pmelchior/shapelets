#ifndef GRID_H
#define GRID_H

/// Grid class.

#include <math.h>
#include <Typedef.h>
#include <frame/Point2D.h>

class Grid {
 public:
  /// Default constructor.
  Grid();
  /// Argumented constructor for 1D grid.
  Grid(data_t start, data_t stop, data_t stepsize0);
  /// Argumented constructor 2D Grid.
  /// It's assumed that start0 < stop0 and start1 < stop1 holds.
  Grid(data_t start0, data_t stop0, data_t stepsize0, data_t start1, data_t stop1, data_t stepsize1);
  /// Index operator for const Grid.
  data_t operator() (unsigned int index, bool direction) const;
  /// Return the ith point as a Point2D.
  Point2D operator() (unsigned int index) const;
  /// Return stepsize between grid points in given direction.
  data_t getStepsize(bool direction) const;
  /// Return starting position in given direction.
  data_t getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  data_t getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned int size() const;
  /// Set the x and y coordinate from the pixel number.
  /// The coordinate system is defined such, that the left lower corder of the image
  /// has the coordinates (0,0) and each pixel has unit length.
  void getCoords(unsigned int pixel, unsigned int& x, unsigned int& y) const;
  /// Get the pixel number from the x and y corrdinates.
  unsigned int getPixel(unsigned int x, unsigned int y) const;
  /// Get the pixel number of the neighbor pixel using image coordinates.
  /// Returns the pixel index or -1 if the pixel is outside the image area.\n
  /// The directions (0..8) go clockwise from top(1) to top-left(8); 
  /// direction 0 is the pixel itself.
  int getNeighborPixel(unsigned int pixel, unsigned int x, unsigned int y, unsigned int direction) const;
  /// Get the pixel number of the neighbor pixel using its pixel number.
  int getNeighborPixel(unsigned int pixel, unsigned int direction) const;

 private:
  unsigned int computeSize(bool direction) const;
  unsigned int axsize0, axsize1;
  Point2D start,stop,stepsize;
};

#endif
