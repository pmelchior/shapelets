#ifndef SHAPELENS_GRID_H
#define SHAPELENS_GRID_H

#include <math.h>
#include "../Typedef.h"
#include "Point.h"

namespace shapelens {

/// Grid class.
/// The class defines the Grid on which an Image entity is defined.\n\n
/// Grid positions are integer pixel numbers (starting from 0,0
/// by default).
///
/// \todo operator() could include WCS transformations

class Grid {
 public:
  /// Default constructor.
  Grid();
  /// Argumented constructor.
  /// Construct a Grid, starting at coordinates \f$(start_0,start_1)\f$ and traversing 
  /// \f$N_0\ (N_1)\f$ steps of given \p stepsize into positive direction.
  Grid(int start0, int start1, unsigned int N0, unsigned int N1);
  /// Index operator for const Grid.
  data_t operator() (unsigned long index, bool direction) const;
  /// Return the ith point as a Point.
  Point<data_t> operator() (unsigned long index) const;
  /// Return starting position in given direction.
  int getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  int getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned long size() const;
  /// Find the image coordinates from the pixel index.
  /// The coordinate system is defined such, that the left lower corner of 
  /// the image has the coordinates \f$(start_0,start_1)\f$ (defined at 
  /// construction time).
  Point<int> getCoords(unsigned long pixel) const;
  /// Get the pixel index from the \p P with image coordinates.
  unsigned long getPixel(const Point<int>& P) const;
  /// Get the pixel index of the neighbor pixel using image coordinates \p P.
  /// Returns the pixel index or -1 if the pixel is outside the image area.\n
  /// The directions (0..8) go clockwise from top(1) to top-left(8); 
  /// direction 0 is the pixel itself.
  long getNeighborPixel(const Point<int>& P, unsigned char direction) const;
  /// Get the pixel number of the neighbor pixel using its pixel number.
  long getNeighborPixel(unsigned long pixel, unsigned char direction) const;

 private:
  unsigned int N0, N1;
  int start0,start1;
};
} // end namespace

#endif
