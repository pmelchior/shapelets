#ifndef SHAPELENS_GRID_H
#define SHAPELENS_GRID_H

#include "../Typedef.h"
#include "Point.h"
#include "Shapes.h"
#include "CoordinateTransformation.h"
#include <boost/shared_ptr.hpp>

namespace shapelens {

/// Grid class.
/// The class defines the Grid on which an Image entity is defined.\n\n
/// The class distinguishes between 
/// - \b pixel coordinates (\p int numbers, defined at construction time) and
/// - \b World coordinates (\p data_t number, computed from the transformation set with apply()).

class Grid {
 public:
  /// Default constructor.
  Grid();
  /// Argumented constructor.
  /// Construct a Grid, starting at coordinates \f$(start_0,start_1)\f$ and traversing 
  /// \f$N_0\ (N_1)\f$ steps of given \p stepsize into positive direction.
  Grid(int start0, int start1, unsigned int N0, unsigned int N1);
  /// Set size of grid.
  void setSize(int start0, int start1, unsigned int N0, unsigned int N1);
  /// Return single World coordinate from pixel \p index.
  /// The World coordinate system ist set by calling apply().
  data_t operator() (long index, bool direction) const;
  /// Return World coordinates from pixel \p index.
  /// The World coordinate system ist set by calling apply().
  Point<data_t> operator() (long index) const;
  /// Set a coordinate transformation to grid points returned 
  /// by Grid::operator().
  void setWCS(const CoordinateTransformation& C);
  /// Reset to coordinate transformation.
  /// After a call to this function, World coordinates are image coordinates.
  void resetWCS();
  /// Get WCS transformation.
  const CoordinateTransformation& getWCS() const;
  /// Get effective pixel scale in WCS units.
  /// This is averaged over the entire footprint of the Grid.
  data_t getScaleFactor() const;
  /// Return starting position in given direction.
  int getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  int getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned long size() const;
  /// Get Polygon which contains the grid points in World Coordinates.
  Polygon<data_t> getSupport() const;
  /// Return rectangular bounding box in image coordinates.
  Rectangle<int> getBoundingBox() const;
  /// Find the image coordinates from the pixel index.
  /// The coordinate system is defined such, that the left lower corner of 
  /// the image has the coordinates \f$(start_0,start_1)\f$ (defined at 
  /// construction time).
  Point<int> getCoords(long pixel) const;
  /// Find the image coordinates from the World coordinates.
  /// The pixel coordinates are converted to \p int by rounding off,
  /// so the denote the left-lower corner of the pixel associated with \p P.\n\n
  /// The World coordinate system ist set by calling apply().
  Point<int> getCoords(const Point<data_t>& P) const;
  /// Get the pixel index from the \p P with image coordinates.
  /// Returns <tt>-1</tt> if the pixel is outside the image area.\n
  long getPixel(const Point<int>& P) const;
  /// Get the pixel index of the neighbor pixel using image coordinates \p P.
  /// Returns the pixel index or -1 if the pixel is outside the image area.\n
  /// For performance issues, the directions are defined as:
  /// \f[\begin{bmatrix}6 & 7 & 8 \\ 1 & 0 & 2 \\ 3 & 4 & 5\end{bmatrix}\f]
  long getNeighborPixel(const Point<int>& P, unsigned char direction) const;
  /// Get the pixel number of the neighbor pixel using its pixel number.
  /// Returns <tt>-1</tt> if the pixel is outside the image area.\n
  /// For performance issues, the directions are defined as:
  /// \f[\begin{bmatrix}6 & 7 & 8 \\ 1 & 0 & 2 \\ 3 & 4 & 5\end{bmatrix}\f]
  long getNeighborPixel(long pixel, unsigned char direction) const;

 private:
  int N0, N1, start0,start1;
  boost::shared_ptr<CoordinateTransformation> ct;
  NullTransformation nt;
  int round(data_t x) const;
};
} // end namespace

#endif
