#ifndef SHAPELENS_GRID_H
#define SHAPELENS_GRID_H

#include <math.h>
#include "../Typedef.h"
#include "Point.h"
#include "Shapes.h"
#include "CoordinateTransformation.h"

namespace shapelens {

  /// WCS class.
  class WCS {
  public:
    WCS();
    /// CoordindateTransformation from pixel to World coordinates
    boost::shared_ptr<CoordinateTransformation<data_t> > CT;
    /// Inverse transformation, from World to pixel coordinates
    boost::shared_ptr<CoordinateTransformation<data_t> > CT_;
  };


/// Grid class.
/// The class defines the Grid on which an Image entity is defined.\n\n
/// The class distinguishes between 
/// - \b pixel coordinates (\p int numbers, defined at construction time) and
/// - \b World coordinates (\p data_t number, computed from eh transformation set with apply()).

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
  data_t operator() (unsigned long index, bool direction) const;
  /// Return World coordinates from pixel \p index.
  /// The World coordinate system ist set by calling apply().
  Point<data_t> operator() (unsigned long index) const;
  /// Apply coordinate transformation to grid points returned 
  /// by Grid::operator().
  void apply(const CoordinateTransformation<data_t>& C);
  /// Return whether WCS information is set.
  bool hasWCS() const;
  /// Get WCS transformation.
  const WCS& getWCS() const;
  /// Return starting position in given direction.
  int getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  int getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned long size() const;
  /// Return rectangular bounding box in Woorld Coordinates.
  Rectangle<data_t> getBoundingBox() const;
  /// Find the image coordinates from the pixel index.
  /// The coordinate system is defined such, that the left lower corner of 
  /// the image has the coordinates \f$(start_0,start_1)\f$ (defined at 
  /// construction time).
  Point<int> getCoords(unsigned long pixel) const;
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
  /// The directions (0..8) go clockwise from top(1) to top-left(8); 
  /// direction 0 is the pixel itself.
  long getNeighborPixel(const Point<int>& P, unsigned char direction) const;
  /// Get the pixel number of the neighbor pixel using its pixel number.
  /// Returns <tt>-1</tt> if the pixel is outside the image area.
  long getNeighborPixel(unsigned long pixel, unsigned char direction) const;

 private:
  unsigned int N0, N1;
  int start0,start1;
  WCS wcs;
  int round(data_t x) const;
};
} // end namespace

#endif
