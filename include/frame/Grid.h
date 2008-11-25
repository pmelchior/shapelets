#ifndef GRID_H
#define GRID_H

/// Templated grid class.
/// The class defines the Grid on which an Image entity is defined.\n\n
/// Traditionally, grid positions are integer pixel numbers (starting from 0,0
/// by default). Therefore \p Typedef.h defines \p GridT<int> as \p Grid.\n
/// If its required to work on floating point grids,
/// a single change in the file \p Typedef.h is thus sufficient.

#include <math.h>
#include <frame/Point2D.h>

template <class T>
class GridT {
 public:
  /// Default constructor.
  GridT();
  /// Argumented constructor.
  /// Construct a GridT, starting at coordinates \f$(start_0,start_1)\f$ and traversing 
  /// \f$N_0\ (N_1)\f$ steps of given \p stepsize into positive direction.
  GridT(T start0, T start1, T N0, T N1, T stepsize0 = 1, T stepsize1 = 1);
  /// Index operator for const GridT.
  T operator() (unsigned int index, bool direction) const;
  /// Return the ith point as a Point2D.
  Point2D<T> operator() (unsigned int index) const;
  /// Return stepsize between grid points in given direction.
  T getStepsize(bool direction) const;
  /// Return starting position in given direction.
  T getStartPosition(bool direction) const;
  /// Return stopping position in given direction.
  T getStopPosition(bool direction) const;
  /// Return numer of grid points in given direction.
  unsigned int getSize(bool direction) const;
  /// Return the number of grid points.
  unsigned int size() const;
  /// Set the x and y coordinate from the pixel number.
  /// The coordinate system is defined such, that the left lower corder of the image
  /// has the coordinates (0,0) and each pixel has unit length.
  void getCoords(unsigned int pixel, int& x, int& y) const;
  /// Get the pixel number from the x and y corrdinates.
  unsigned int getPixel(int x, int y) const;
  /// Get the pixel number of the neighbor pixel using image coordinates.
  /// Returns the pixel index or -1 if the pixel is outside the image area.\n
  /// The directions (0..8) go clockwise from top(1) to top-left(8); 
  /// direction 0 is the pixel itself.
  int getNeighborPixel(unsigned int pixel, int x, int y, unsigned int direction) const;
  /// Get the pixel number of the neighbor pixel using its pixel number.
  int getNeighborPixel(unsigned int pixel, unsigned int direction) const;

 private:
  int N0, N1;
  T start0,start1,stepsize0,stepsize1;
};

template <class T> 
inline GridT<T>::GridT() :
  N0(0),
  N1(0)
{
}

template <class T> 
inline GridT<T>::GridT(T start0, T start1, T N0, T N1, T stepsize0, T stepsize1) :
  start0(start0),
  start1(start1),
  stepsize0(stepsize0),
  stepsize1(stepsize1),
  N0(N0),
  N1(N1)
{
}

template <class T> 
inline T GridT<T>::operator() (unsigned int index, bool direction) const {
  int offset;
  if (direction) {
    offset = index/N0;
    return start1 + offset*stepsize1;
  }
  else {
    offset = index%N0;
    return start0 + offset*stepsize0;
  }
}

template <class T> 
inline Point2D<T> GridT<T>::operator() (unsigned int i) const {
  return Point2D<T>(operator()(i,0),operator()(i,1));
}

template <class T> 
inline T GridT<T>::getStepsize(bool direction) const {
  if (direction)
    return stepsize1;
  else
    return stepsize0;
}

template <class T> 
inline T GridT<T>::getStartPosition(bool direction) const {
  if (direction)
    return start1;
  else
    return start0;
}

template <class T> 
inline T GridT<T>::getStopPosition(bool direction) const {
  if (direction)
    return start1 + N1*stepsize1;
  else
    return start0 + N0*stepsize0;
}

template <class T> 
inline unsigned int GridT<T>::getSize(bool direction) const {
  if (direction)
    return N1;
  else
    return N0;
}

template <class T> 
inline unsigned int GridT<T>::size() const {
  return N0*N1;
}

template <class T> 
inline void GridT<T>::getCoords(unsigned int pixel, int& x, int& y) const {
  x = start0 + pixel%N0;
  y = start1 + pixel/N0;
}

template <class T> 
inline unsigned int GridT<T>::getPixel(int x, int y) const {
  return (unsigned int) (x-start0) + (y-start1)*N0;
}

template <class T> 
inline int GridT<T>::getNeighborPixel(unsigned int pixel, int x, int y, unsigned int direction) const {
  int index;
  switch(direction) {
  case 0: 
    // the pixel itself
    index = pixel;
    break;
  case 1: 
    if (y<N1) index = (y+1)*N0 + x ;  // top
    else index = -1;
    break;
  case 2:
    if (y<N1 && x<N0) index = (y+1)*N0 + x + 1;  // top right
    else index = -1;
    break;
  case 3:
    if (x<N0) index = y*N0 + x + 1;  // right neighbour
    else index = -1;
    break;
  case 4: 
    if (y>0 && x<N0) index = (y-1)*N0 + x + 1;  // bottom right
    else index = -1;
    break;  
  case 5: 
    if (y>0) index = (y-1)*N0 + x;  // bottom
    else index = -1;
    break;
  case 6: 
    if (y>0 && x>0) index = (y-1)*N0 + x - 1;  // bottom left
    else index = -1;
    break;   
  case 7: 
    if (x>0) index = y*N0 + x - 1; // left
    else index = -1;
    break;
  case 8: 
    if (y<N1 && x>0) index = (y+1)*N0 + x - 1;  // top left
    else index = -1;
    break;  
  }
  return index;
}

template <class T> 
inline int GridT<T>::getNeighborPixel(unsigned int pixel, unsigned int direction) const {
  int x,y;
  getCoords(pixel,x,y);
  return getNeighborPixel(pixel,x,y,direction);
}


#endif
