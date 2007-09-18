#ifndef PROFILE_H
#define PROFILE_H

/// Class for Profile plots from ShapeletsImage.
/// \todo implement profile thru center pixel in 4 directions

#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <NumVector.h>
#include <Typedef.h>
#include <frame/Point2D.h>

class Profile {
 public:
  /// Default constructor.
  Profile();
  /// Argumented constructor.
  /// Defines the starting and the ending point.
  /// The center will be exactly in the middle between the starting and ending point.
  Profile(const Point2D& start, const Point2D& end);
  /// Return distance of the ith pixel of the Profile to the center.
  data_t getDistance(unsigned int i);
  /// Return the value of the ith pixel of the Profile.
  data_t getValue(unsigned int i);
  /// Return starting point.
  Point2D& getStart();
  /// Return ending point.
  Point2D& getEnd();
  /// Return center point.
  Point2D& getCenter();
  /// Calculate the entries of the profile plot.
  /// The distance will be set to the signed radius from the center to the ith pixel,
  /// the value is simply the value of the ith pixel.\n
  /// The axsize is needed for converting positions to indices, because data 
  /// is arranged as a vector.
  void calculate(NumVector<data_t>& data, int axsize);
  /// Return the number of pixels contained in the Profile.
  int size();
 private:
  Point2D start, stop, center;
  int pixels;
  boost::numeric::ublas::matrix<data_t> distance_value;
};

#endif
