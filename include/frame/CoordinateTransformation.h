#ifndef SHAPELENS_COORDINATETRANSFORMATION_H
#define SHAPELENS_COORDINATETRANSFORMATION_H

#include "Point.h"
#include <numla/NumMatrix.h>

namespace shapelens {
  /// Base class for coordinate transformations in 2D.
  /// \b CAUTION: Due to calculation of the inverse transformation,
  /// the correct behaviour can only be guaranteed for floating-point
  /// data types, both real and complex.
  template <class T>
    class CoordinateTransformation {
  public:
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const = 0;
    /// Multiply with another transformation from the right.
    /// This allows to form arbitrary products of transformations.
    CoordinateTransformation& operator*=(const CoordinateTransformation& C) {
      if (right == NULL)
	right = const_cast<CoordinateTransformation*>(&C);
      else
	right->right = const_cast<CoordinateTransformation*>(&C);
    }
  protected:
    CoordinateTransformation* right;
  };

  /// Class for rescaling transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P=s\cdot P\f$ with some scalar value \f$s\f$.
  template <class T>
    class ScalarTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    ScalarTransformation(T scale): s(scale) {
      CoordinateTransformation<T>::right = NULL;
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      if (CoordinateTransformation<T>::right != NULL)
	CoordinateTransformation<T>::right->transform(P);
      P *= s;
    }
    /// Get inverse transformation.
    ScalarTransformation getInverse() const {
      return ScalarTransformation(1./s);
    }
  private:
    T s;
  };

  /// Class for affine transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f[P=R_{out} + M\cdot\ (P-R_{in}),\f]
  /// where \f$M\f$ is a 2x2 matrix 
  template <class T>
    class AffineTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
  AffineTransformation(const NumMatrix<T>& M, const Point<T>& Rin, const Point<T>& Rout = Point<T>(0,0)) : M(M), Rin(Rin), Rout(Rout) {
      CoordinateTransformation<T>::right = NULL; 
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      if (CoordinateTransformation<T>::right != NULL)
	CoordinateTransformation<T>::right->transform(P);
      // store temporarily
      T p0 = P(0);
      P(0) = Rout(0) + M(0,0)*(P(0)-Rin(0)) + M(0,1)*(P(1)-Rin(1));
      P(1) = Rout(1) + M(1,0)*( p0 -Rin(0)) + M(1,1)*(P(1)-Rin(1));
    }
    /// Get inverse transformation.
    AffineTransformation getInverse() const {
      return AffineTransformation(M.invert(),Rout,Rin);
    }
  private:
    NumMatrix<T> M;
    Point<T> Rin, Rout;
  };
} // end namespace
#endif
