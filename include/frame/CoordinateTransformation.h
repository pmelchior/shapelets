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
    virtual void transform(Point<T>& P) const {} // do nothing
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const = 0;
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const  = 0;
  };

  /// Empty/Null transformations in 2D.
  /// This transformations will leave all points unchanged.
  template <class T>
    class NullTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    NullTransformation() {
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return new NullTransformation();
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return new NullTransformation();
    }
  };

  /// Class for rescaling transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=s\cdot P\f$ with some scalar value \f$s\f$.
  template <class T>
    class ScalarTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    ScalarTransformation(T scale): s(scale) {
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      P *= s;
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return new ScalarTransformation(s);
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return new ScalarTransformation(1./s);
    }
  private:
    T s;
  };

  /// Class for translation transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=P- \Delta P\f$ with some shift direction \f$\Delta P\f$.
  template <class T>
    class ShiftTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    ShiftTransformation(const Point<T>& delta): dP(delta) {
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      P -= dP;
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return new ShiftTransformation(dP);
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return new ShiftTransformation(-dP);
    }
  private:
    Point<T> dP;
  };

  /// Class for linear transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f[P^\prime= M\cdot\ P,\f]
  /// where \f$M\f$ is a 2x2 matrix 
  template <class T>
    class LinearTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    LinearTransformation(const NumMatrix<T>& M_): M(M) {
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      // store temporarily
      T p0 = P(0);
      P(0) = M(0,0)*P(0) + M(0,1)*P(1);
      P(1) = M(1,0)* p0  + M(1,1)*P(1);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return new LinearTransformation(M);
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return new LinearTransformation(M.invert());
    }
  private:
    NumMatrix<T> M;
  };

  /// Class for affine transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f[P^\prime=R_{out} + M\cdot\ (P-R_{in}),\f]
  /// where \f$M\f$ is a 2x2 matrix 
  template <class T>
    class AffineTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
  AffineTransformation(const NumMatrix<T>& M, const Point<T>& Rin, const Point<T>& Rout = Point<T>(0,0)) : M(M), Rin(Rin), Rout(Rout) {
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      // store temporarily
      T p0 = P(0);
      P(0) = Rout(0) + M(0,0)*(P(0)-Rin(0)) + M(0,1)*(P(1)-Rin(1));
      P(1) = Rout(1) + M(1,0)*( p0 -Rin(0)) + M(1,1)*(P(1)-Rin(1));
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return new AffineTransformation(*this);
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return new AffineTransformation(M.invert(),Rout,Rin);
    }
  private:
    NumMatrix<T> M;
    Point<T> Rin, Rout;
  };
} // end namespace
#endif
