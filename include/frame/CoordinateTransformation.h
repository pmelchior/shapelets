#ifndef SHAPELENS_COORDINATETRANSFORMATION_H
#define SHAPELENS_COORDINATETRANSFORMATION_H

#include "Point.h"
#include "../Typedef.h"
#include <numla/NumMatrix.h>
#include <boost/shared_ptr.hpp>

namespace shapelens {
  /// Base class for coordinate transformations in 2D.
  class CoordinateTransformation {
  public:
    /// Destructor.
    ~CoordinateTransformation();
    /// Chain a succeding transformation.
    /// The new transformation will be made of the product
    /// <tt>(*this) -> (*this) * C</tt>
    void operator*=(const CoordinateTransformation& C);
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const = 0;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const = 0;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const = 0;
    /// Stack of transformation to be applied after another
    std::list<boost::shared_ptr<CoordinateTransformation> > stack;
  protected:
    /// Apply all transformations from the stack.
    void stack_transform(Point<data_t>& P) const;
    /// Apply all inverse transformation from the stack in reverse order.
    void stack_inverse_transform(Point<data_t>& P) const;
    /// Set the stack of \p ct to \p this.stack .
    //CoordinateTransformation* clone_stack(CoordinateTransformation* ct) const;
  };

  /// Empty/Null transformations in 2D.
  /// This transformations will leave all points unchanged.
  class NullTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    NullTransformation();
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const;
  };

  /// Class for rescaling transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=s\cdot P\f$ with some scalar value \f$s\f$.
  class ScalarTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    ScalarTransformation(data_t scale);
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const;
  private:
    data_t s;
  };

  /// Class for translation transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=P+ \Delta P\f$ with some shift direction \f$\Delta P\f$.
  class ShiftTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    ShiftTransformation(const Point<data_t>& dP_);
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const;
  private:
    Point<data_t> dP;
  };

  /// Class for linear transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f[P^\prime= M\cdot\ P,\f]
  /// where \f$M\f$ is a 2x2 matrix 
  class LinearTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    LinearTransformation(const NumMatrix<data_t>& M_);
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const;
  private:
    NumMatrix<data_t> M, M_1;
  };

  /// Class for lensing transformations in 2D.
  class LensingTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    LensingTransformation(data_t kappa_, complex<data_t> gamma);
    /// Constructor for 2nd order lensing transformation.
    LensingTransformation(data_t kappa_, complex<data_t> gamma, complex<data_t> F, complex<data_t> G);
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation* clone() const;
  private:
    bool flex;
    data_t kappa, gamma1, gamma2, D111, D112, D121, D122, D211, D212, D221, D222;
    
  };
} // end namespace
#endif
