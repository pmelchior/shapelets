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
    /// Destructor.
    ~CoordinateTransformation() {
      if (next != NULL)
        delete next;
    }
    /// Chain a succeding transformation.
    /// The new transformation will be made of the product
    /// <tt>(*this) -> (*this) * C</tt>
    void operator*=(const CoordinateTransformation<T>& C) {
      CoordinateTransformation<T>* tmp = this;
      while (tmp->next != NULL) {
	tmp = this->next;
	}
      tmp->next = C.clone();
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const = 0;
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const = 0;
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const  = 0;
  protected:
    /// Points to succeding transformation.
    CoordinateTransformation<T>* next;
    /// Prepare correct deep copy from clone pointer.
    CoordinateTransformation<T>* clone_pointer(CoordinateTransformation<T>* ct)  const {
      if (next != NULL)
        ct->next = next->clone();
      else
        ct->next = NULL;
      return ct;	
    }
    /// Prepare correct deep copy of inverse pointer.
    CoordinateTransformation<T>* inverse_pointer(CoordinateTransformation<T>* inv)  const {
      CoordinateTransformation<T>* tmp = this->clone();
      while(tmp->next != NULL) {
	CoordinateTransformation<T>* next_inv = tmp->next->getInverse();
	next_inv->next = inv;
	// traverse linked list
	tmp = tmp->next;
        inv = next_inv;
      }
      return inv;
    }
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
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return clone_pointer(new NullTransformation());
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return inverse_pointer(new NullTransformation());
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
      this->next = NULL;
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      P *= s;
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return clone_pointer(new ScalarTransformation(s));
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return inverse_pointer(new ScalarTransformation(1./s));
    }
  private:
    T s;
  };

  /// Class for translation transformations in 2D.
  /// A Point \f$P\f$ will be transformed according to
  /// \f$P^\prime=P+ \Delta P\f$ with some shift direction \f$\Delta P\f$.
  template <class T>
    class ShiftTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
    ShiftTransformation(const Point<T>& delta): dP(delta) {
      this->next = NULL;
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      P += dP;
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return clone_pointer(new ShiftTransformation(dP));
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return inverse_pointer(new ShiftTransformation(-dP));
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
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return clone_pointer(new LinearTransformation(M));
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return inverse_pointer(new LinearTransformation(M.invert()));
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
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      return clone_pointer(new AffineTransformation(*this));
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      return inverse_pointer(new AffineTransformation(M.invert(),Rout,Rin));
    }
  private:
    NumMatrix<T> M;
    Point<T> Rin, Rout;
  };

  /// Class for lensing transformations in 2D.
  template <class T>
    class LensingTransformation : public CoordinateTransformation<T> {
  public:
    /// Constructor.
  LensingTransformation(data_t kappa_, complex<data_t> gamma) : 
    kappa(kappa_), flex(false) {
      gamma1 = real(gamma);
      gamma2 = imag(gamma);
    }
    /// Constructor for 2nd order lensing transformation.
  LensingTransformation(data_t kappa_, complex<data_t> gamma, complex<data_t> F, complex<data_t> G) : 
    kappa(kappa_), flex(true) {
      gamma1 = real(gamma);
      gamma2 = imag(gamma);
      // invert eq. (14) in Bacon et al. (2006)
      double Gamma1_1 = 0.5*(real(F) + real(G));
      double Gamma1_2 = 0.5*(imag(G) - imag(F));
      double Gamma2_1 = 0.5*(imag(F) + imag(G));
      double Gamma2_2 = 0.5*(real(F) - real(G));
      // 1/2 of eq. (5)
      D111 = -Gamma1_1 - 0.5*Gamma2_2;
      D121 = -0.5*Gamma2_1;
      D211 = -0.5*Gamma2_1;
      D221 = -0.5*Gamma2_2;
      D112 = -0.5*Gamma2_1;
      D122 = -0.5*Gamma2_2;
      D212 = -0.5*Gamma2_2;
      D222 = Gamma1_2 - 0.5*Gamma2_1;
      complex<data_t> F_(-D111-D221, -D121-D222);
      complex<data_t> G_(-D111+3*D221, -3*D121+D222);
    }
    /// Apply transformation to \p P.
    virtual void transform(Point<T>& P) const {
      // store temporarily
      T p0 = P(0), p1 = P(1);
      // apply eq. (3)
      P(0) = (1-kappa-gamma1)*p0 - gamma2*p1;
      P(1) = -gamma2*p0 + (1-kappa+gamma1)*p1;
      if (flex) {
	P(0) += D111*p0*p0 + (D112 + D121)*p0*p1 + D122*p1*p1;
	P(1) += D211*p0*p0 + (D212 + D221)*p0*p1 + D222*p1*p1;
      }
      if (this->next != NULL)
	this->next->transform(P);
    }
    /// Get a deep copy of \p this.
    virtual CoordinateTransformation<T>* clone() const {
      if (flex) {
	complex<data_t> F(-D111-D221, -D121-D222);
	complex<data_t> G(-D111+3*D221, -3*D121+D222);
	return clone_pointer(new LensingTransformation(kappa,complex<data_t>(gamma1,gamma2), F, G));
      }
      else
	return clone_pointer(new LensingTransformation(kappa,complex<data_t>(gamma1,gamma2)));
    }
    /// Get a deep copy of this inverse of \p this.
    virtual CoordinateTransformation<T>* getInverse() const {
      if (flex) {
	complex<data_t> F(D111+D221, D121+D222);
	complex<data_t> G(D111-3*D221, 3*D121-D222);
	return inverse_pointer(new LensingTransformation(-kappa,complex<data_t>(-gamma1,-gamma2), F, G));
      } else
	return inverse_pointer(new LensingTransformation(-kappa,complex<data_t>(-gamma1,-gamma2)));
    }
  private:
    bool flex;
    data_t kappa, gamma1, gamma2, D111, D112, D121, D122, D211, D212, D221, D222;
    
  };
} // end namespace
#endif
