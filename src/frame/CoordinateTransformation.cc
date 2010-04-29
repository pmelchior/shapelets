#include "../../include/frame/CoordinateTransformation.h"

namespace shapelens {

  CoordinateTransformation::~CoordinateTransformation() {
  }

  void CoordinateTransformation::operator*=(const CoordinateTransformation& C) {
    stack.push_back(boost::shared_ptr<CoordinateTransformation>(C.clone()));
  }

  void CoordinateTransformation::stack_transform(Point<data_t>& P) const {
    for(std::list<boost::shared_ptr<CoordinateTransformation> >::const_iterator iter = stack.begin(); iter != stack.end(); iter++)
      (*iter)->transform(P);
  }

  void CoordinateTransformation::stack_inverse_transform(Point<data_t>& P) const {
    for(std::list<boost::shared_ptr<CoordinateTransformation> >::const_reverse_iterator iter = stack.rbegin(); iter != stack.rend(); iter++)
      (*iter)->inverse_transform(P);
  }
    
  /// NullTransformation ...
  NullTransformation::NullTransformation() {
  }
  void NullTransformation::transform(Point<data_t>& P) const {
    // do nothing, just go through the remainder of the stack
    stack_transform(P);
  }
  void NullTransformation::inverse_transform(Point<data_t>& P) const {
    // do nothing, just go through the remainder of the stack
    stack_inverse_transform(P);
  }
  boost::shared_ptr<CoordinateTransformation> NullTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new NullTransformation(*this));
  }

  /// ScalarTransformation ...
  ScalarTransformation::ScalarTransformation(data_t scale): s(scale) {}
    
  void ScalarTransformation::transform(Point<data_t>& P) const {
    P *= s;
    stack_transform(P);
  }
  void ScalarTransformation::inverse_transform(Point<data_t>& P) const {
    stack_inverse_transform(P);
    P /= s; // inverse: this trafo comes latest
  }
  boost::shared_ptr<CoordinateTransformation> ScalarTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new ScalarTransformation(*this));
  }

  /// ShiftTransformation ...
  ShiftTransformation::ShiftTransformation(const Point<data_t>& dP_): dP(dP_) {}
    
  void ShiftTransformation::transform(Point<data_t>& P) const {
    P += dP;
    stack_transform(P);
  }
  void ShiftTransformation::inverse_transform(Point<data_t>& P) const {
    stack_inverse_transform(P);
    P -= dP;
  }
  boost::shared_ptr<CoordinateTransformation> ShiftTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new ShiftTransformation(*this));
  }

  /// LinearTransformation ...
  LinearTransformation::LinearTransformation(const NumMatrix<data_t>& M_): M(M_), M_1(M_.invert()) {}
    
  void LinearTransformation::transform(Point<data_t>& P) const {
    data_t p0 = P(0);
    P(0) = M(0,0)*P(0) + M(0,1)*P(1);
    P(1) = M(1,0)* p0  + M(1,1)*P(1);
    stack_transform(P);
  }
  void LinearTransformation::inverse_transform(Point<data_t>& P) const {
    stack_inverse_transform(P);
    data_t p0 = P(0);
    P(0) = M_1(0,0)*P(0) + M_1(0,1)*P(1);
    P(1) = M_1(1,0)* p0  + M_1(1,1)*P(1);
  }
  boost::shared_ptr<CoordinateTransformation> LinearTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new LinearTransformation(*this));
  }

  /// LensingTransformation ...
  LensingTransformation::LensingTransformation(data_t kappa_, complex<data_t> gamma) : 
    kappa(kappa_), flex(false) {
    gamma1 = real(gamma);
    gamma2 = imag(gamma);
  }
  LensingTransformation::LensingTransformation(data_t kappa_, complex<data_t> gamma, complex<data_t> F, complex<data_t> G) : 
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

  void LensingTransformation::transform(Point<data_t>& P) const {
    // store temporarily
    data_t p0 = P(0), p1 = P(1);
    // apply eq. (3)
    P(0) = (1-kappa-gamma1)*p0 - gamma2*p1;
    P(1) = -gamma2*p0 + (1-kappa+gamma1)*p1;
    if (flex) {
      P(0) += D111*p0*p0 + (D112 + D121)*p0*p1 + D122*p1*p1;
      P(1) += D211*p0*p0 + (D212 + D221)*p0*p1 + D222*p1*p1;
    }
      
    stack_transform(P);
  }

  void LensingTransformation::inverse_transform(Point<data_t>& P) const {
    stack_inverse_transform(P);

    // store temporarily
    data_t p0 = P(0), p1 = P(1);
    // apply inverse eq. (3): gamma -> -gamma etc.
    P(0) = (1+kappa+gamma1)*p0 + gamma2*p1;
    P(1) = gamma2*p0 + (1+kappa-gamma1)*p1;
    if (flex) {
      P(0) -= D111*p0*p0 + (D112 + D121)*p0*p1 + D122*p1*p1;
      P(1) -= D211*p0*p0 + (D212 + D221)*p0*p1 + D222*p1*p1;
    }
  }

  boost::shared_ptr<CoordinateTransformation> LensingTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new LensingTransformation(*this));
  }

} // end namespace
