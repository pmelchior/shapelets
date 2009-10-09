#include "../../include/lensing/HOLICS.h"
#include "../../include/frame/WeightFunction.h"

namespace shapelens {
  HOLICS::HOLICS(const Object& obj_) {
    // for changing obj.w, we need to have write permissions...
    Object& obj = const_cast<Object&>(obj_);
    if (obj.w.getType() == 0)
      throw std::invalid_argument("HOLICS: usage of flat weights is not permitted");
    sigma = obj.w.getScale();

    M = obj.flux;
    Moment2 Q2(obj);
    trQ = __trQ(Q2);
    Moment4 Q4(obj);
    xi = __xi(Q4);
    // zeta and delta need xi
    Moment3 Q3(obj);
    zeta = __zeta(Q3);
    delta = __delta(Q3);
    
    // measure moments with W'
    obj.w.setDerivative(1);
    Moment0 Q0_(obj);
    M_ = Q0_(0);
    Moment2 Q2_(obj);
    trQ_ = __trQ(Q2_);
    Moment4 Q4_(obj);
    xi_ = __xi(Q4_);
    Moment6 Q6_(obj);
    v0_ = __v0(Q6_);
    
    // measure moments with W''
    obj.w.setDerivative(2);
    Moment2 Q__(obj);
    trQ__ = __trQ(Q__);
    Moment4 Q4__(obj);
    xi__ = __xi(Q4__);

    // measure moments with W'''
    obj.w.setDerivative(3);
    Moment4 Q4___(obj);
    xi___ = __xi(Q4___);
    Moment6 Q6___(obj);
    v0___ = __v0(Q6___);
    // reset obj.w
    obj.w.setDerivative(0);

    // form high-order quantities...
    // eq. (C27)
    C_0_zeta = 9./4 + 3./(4*sigma*sigma)*v0_/xi;
    C_Delta_zeta = -2*trQ/xi - xi_/(sigma*sigma*xi);
    C_0_delta = 3./4 + v0_/(4*sigma*sigma*xi);
    // eq. (B8)
    Delta_0_L = 1.5 * trQ/M + (3*xi_)/(4*sigma*sigma*M);
    Delta_0_L /= (1 + trQ_/(sigma*sigma*M));
    // eq. (C18)
    P_0_zeta =  1./xi *( M + 5*trQ_/(sigma*sigma) + 7*xi__/(2*gsl_pow_4(sigma)) + v0___/(2*gsl_pow_6(sigma)));
    P_0_D = M_/(M*sigma*sigma) + 2*trQ__/(M*gsl_pow_4(sigma)) + xi___/(2*M*gsl_pow_6(sigma));
    P_Delta_D = 1 + trQ_/(M*sigma*sigma);
    P_Delta_zeta = 2*trQ/xi + xi_/(xi*sigma*sigma);
    P_0_delta = 1./xi *( M + 3*trQ_/(sigma*sigma) + 3*xi__/(2*gsl_pow_4(sigma)) + v0___/(6*gsl_pow_6(sigma)));
  }

  // eq. (20)
  complex<data_t> HOLICS::__zeta(const Moment3& Q) const {
    return complex<data_t>(Q(0,0,0) + Q(0,1,1), Q(0,0,1) + Q(1,1,1))/xi;
  }
  // eq. (20)
  complex<data_t> HOLICS::__delta(const Moment3& Q) const {
    return complex<data_t>(Q(0,0,0) - 3*Q(0,1,1), 3*Q(0,0,1) - Q(1,1,1))/xi;
  }
  // trace of eq. (12)
  data_t HOLICS::__trQ(const Moment2& Q) const {
    return Q(0,0) + Q(1,1);
  }
  // eq. (21)
  data_t HOLICS::__xi(const Moment4& Q) const {
    return Q(0,0,0,0) + 2*Q(0,0,1,1) + Q(1,1,1,1);
  }
  // eq. (A6)
  data_t HOLICS::__v0(const Moment6& Q) const {
    return Q(0,0,0,0,0,0) + 3*Q(0,0,1,1,1,1) + 3*Q(0,0,0,0,1,1) + Q(1,1,1,1,1,1);
  }

  
  // eq. (34)
  complex<data_t> HOLICS::F() const {
    return zeta/(9./4 + 3*v0_/(4*xi*sigma*sigma) - 
		 (2*trQ/xi + xi_/(xi*sigma*sigma))*Delta_0_L);
  }
  // eq. (35)
  complex<data_t> HOLICS::G() const {
    return delta/(3./4 + v0_/(4*xi*sigma*sigma));
  }
  // eqs. (52) & (78)
  complex<data_t> HOLICS::F(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> zeta_iso = g.zeta - 
      (g.P_0_zeta - s.P_0_D/s.P_Delta_D*s.P_Delta_zeta)*__zeta_q(s);
    return zeta_iso/((g.C_0_zeta + g.C_Delta_zeta*g.Delta_0_L) -
		     (s.C_0_zeta + s.C_Delta_zeta*s.Delta_0_L)*
		     (s.P_0_zeta - s.P_0_D*s.P_Delta_zeta/s.P_Delta_D)*
		     (g.P_0_zeta - g.P_0_D*g.P_Delta_zeta/g.P_Delta_D));
  }
  // eqs. (53) & (79)
  complex<data_t> HOLICS::G(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> delta_iso = g.delta - g.P_0_delta*__delta_q(s);
    return delta_iso/(g.C_0_delta - 
		      (s.C_0_delta/s.P_0_delta)*g.P_0_delta);
  }
  // eq. (56)
  complex<data_t> HOLICS::__zeta_q(const HOLICS& s) const {
    return s.zeta/(s.P_0_zeta - 
		   s.P_0_D/s.P_Delta_D*s.P_Delta_zeta);
  }
  // eq. (57)
  complex<data_t> HOLICS::__delta_q(const HOLICS& s) const {
    return s.delta/s.P_0_delta;
  }


  
  

} // end namespace
