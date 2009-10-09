#include "../../include/lensing/KSB.h"

namespace shapelens {
  // helper function
  complex<data_t> operator*(const NumMatrix<data_t>& M, complex<data_t> d) {
    return complex<data_t> (M(0,0)*real(d)*M(0,1)*imag(d),M(1,0)*real(d)+M(1,1)*imag(d));
  }
    
  KSB::KSB(const Object& obj_) {
    // for changing obj.w, we need to have write permissions...
    Object& obj = const_cast<Object&>(obj_);
    M = obj.flux;
    if (obj.w.getType() == 0)
      throw std::invalid_argument("KSB: usage of flat weights is not permitted");
    Moment2 Q2(obj);
    trQ = __trQ(Q2);
    chi = __chi(Q2);

    // measure moments with W' (w.r.t r^2)
    obj.w.setDerivative(-1);
    Moment1 Q1_(obj);
    Moment2 Q2_(obj);
    Moment3 Q3_(obj);
    Moment4 Q4_(obj);
    trQ_ = __trQ(Q2_);
    a_ = __a(Q2_,Q4_);
    b_ = __b(Q1_,Q3_);
    c_ = __c(Q2_);
    d_ = __d(Q2_);
    e_ = __e(Q1_);
    f_ = __f(Q1_,Q3_);
    g_ = __g(Q4_);

    // measure moments with W'' (w.r.t r^2)
    obj.w.setDerivative(-2);
    Moment1 Q1__(obj);
    Moment2 Q2__(obj);
    Moment3 Q3__(obj);
    Moment4 Q4__(obj);
    a__ = __a(Q2__,Q4__);
    b__ = __b(Q1__,Q3__);
    c__ = __c(Q2__);
    f__ = __f(Q1__,Q3__);
    g__ = __g(Q4__);
    
    // reset obj.w
    obj.w.setDerivative(0);

    // compute P_sh
    P_sh.resize(2,2);
    P_sh(0,0) = 2*trQ + 2*a_;
    P_sh(0,1) = 4*b_;
    P_sh(1,0) = P_sh(0,1);
    P_sh(1,1) = 2*trQ + 8*c_;
    P_sh /= trQ;
    NumMatrix<data_t> e_sh(1,2);
    e_sh(0,0) = 2*real(chi) + 2*g_/trQ;
    e_sh(0,1) = 2*imag(chi) + f_/trQ;
    P_sh(0,0) -= real(chi)*e_sh(0,0);
    P_sh(0,1) -= real(chi)*e_sh(0,1);
    P_sh(1,0) -= imag(chi)*e_sh(0,0);
    P_sh(1,1) -= imag(chi)*e_sh(0,1);

    // compute P_sm
    P_sm.resize(2,2);
    P_sm(0,0) = M + 2*trQ_ + a__;
    P_sm(0,1) = 2*b__;
    P_sm(1,0) = P_sm(0,1);
    P_sm(1,1) = M + 2*trQ_ + 4*c__;
    P_sm /= trQ;
    NumMatrix<data_t> e_sm(1,2);
    e_sm(0,0) = (2*d_ + g__)/trQ;
    e_sm(0,1) = (2*e_ + 2*f__)/trQ;
    P_sm(0,0) -= real(chi)*e_sm(0,0);
    P_sm(0,1) -= real(chi)*e_sm(0,1);
    P_sm(1,0) -= imag(chi)*e_sm(0,0);
    P_sm(1,1) -= imag(chi)*e_sm(0,1);
  }

  complex<data_t> KSB::__chi(const Moment2& Q) const {
    return complex<data_t>(Q(0,0)-Q(1,1),2*Q(0,1))/trQ;
  }
  data_t KSB::__trQ(const Moment2& Q) const {
    return Q(0,0) + Q(1,1);
  }
  data_t KSB::__a(const Moment2& Q2, const Moment4& Q4) const {
    return Q4(0,0,0,0) - 2*Q2(0,0)*Q2(1,1) + Q4(1,1,1,1);
  }
  data_t KSB::__b(const Moment1& Q1, const Moment3& Q3) const {
    return Q3(0,0,0)*Q1(1) - Q1(0)*Q3(1,1,1);
  }
  data_t KSB::__c(const Moment2& Q) const {
    return Q(0,0)*Q(1,1);
  }
  data_t KSB::__d(const Moment2& Q) const {
    return Q(0,0) - Q(1,1);
  }
  data_t KSB::__e(const Moment1& Q) const {
    return 2*Q(0)*Q(1);
  }
  data_t KSB::__f(const Moment1& Q1, const Moment3& Q3) const {
    return Q3(0,0,0)*Q1(1) + Q1(0)*Q3(1,1,1);
  }
  data_t KSB::__g(const Moment4& Q) const {
    return Q(0,0,0,0) - Q(1,1,1,1);
  }

  // eq. (22)
  complex<data_t> KSB::__p(const KSB& star) const {
    return star.P_sm.invert() * star.chi;
  }
  // eq. (12)
  complex<data_t> KSB::gamma() const {
    return P_sh.invert()*chi;
  }
  // eqs. (23) and (27, no delta)
  complex<data_t> KSB::gamma(const KSB& psf) const {
    NumMatrix<data_t> P_gamma = P_sh - psf.P_sh*psf.P_sm.invert()*P_sm;
    complex<data_t> p = __p(psf);
    return P_gamma.invert()*(chi - P_sm*p);
  }
}
