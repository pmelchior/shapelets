#include "../../include/lensing/KSB.h"

namespace shapelens {
  // helper function
  complex<data_t> operator*(const NumMatrix<data_t>& M, complex<data_t> d) {
    return complex<data_t> (M(0,0)*real(d)+M(0,1)*imag(d),M(1,0)*real(d)+M(1,1)*imag(d));
  }
    
  KSB::KSB(const Object& obj_) {
    // for changing obj.w, we need to have write permissions...
    Object& obj = const_cast<Object&>(obj_);
    Moment0 M0(obj);
    M = M0(0);
    if (obj.w.getType() == 0)
      throw std::invalid_argument("KSB: usage of flat weights is not permitted");
    Moment2 Q2(obj);
    trQ = __trQ(Q2);
    chi = __chi(Q2);

    // measure moments with W' (w.r.t r^2)
    obj.w.setDerivative(-1);
    Moment2 Q2_(obj);
    Moment4 Q4_(obj);
    trQ_ = __trQ(Q2_);
    psi_ = __psi(Q4_);
    mu_ = __mu(Q4_);
    pi_ = __pi(Q4_);
    nu_ = __nu(Q4_);

    // measure moments with W'' (w.r.t r^2)
    obj.w.setDerivative(-2);
    Moment4 Q4__(obj);
    mu__ = __mu(Q4__);
    psi__ = __psi(Q4__);
    pi__ = __pi(Q4__);
    nu__ = __nu(Q4__);
    
    // reset obj.w
    obj.w.setDerivative(0);

    // compute P_sh
    P_sh.resize(2,2);
    P_sh(0,0) = 2*trQ + 2*psi_;
    P_sh(0,1) = 4*mu_;
    P_sh(1,0) = P_sh(0,1);
    P_sh(1,1) = 2*trQ + 8*Q4_(0,0,1,1);
    P_sh /= trQ;
    e_sh.resize(1,2);
    e_sh(0,0) = 2*real(chi) + 2*pi_/trQ;//2*real(chi) -2*real(chi)*real(chi) + 2*pi_/trQ;
    e_sh(0,1) = 2*imag(chi) + 4*nu_/trQ;//2*imag(chi) -2*imag(chi)*imag(chi) + 4*nu_/trQ;
    P_sh(0,0) -= real(chi)*e_sh(0,0);
    P_sh(0,1) -= real(chi)*e_sh(0,1);
    P_sh(1,0) -= imag(chi)*e_sh(0,0);
    P_sh(1,1) -= imag(chi)*e_sh(0,1);

    // compute P_sm
    P_sm.resize(2,2);
    P_sm(0,0) = M + 2*trQ_ + psi__;
    P_sm(0,1) = 2*mu__;
    P_sm(1,0) = P_sm(0,1);
    P_sm(1,1) = M + 2*trQ_ + 4*Q4__(0,0,1,1);
    P_sm /= trQ;
    e_sm.resize(1,2);
    e_sm(0,0) = (2*(Q2_(0,0)-Q2_(1,1)) + pi__)/trQ;
    e_sm(0,1) = (4*Q2_(0,1) + 2*nu__)/trQ;
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
  data_t KSB::__psi(const Moment4& Q) const {
    return Q(0,0,0,0) - 2*Q(0,0,1,1) + Q(1,1,1,1);
  }
  data_t KSB::__mu(const Moment4& Q) const {
    return Q(0,0,0,1) - Q(0,1,1,1);
  }
  data_t KSB::__nu(const Moment4& Q) const {
    return Q(0,0,0,1) + Q(0,1,1,1);
  }
  data_t KSB::__pi(const Moment4& Q) const {
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
