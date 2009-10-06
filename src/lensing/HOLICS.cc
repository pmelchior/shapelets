#include "../../include/lensing/HOLICS.h"

namespace shapelens {
  HOLICS::HOLICS(const Object& obj_, data_t sigma_) : sigma(sigma_) {
    // copy for not changing original with modified weights
    Object obj(obj_);
    if (obj.weight.size()==0)
      obj.weight.resize(obj.size());
    // weight with W and recompute flux and centroid
    W w(sigma,obj.centroid); // centroid still unweighted.
    for (long i=0; i < obj.size(); i++) {
      if (obj_.weight.size())
	obj.weight(i) *= w.f(obj.grid(i));
      else
	obj.weight(i) = w.f(obj.grid(i));  
    }
    // FIXME: check if centroid moves a lot here, see below...
    obj.computeFluxCentroid();
    w = W(sigma,obj.centroid); // weighted centroid
    
    // multiply obj.weight with W ...
    for (long i=0; i < obj.size(); i++) {
      if (obj_.weight.size())
	obj.weight(i) = obj_.weight(i) * w.f(obj.grid(i));
      else
	obj.weight(i) = w.f(obj.grid(i));
    }
    // FIXME: recomputing centroid is not necessary if it does not change
    obj.computeFluxCentroid();
    Moment0 M(obj);
    Moment2 Q2(obj);
    trQ = __trQ(Q2);
    Moment4 Q4(obj);
    xi = __xi(Q4);
    // zeta and delta need xi
    Moment3 Q3(obj);
    zeta = __zeta(Q3);
    delta = __delta(Q3);
    
    // reset obj.weight and mutiply with W' ...
    for (long i=0; i < obj.size(); i++) {
      if (obj_.weight.size())
	obj.weight(i) = obj_.weight(i) * w.f_(obj.grid(i));
      else
	obj.weight(i) = w.f_(obj.grid(i));
    }
    Moment0 M_(obj);
    Moment2 Q2_(obj);
    trQ_ = __trQ(Q2_);
    Moment4 Q4_(obj);
    xi_ = __xi(Q4_);
    Moment6 Q6_(obj);
    v0_ = __v0(Q6_);
    
    // reset obj.weight and mutiply with W'' ...
    for (long i=0; i < obj.size(); i++) {
      if (obj_.weight.size())
	obj.weight(i) = obj_.weight(i) * w.f__(obj.grid(i));
      else
	obj.weight(i) = w.f__(obj.grid(i));
    }
    Moment2 Q__(obj);
    trQ__ = __trQ(Q__);
    Moment4 Q4__(obj);
    xi__ = __xi(Q4__);

    // reset obj.weight and mutiply with W''' ...
    for (long i=0; i < obj.size(); i++) {
      if (obj_.weight.size())
	obj.weight(i) = obj_.weight(i) * w.f___(obj.grid(i));
      else
	obj.weight(i) = w.f___(obj.grid(i));
    }
    Moment4 Q4___(obj);
    xi___ = __xi(Q4___);
    Moment6 Q6___(obj);
    v0___ = __v0(Q6___);
    
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
  // eq. (C27)
  data_t HOLICS::__C_0_zeta() const {
    return 9./4 + 3./(4*sigma*sigma)*v0_/xi;
  }
  data_t HOLICS::__C_Delta_zeta() const {
    return -2*trQ/xi - xi_/(sigma*sigma*xi);
  }
  data_t HOLICS::__C_0_delta() const {
    return 3./4 + v0_/(4*sigma*sigma*xi);
  }
  // eq. (B8)
  data_t HOLICS::__Delta_0_L() const {
    data_t num = 1.5 * trQ/M + (3*xi_)/(4*sigma*sigma*M);
    return num/(1 + trQ_/(sigma*sigma*M));
  }
  // eq. (C18)
  data_t HOLICS::__P_0_zeta() const {
    return 1./xi *( M + 5*trQ_/(sigma*sigma) + 7*xi__/(2*gsl_pow_4(sigma)) + v0___/(2*gsl_pow_6(sigma)));
  }
  data_t HOLICS::__P_0_D() const {
    return M_/(M*sigma*sigma) + 2*trQ__/(M*gsl_pow_4(sigma)) + xi___/(2*M*gsl_pow_6(sigma));
  }
  data_t HOLICS::__P_Delta_D() const {
    return 1 + trQ_/(M*sigma*sigma);
  }
  data_t HOLICS::__P_Delta_zeta() const {
    return 2*trQ/xi + xi_/(xi*sigma*sigma);
  }
  data_t HOLICS::__P_0_delta() const {
    return 1./xi *( M + 3*trQ_/(sigma*sigma) + 3*xi__/(2*gsl_pow_4(sigma)) + v0___/(6*gsl_pow_6(sigma)));
  }
  
  // eq. (34)
  complex<data_t> HOLICS::F() const {
    return zeta/(9./4 + 3*v0_/(4*xi*sigma*sigma) - 
		 (2*trQ/xi + xi_/(xi*sigma*sigma))*__Delta_0_L());
  }
  // eq. (35)
  complex<data_t> HOLICS::G() const {
    return delta/(3./4 + v0_/(4*xi*sigma*sigma));
  }
  // eqs. (52) & (78)
  complex<data_t> HOLICS::F(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> zeta_iso = g.zeta - 
      (g.__P_0_zeta() - s.__P_0_D()/s.__P_Delta_D()*s.__P_Delta_zeta())*__zeta_q(s);
    return zeta_iso/((g.__C_0_zeta()+g.__C_Delta_zeta()*g.__Delta_0_L()) -
		     (s.__C_0_zeta()+s.__C_Delta_zeta()*s.__Delta_0_L())*
		     (s.__P_0_zeta()-s.__P_0_D()*s.__P_Delta_zeta()/s.__P_Delta_D())*
		     (g.__P_0_zeta()-g.__P_0_D()*g.__P_Delta_zeta()/g.__P_Delta_D()));
  }
  // eqs. (53) & (79)
  complex<data_t> HOLICS::G(const HOLICS& s) const {
    const HOLICS& g = *this;
    complex<data_t> delta_iso = g.delta - g.__P_0_delta()*__delta_q(s);
    return delta_iso/(g.__C_0_delta() - 
		      (s.__C_0_delta()/s.__P_0_delta())*g.__P_0_delta());
  }
  // eq. (56)
  complex<data_t> HOLICS::__zeta_q(const HOLICS& s) const {
    return s.zeta/(s.__P_0_zeta() - 
		   s.__P_0_D()/s.__P_Delta_D()*s.__P_Delta_zeta());
  }
  // eq. (57)
  complex<data_t> HOLICS::__delta_q(const HOLICS& s) const {
    return s.delta/s.__P_0_delta();
  }


  // Gaussian weight function
  HOLICS::W::W(data_t sigma, const Point<data_t>& centroid) :
    sigma2(sigma*sigma), C(centroid) {
  }
  // W(P)
  data_t HOLICS::W::f(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return f(r);
  }
  // W'(P)
  data_t HOLICS::W::f_(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return f(r)*(-r/sigma2);
  }
  // W''(P)
  data_t HOLICS::W::f__(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return f(r)*(gsl_pow_2(r/sigma2) - 1./sigma2);
  }
  // W'''(P)
  data_t HOLICS::W::f___(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return f(r)*(gsl_pow_3(-r/sigma2) + 3*r/(sigma2*sigma2));
  }
  data_t HOLICS::W::f(data_t r) const {
    return exp(-r*r/(2*sigma2));
  }
  

} // end namespace
