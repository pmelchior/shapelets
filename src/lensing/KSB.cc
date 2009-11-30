#include "../../include/lensing/KSB.h"
#include <gsl/gsl_multimin.h>

namespace shapelens {
  // helper function
  complex<data_t> operator*(const NumMatrix<data_t>& M, complex<data_t> d) {
    return complex<data_t> (M(0,0)*real(d)+M(0,1)*imag(d),M(1,0)*real(d)+M(1,1)*imag(d));
  }
  // helper class
  KSB::NumTensor::NumTensor() : M0(2,2), M1(2,2) {
  }
  data_t& KSB::NumTensor::operator()(bool i, bool j, bool k) {
    if (i)
      return M1(j,k);
    else
      return M0(j,k);
  }
  const data_t& KSB::NumTensor::operator()(bool i, bool j, bool k) const {
    if (i)
      return M1(j,k);
    else
      return M0(j,k);
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

    Moment6 Q6__(obj);
    lambda__ = __lambda(Q6__);
    omega__  = __omega(Q6__);
    sigma__  = __sigma(Q6__);
    rho__    = __rho(Q6__);
    ix__     = __ix(Q6__);
    delta__  = __delta(Q6__);

    // reset to zero-th derivative since obj is const
    obj.w.setDerivative(0);

    R(0,0,0) = (real(chi)/trQ)*(2*rho__ + 4*psi_) - 2*lambda__/trQ - 4*pi_/trQ;
    R(0,0,1) = (real(chi)/trQ)*(4*ix__  + 8*mu_ ) - 4*omega__/trQ;
    R(0,1,0) = R(0,0,1) - 8*nu_/trQ;
    R(0,1,1) = (real(chi)/trQ)*(8*delta__ + 16*Q4_(0,0,1,1)) - 8*sigma__/trQ;
    
    R(1,0,0) = (imag(chi)/trQ)*(2*rho__ + 4*psi_) - 4*omega__/trQ;
    R(1,0,1) = (imag(chi)/trQ)*(4*ix__ + 8*mu_) - 8*sigma__/trQ - 4*pi_/trQ;
    R(1,1,0) = (imag(chi)/trQ)*(4*ix__ + 8*mu_) - 8*sigma__/trQ; 
    R(1,1,1) = (imag(chi)/trQ)*(8*delta__ + 16*Q4_(0,0,1,1)) - 16*Q6__(0,0,0,1,1,1)/trQ - 8*nu_/trQ;

    //Compute P_sh 
    P_sh.resize(2,2);
    P_sh(0,0) = -2*real(chi)*pi_/trQ-2*real(chi)*real(chi)+2*psi_/trQ + 2 ;
    P_sh(0,1) = -4*real(chi)*nu_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,0) = -2*imag(chi)*pi_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,1) = -4*imag(chi)*nu_/trQ-2*imag(chi)*imag(chi)+8*Q4_(0,0,1,1)/trQ+2;
					 
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
  data_t KSB::__lambda(const Moment6& Q) const {
    return Q(0,0,0,0,0,0) - 3*Q(0,0,0,0,1,1) + 3*Q(0,0,1,1,1,1) - Q(1,1,1,1,1,1);
  }
  data_t KSB::__omega(const Moment6& Q) const {
    return Q(0,0,0,0,0,1) + Q(0,1,1,1,1,1) - 2*Q(0,0,0,1,1,1);
  }
  data_t KSB::__sigma(const Moment6& Q) const {
    return Q(0,0,0,0,1,1) - Q(0,0,1,1,1,1);
  }
  data_t KSB::__rho(const Moment6& Q) const {
    return Q(0,0,0,0,0,0) - Q(0,0,0,0,1,1) - Q(0,0,1,1,1,1) + Q(1,1,1,1,1,1);
  }
  data_t KSB::__ix(const Moment6& Q) const {
    return Q(0,0,0,0,0,1) - Q(0,1,1,1,1,1);
  }
  data_t KSB::__delta(const Moment6& Q) const {
    return Q(0,0,0,0,1,1) + Q(0,0,1,1,1,1);
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


  // *** KSB goes non-linear here ! ***

  complex<data_t> __chi_sh_nl(const KSB& ksb, const complex<data_t>& gamma) {
    return complex<data_t>(// 1st component
			   real(gamma)*(ksb.P_sh(0,0)+real(gamma)*ksb.R(0,0,0)+1*imag(gamma)*ksb.R(0,0,1))+
			   imag(gamma)*(ksb.P_sh(0,1)+1*real(gamma)*ksb.R(0,1,0)+imag(gamma)*ksb.R(0,1,1))+
			   2*real(ksb.chi)*(imag(gamma)*imag(gamma))-
			   2*real(gamma)*imag(gamma)*imag(ksb.chi),
			   // 2nd component
			   real(gamma)*(ksb.P_sh(1,0)+real(gamma)*ksb.R(1,0,0)+1*imag(gamma)*ksb.R(1,0,1))+
			   imag(gamma)*(ksb.P_sh(1,1)+1*real(gamma)*ksb.R(1,1,0)+imag(gamma)*ksb.R(1,1,1))+
			   2*imag(ksb.chi)*(real(gamma)*real(gamma))-
			   2*real(gamma)*imag(gamma)*real(ksb.chi));
  }
  
  complex<data_t> __chi_g(const KSB& ksb, const KSB& star, const complex<data_t>& gamma) {
    complex<data_t> chi_sh_star = __chi_sh_nl(star,gamma);
    return ksb.P_sm*star.P_sm.invert()*chi_sh_star;
  }
  
  complex<data_t> __chi_sm(const KSB& ksb, const KSB& star) {
    NumMatrix<data_t> P_sm_star_1 = star.P_sm.invert();
    return complex<data_t>(// 1st component
			   ksb.P_sm(0,0)*(P_sm_star_1(0,0)*real(star.chi) + P_sm_star_1(0,1)*imag(star.chi)) +
			   ksb.P_sm(0,1)*(P_sm_star_1(1,0)*real(star.chi) + P_sm_star_1(1,1)*imag(star.chi)),
			   // 2nd component
			   ksb.P_sm(1,0)*(P_sm_star_1(0,0)*real(star.chi) + P_sm_star_1(0,1)*imag(star.chi)) +
			   ksb.P_sm(1,1)*(P_sm_star_1(1,0)*real(star.chi) + P_sm_star_1(1,1)*imag(star.chi)));
  }

  // helper struct
  struct KSBTuple {
    const KSB& ksb1;
    const KSB& ksb2;
  };

  double deltaChi_sh(const gsl_vector* v, void* p) {
    KSBTuple* kt = reinterpret_cast<KSBTuple*>(p);
    complex<data_t>* gamma = reinterpret_cast<complex<data_t>*>(v->data);
    complex<data_t> chi_sh = __chi_sh_nl(kt->ksb1,*gamma);
    return abs(kt->ksb1.chi - chi_sh);
  }

  double deltaChi(const gsl_vector* v, void* p) {
    KSBTuple* kt = reinterpret_cast<KSBTuple*>(p);
    complex<data_t>* gamma = reinterpret_cast<complex<data_t>*>(v->data);
    complex<data_t> chi_sh = __chi_sh_nl(kt->ksb1,*gamma);
    complex<data_t> chi_sm = __chi_sm(kt->ksb1,kt->ksb2);
    complex<data_t> chi_g = __chi_g(kt->ksb1,kt->ksb2,*gamma);
    //std::cout << real(*gamma) << "\t" << imag(*gamma) << "\t" << abs(kt->ksb1.chi - (chi_sh + chi_sm - chi_g)) << std::endl;
    return abs(kt->ksb1.chi - (chi_sh + chi_sm - chi_g));
  }

  complex<data_t> __gamma_nl(data_t accuracy, const KSBTuple& kt, bool convolved) {
    // initial value from linear treatment
    complex<data_t> gamma0;
    if (convolved)
      gamma0 = kt.ksb1.gamma(kt.ksb2);
    else
      gamma0 = kt.ksb1.gamma();
    // accuracy relative w.r.t. gamma0
    accuracy*=abs(gamma0);
  
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x, *min;
    gsl_multimin_function F;
    int status;
    data_t size;

    // define initial vertex vector
    ss = gsl_vector_alloc (2);
    // define primary stepsize
    // quadratic dependence on gamma0 and in opposite direction
    gsl_vector_set (ss, 0, -0.5*real(gamma0)*fabs(real(gamma0)));
    gsl_vector_set (ss, 1, -0.5*imag(gamma0)*fabs(imag(gamma0)));

    // Starting point
    x = gsl_vector_alloc(2);
    gsl_vector_set (x, 0, real(gamma0));
    gsl_vector_set (x, 1, imag(gamma0));

    // define the function which should be minimized and its parameters
    if (convolved)
      F.f = &deltaChi;
    else
      F.f = &deltaChi_sh;
    F.n = 2;
    F.params = reinterpret_cast<void*>(const_cast<KSBTuple*>(&kt));
    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &F, x, ss);

  
    do {
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	break;

      // the accuracy comes in here:
      // when size is smaller than accuracy we have convergence
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size,accuracy);
    }
    while (status == GSL_CONTINUE);

    min = gsl_multimin_fminimizer_x(s);
    complex<data_t> gamma_nl(gsl_vector_get(min,0), gsl_vector_get(min, 1));
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return gamma_nl;
  }

  complex<data_t> KSB::gamma_nl(data_t accuracy) const {
    const KSBTuple kt = {*this,*this};
    return __gamma_nl(accuracy, kt, false); 
  }

//   data_t devChi(const KSB& ksb, const KSB& star, complex<data_t>& gamma) {
//     complex<data_t> chi_sh = __chi_sh_nl(ksb,gamma);
//     complex<data_t> chi_sm = __chi_sm(ksb,star);
//     complex<data_t> chi_g = __chi_g(ksb,star,gamma);
//     return abs(ksb.chi - (chi_sh + chi_sm - chi_g));
//   }

  complex<data_t> KSB::gamma_nl(const KSB& psf, data_t accuracy) const {
    const KSBTuple kt = {*this,psf};
    // complex<data_t> gamma;
//     for (data_t g1 = -1; g1 <= 1; g1+=0.005) {
//       real(gamma) = g1;
//       for (data_t g2 = -1; g2 <= 1; g2+=0.005) {
// 	imag(gamma) = g2;
// 	if (g2 >= -sqrt(1-g1*g1) && g2 <= sqrt(1-g1*g1))
// 	  std::cout << g1 << "\t" << g2 << "\t" << devChi(*this,psf,gamma) << std::endl;
// 	else
// 	  std::cout << g1 << "\t" << g2 << "\t" << 0 << std::endl;
//       }
//     }
    return __gamma_nl(accuracy, kt, true); 
  }
}
