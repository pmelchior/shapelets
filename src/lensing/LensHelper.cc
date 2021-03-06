#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/MathHelper.h"
#include <gsl/gsl_randist.h>

// All equations taken from Bartelmann & Schneider (2001)
namespace shapelens {
  
  void eps2chi(const std::complex<data_t>& eps, std::complex<data_t>& chi) {
    chi = 2.*eps/(1+pow2(abs(eps)));
  }
  void chi2eps(const std::complex<data_t>& chi, std::complex<data_t>& eps) {
    eps = chi/(1. + sqrt(complex<data_t>(1 - pow2(abs(chi)))));
  }
  void lensEps(const std::complex<data_t>& gamma, std::complex<data_t>& eps) {
    eps = (eps + gamma)/(1. + conj(gamma)*eps);
  }

  data_t sampleEllipticity(const gsl_rng * r, data_t sigma_e, complex<data_t>& eps, data_t limit) {
    complex<data_t> I(0,1);
    data_t e = 1;
    while (e > limit) {
      e = gsl_ran_rayleigh (r,sigma_e);
    }
    eps = e * exp(2*M_PI*I*gsl_rng_uniform(r));
    return gsl_ran_rayleigh_pdf (e, sigma_e);
  }

  std::complex<data_t> epsilon(const Moments& mo) {
    if (mo.getOrder() >= 2) {
      if (mo(2,0) == 0. && mo(1,1) == 0. && mo(0,2) == 0.)
	return std::complex<data_t>(0,0);
      else {
	std::complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
	e/= (std::complex<data_t>(mo(2,0) + mo(0,2)) + 2.*sqrt(std::complex<data_t>(mo(0,2)*mo(2,0) - mo(1,1)*mo(1,1))));
	return e;
      }
    } else
      return std::complex<data_t>(0,0);
  }

  std::complex<data_t> chi(const Moments& mo) {
    if (mo.getOrder() >= 2) {
      if (mo(2,0) == 0. && mo(1,1) == 0. && mo(0,2) == 0.)
	return std::complex<data_t>(0,0);
      else {
	std::complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
	e/= mo(2,0) + mo(0,2);
	return e;
      }
    } else
      return std::complex<data_t>(0,0);
  }

  // HOLICS equations from OMU 2007
  complex<data_t> zeta(const Moments& mo) {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4); // eq. 25
      complex<data_t> zeta(mo(3,0) + mo(1,2),    // eq. 26
			   mo(2,1) + mo(0,3));
      zeta /= xi;
      return zeta;
    } else
      return complex<data_t>(0,0);
  }

  complex<data_t> delta(const Moments& mo) {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4);  // eq. 25
      complex<data_t> delta(mo(3,0) - 3*mo(1,2),  // eq. 27
			    3*mo(2,1) - mo(0,3));
      delta /= xi;
      return delta;
    } else
      return complex<data_t>(0,0);
  }

  data_t trQ(const Moments& mo) {
    if (mo.getOrder() >= 2)
      return (mo(2,0)+mo(0,2))/mo(0,0);
    else
      return 0;
  }
  
  data_t R2(const Moments& g, const Moments& p) {
    if (g.getOrder() >= 2)
      return 1 - (trQ(p)/trQ(g));
    else
      return 0;
  } 

} // end namespace
