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
} // end namespace
