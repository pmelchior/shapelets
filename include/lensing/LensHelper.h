#ifndef SHAPELENS_LENSHELPER
#define SHAPELENS_LENSHELPER
#include <complex>
#include <gsl/gsl_rng.h>
#include "../Typedef.h"
#include "../../include/frame/Moments.h"

namespace shapelens {
  
  /// Convert \f$\epsilon\f$ to \f$\chi\f$.
  void eps2chi(const std::complex<data_t>& eps, std::complex<data_t>& chi);
  /// Convert \f$\chi\f$ to \f$\epsilon\f$.
  void chi2eps(const std::complex<data_t>& chi, std::complex<data_t>& eps);
  /// Apply shear \f$\gamma\f$ to unlensed ellipticity \f$\epsilon\f$.
  void lensEps(const std::complex<data_t>& gamma, std::complex<data_t>& eps);
  /// Sample the intrinsic ellipiticity distribution.
  /// The distribution is described by a Rayleigh distribution with mode
  /// \f$\sigma_e\f$ (the dispersion in each components of the ellipticity),
  /// truncated at \p limit.\n
  /// Returns probability of finding the sampled value of \f$\epsilon\f$.
  data_t sampleEllipticity(const gsl_rng * r, data_t sigma_e, complex<data_t>& eps, data_t limit=0.95);
  /// Get complex ellipticity \f$\epsilon\f$ from Moments.
  std::complex<data_t> epsilon(const Moments& mo);
  /// Get complex ellipticity \f$\chi\f$ from Moments.
  std::complex<data_t> chi(const Moments& mo);
} // end namespace

#endif
