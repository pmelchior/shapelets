#ifndef SHAPELENS_KSB_H
#define SHAPELENS_KSB_H

#include "../frame/Object.h"
#include "../frame/Moments.h"
#include <set>

namespace shapelens {
  /// Class for Kaiser, Squires & Broadhurst (1995) shear estimates.
  class KSB {
  public:
    /// Constructor.
    KSB();
    /// Constructor for shear estimates from \p obj.
    /// The measurement employs a circular Gaussian weight function with
    /// a width given by \p scale.
    KSB(const Object& obj, data_t scale);
    /// Resulting complex ellipticity \f$\chi\f$.
    std::complex<data_t> chi;
    /// S/N of the measurement.
    data_t SN;

    // shear estimators without PSF correction
    // for stellar shapes
    /// Original KSB shear estimator.
    std::complex<data_t> gamma() const;
    /// First-order shear estimator. 
    std::complex<data_t> gamma_first() const;
    /// Shear estimator with the \i trace-trick.
    std::complex<data_t> gammaTr() const;
    /// First-order shear estimator with the \i trace-trick.
    std::complex<data_t> gammaTr1() const;

    // shear estimators with PSF correction: for galaxies
    /// Original KSB shear estimator with PSF correction.
    std::complex<data_t> gamma(const KSB& psf) const;
    /// First-order shear estimator with PSF correction.
    std::complex<data_t> gamma_first(const KSB& psf) const;
    /// Shear estimator with the \i trace-trick and PSF correction.
    std::complex<data_t> gammaTr(const KSB& psf) const;
    /// First-order shear estimator with the \i trace-trick and PSF correction.
    std::complex<data_t> gammaTr1(const KSB& psf) const;

    // non-linear estimator
    /// Second-order shear estimator.
    std::complex<data_t> gamma_second(data_t accuracy = 0.001) const;
    /// Exact shear estimator.
    std::complex<data_t> gamma_exact(data_t accuracy = 0.0001) const;
    /// Non-linear solver for shear given ellipticity.
    std::complex<data_t> gamma_nl(data_t accuracy = 0.000001) const;
    /// Non-linear solver for shear given ellipticity with PSF correction.
    std::complex<data_t> gamma_nl(const KSB& psf, data_t accuracy = 0.0001) const;

  private:
    std::complex<data_t> __chi(const Moment2& Q) const;
    data_t __trQ(const Moment2& Q) const;
    data_t __psi(const Moment4& Q) const;
    data_t __mu(const Moment4& Q) const;
    data_t __nu(const Moment4& Q) const;
    data_t __pi(const Moment4& Q) const;
    data_t __lambda(const Moment6& Q) const;
    data_t __omega(const Moment6& Q) const;
    data_t __sigma(const Moment6& Q) const;
    data_t __rho(const Moment6& Q) const;
    data_t __ix(const Moment6& Q) const;
    data_t __delta(const Moment6& Q) const;
    
    data_t ___a1(const Moment8& Q) const;
    data_t ___a2(const Moment8& Q) const;
    data_t ___a3(const Moment8& Q) const;
    data_t ___a4(const Moment8& Q) const;
    data_t ___a5(const Moment8& Q) const;
    
    std::complex<data_t> __p(const KSB& star) const;

  public:
    // various moment combinations measured with different derivatives 
    // of the weight function
    data_t trQ, trQ_, M,
      mu_, mu__, psi_, psi__, pi_, pi__, nu_, nu__, lambda__, omega__, sigma__,rho__,ix__,delta__, a1___,a2___,a3___,a4___,a5___;
    NumMatrix<data_t> P_sh, P1,P2,P_sm, e_sh, e_sm,K,B,Pa,Pb,Pc;

    // helper class for tensor of rank 3:
    // glue two NumMatrices togehter...
    class NumTensor {
    public:
      NumTensor();
      data_t& operator()(bool i, bool j, bool k);
      const data_t& operator()(bool i, bool j, bool k) const;
    private:
      NumMatrix<data_t> M0, M1;
    };
    NumTensor R,R1,U1,U2;
    
  };

  /// Practical implementation of KSB.
  /// Centroid position and optimum weight function scale are determined
  /// iteratively by requiring that the dipole moment of the weighted
  /// brightness distribution vanishes and that the measurement S/N is
  /// maximized.
  class KSBIterative : public KSB {
  public:
    /// Constructor from an Object and a list of scales.
    KSBIterative(const Object& obj, const std::set<data_t>& scales);
    /// Scale of the maximum S/N.
    data_t scale;
    /// Centroid measued at scale.
    Point<data_t> centroid;
  private:
    void findCentroid(Object& obj, data_t scale) const;
  };

} // end namespace shapelens

#endif
