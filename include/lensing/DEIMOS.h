#ifndef SHAPELENS_DEIMOS_H
#define SHAPELENS_DEIMOS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"
#include <bitset>
#include <map>

namespace shapelens {
  /// Base class for moment-based lensing and deconvolution.
  /// DEIMOS implements the method introduced by Melchior et al. (2011), 	
  /// a moment-based approach to weak lensing with an analytic deconvolution
  /// procedure.
  ///
  /// \b NOTE: If ShapeLensConfig::USE_WCS is set to \p true, the moments mo
  /// and all moment-related quantities (foremost scale, centroid, epsilon) are
  /// measured in these units.
  class DEIMOS {
  public:
    /// Default constructor.
    DEIMOS(int N = 0);
    /// Constructor from filename.
    DEIMOS(std::string filename);
    /// Deconvolve from Moments of \p psf.
    void deconvolve(const Moments& psf);
    /// Convolve with Moments of \p psf.
    void convolve(const Moments& psf);
    /// Save to a file.
    void save(std::string filename) const;
    /// Get std::complex ellipticity from mo.
    std::complex<data_t> epsilon() const;
    /// Get std::complex ellipticity from mo.
    std::complex<data_t> chi() const;
    /// Get first flexion distortion from mo.
    std::complex<data_t> zeta() const;
    /// Get second flexion distortion from mo.
    std::complex<data_t> delta() const;
    /// Get marginalized moment errors.
    Moments getMomentErrors() const;
    /// Check if moments are sensical (return 1 if they are not).
    bool flagMoments(const Moments& M) const;
    /// Ordered set of multipole moments.
    Moments mo;
    /// Covariance matrix of mo.
    NumMatrix<data_t> S;
    /// Moment order.
    int N;
  protected:
    void setNoiseImage(const Object& obj);
    Object noise;

  public:
    /// DEIMOS multi-scale PSF container.
    class PSFMultiScale : public std::map<data_t, Moments> {
    public:
      /// Insert Moments \p mo measured with \p scale.
      void insert(data_t scale, const Moments& mo);
      /// Get moments measured with \p scale.
      /// Throws <tt>std::invalid_argument</tt> if \p scale is not available.
      const Moments& getAtScale(data_t scale) const;
      /// Get next smaller scale.
      /// Throws <tt>std::invalid_argument</tt> if \p scale is not available,
      /// and <tt>std::runtime_error</tt> if \p scale is the smallest scale. 
      data_t getScaleSmallerThan(data_t scale) const;
      /// Get minimum scale available.
      data_t getMinimumScale() const;
      /// Get maximum scale available.
      data_t getMaximumScale() const;
      /// Get available scale closest to \p scale.
      data_t getScaleClosestTo(data_t scale) const;
    };
  };
} // end namespace
#endif
