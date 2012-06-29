#ifndef SHAPELENS_DEIMOSCIRCULAR_H
#define SHAPELENS_DEIMOSCIRCULAR_H

#include "DEIMOS.h"
#include "../utils/History.h"

namespace shapelens {
  /// Class for moment-based lensing and deconvolution with an 
  /// elliptical weight function.
  /// This is a simplification to the approach published by
  /// Melchior et al. (2011), which uses a circular instead of an elliptcal
  /// weight function. The optimization for centroid and scale are still done.
  ///
  /// \b NOTE: If ShapeLensConfig::USE_WCS is set to \p true, the moments mo
  /// and all moment-related quantities (foremost scale, centroid, epsilon) are
  /// measured in these units.
  class DEIMOSCircular : public DEIMOS {
  public:
    /// Default constructor.
    DEIMOSCircular();
    /// Constructor from filename.
    DEIMOSCircular(std::string filename);
    /// Constructor for computing moments up to order \p N from \p obj.
    DEIMOSCircular(const Object& obj, int N, int C, data_t scale);
    /// Constructor for computing deconvolved moments up to order \p N from \p obj.
    DEIMOSCircular (const Object& obj, const PSFMultiScale& psf, int N, int C, data_t scale);
    /// Save to a file.
    void save(std::string filename) const;
    /// Correction order.
    int C;
    /// Matching scale of the weight function.
    /// Corresponds to the size of the semi-minor axis [pixel]:
    data_t matching_scale;
    /// Actual width of the weight function [WCS units].
    /// Determined during the matching procedure as
    /// \f$s = s_m/\sqrt{1 + \epsilon^2 - 2|\epsilon|}\f$,
    /// where \f$s_m\f$ denotes DEIMOS::matching_scale and \f$\epsilon\f$ 
    /// the convolved ellipticity. If WCS units are used, \f$s\f$ is rescale
    /// accordingly. 
    data_t scale;
    /// Pixel scale [WCS unit/pixel].
    data_t scale_factor;
    /// Centroid of weight function.
    Point<data_t> centroid;
    /// S/N of moment measurement.
    std::map<data_t, data_t> SN;
    /// Resolution factor \f$R_2\f$.
    /// From Hirata et al. (2004), eq. 8: \f$R_2 = 1 - trQ_p / trQ_{g*}\f$, where \f$trQ\f$
    /// denotes the sum of the symmetric second moments (of either the PSF or the
    /// convolved galaxy).
    data_t R2;
    /// Matching and deconvolution flags.
    /// The flags are ordered as DEIMOS processes objects, 
    /// i.e. matching (first centroid then ellipticity) followed by
    /// deconvolution (if necessary):
    /// - <tt>flags[0]</tt>: centroid determination failed
    /// - <tt>flags[1]</tt>: non-sensical moments or ellipiticty
    /// - <tt>flags[2]</tt>: ellipticity matching failed
    /// - <tt>flags[3]</tt>: deconvolution results in non-sensical moments
    std::bitset<4> flags;

    static bool FIX_CENTROID;
    friend class DEIMOSForward;
  protected:
    void match(const Object& obj, Moments& mo_w);
    void deweight(const Moments& mo_w);
    void computeDeweightingMatrix(const Moments& mo_w);
    void computeCovariances(const Moments& mo_w);
    data_t computeSN(const Moments& mo_w);
    NumMatrix<data_t> D;
    History history;
  };
} // end namespace
#endif
