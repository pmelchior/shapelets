#ifndef SHAPELENS_DEIMOS_H
#define SHAPELENS_DEIMOS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"
#include "../frame/CoordinateTransformation.h"
#include "../utils/SQLiteDB.h"
#include "../utils/History.h"
#include <boost/shared_ptr.hpp>
#include <bitset>

namespace shapelens {
  /// Class for moment-based lensing and deconvolution.
  /// DEIMOS implements the method introduced by Melchior et al. (2011), 	
  /// a moment-based approach to weak lensing with an analytic deconvolution
  /// procedure.
  ///
  /// It uses an iterative matching algorithm to determine the optimal
  /// size, center, and ellipticity of a weight function to measure the
  /// convolved object. If required, a deconvolution can be performed with
  /// a PSF shape, measured with the same parameters as the object.
  ///
  /// \b NOTE: If ShapeLensConfig::USE_WCS is set to \p true, the moments mo
  /// and all moment-related quantities (foremost scale, centroid, epsilon) are
  /// measured in these units.
  class DEIMOS {
  public:
    /// DEIMOS weight function.
    /// DEIMOS uses an elliptical (and possibly G-flexed) Gaussian 
    /// weight function, \f$W(x) = \exp\Bigl(\frac{-x^{\prime 2}}{2s^s}\Bigr)\f$,
    /// where \f$x^\prime\f$ undergoes a suitable LensingTransformation.
    class DEIMOSWeightFunction : public GaussianWeightFunction {
    public:
      /// Constructor for ellipticial weight function.
      DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid, const std::complex<data_t>& eps);
      /// Constructor for elliptical and G-flexed weight function.
      DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid_, const std::complex<data_t>& eps, const std::complex<data_t>& G);
      /// Get value of weight function at position \p P.
      virtual data_t operator()(const Point<data_t>& P) const;
    protected:
      LensingTransformation T;
      const Point<data_t>& centroid;
    };

    /// DEIMOS mulit-scale PSF container.
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
    };

    /// Default constructor.
    DEIMOS();
    /// Constructor from filename.
    DEIMOS(std::string filename);
    /// Constructor for computing moments up to order \p N from \p obj.
    DEIMOS (const Object& obj, int N, int C, data_t scale, bool flexed = false);
    /// Constructor for computing deconvolved moments up to order \p N from \p obj.
    DEIMOS (const Object& obj, const PSFMultiScale& psf, int N, int C, data_t scale, bool flexed = false);
    /// Deconvolve from Moments of \p psf.
    /// This is only relevant for the contructor without PSF specification. 
    void deconvolve(const Moments& psf);
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
    /// Deweighting order.
    int C;
    /// Copy of Object::id.
    unsigned long id;
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
    /// Ellipticity of weight function.
    std::complex<data_t> eps;
    /// 2nd flexion of weight function.
    std::complex<data_t> G;
    /// S/N of moment measurement.
    std::map<data_t, data_t> SN;
    /// Matching and deconvolution flags.
    /// The flags are ordered as DEIMOS processes objects, 
    /// i.e. matching (first centroid then ellipticity) followed by
    /// deconvolution (if necessary):
    /// - <tt>flags[0]</tt>: centroid determination failed
    /// - <tt>flags[1]</tt>: non-sensical moments or ellipiticty
    /// - <tt>flags[2]</tt>: ellipticity matching failed
    /// - <tt>flags[3]</tt>: deconvolution results in non-sensical moments
    std::bitset<4> flags;

    friend class DEIMOSList;
  protected:
    void match(Object& obj);
    void deweight(const Moments& mo_w);
    void setNoiseImage(const Object& obj);
    void computeCovariances();
    data_t getEpsScale() const;
    bool flexed;
    NumMatrix<data_t> D;
    Object noise;
    History history;
  };

  /// Class for collections of DEIMOS instances.
  class DEIMOSList : public std::vector<boost::shared_ptr<DEIMOS> > {
  public:
    /// Default constructor.
    DEIMOSList ();
#ifdef HAS_SQLiteDB
    /// Constructor from SQLite
    DEIMOSList (SQLiteDB& sql, std::string table, std::string where);
    /// Save list to table in sql.
    void save(SQLiteDB& sql, std::string table) const;
#endif
  };

} // end namespace
#endif
