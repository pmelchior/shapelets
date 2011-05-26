#ifndef SHAPELENS_DEIMOS_H
#define SHAPELENS_DEIMOS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"
#include "../frame/CoordinateTransformation.h"
#include "../utils/SQLiteDB.h"
#include "../utils/History.h"
#include <boost/shared_ptr.hpp>

namespace shapelens {
  /// Class for moment-based lensing and deconvolution.
  class DEIMOS {
  public:
    class DEIMOSWeightFunction : public GaussianWeightFunction {
    public:
      DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid, const std::complex<data_t>& eps);
      DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid_, const std::complex<data_t>& eps, const std::complex<data_t>& G);
      virtual data_t operator()(const Point<data_t>& P) const;
    private:
      LensingTransformation T;
      const Point<data_t>& centroid;
    };

    class PSFMultiScale : private std::map<data_t, Moments> {
    public:
      void insert(data_t scale, const Moments& mo);
      const Moments& getAtScale(data_t scale) const;
      data_t getScaleSmallerThan(data_t scale) const;
      data_t getMinimumScale() const;
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
    /// Ordered set of multipole moments.
    Moments mo;
    /// Covariance matrix of mo.
    NumMatrix<data_t> S;
    /// Moment order.
    int N;
    /// Deweighting order.
    int C;
    /// Object id.
    unsigned long id;
    /// Final width of the weight function.
    data_t scale;
    /// Centroid of weight function.
    Point<data_t> centroid;
    /// Ellipticity of weight function
    std::complex<data_t> eps;
    /// 2nd flexion of weight function
    std::complex<data_t> G;
    
    friend class DEIMOSList;
  protected:
    void match(const Object& obj, data_t min_scale);
    void deweight(const Moments& mo_w);
    void deconvolve(const Moments& psf);
    std::complex<data_t> epsilon_limited();
    void setNoiseImage(const Object& obj);
    void computeCovariances();
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
