#ifndef SHAPELENS_DEIMOS_H
#define SHAPELENS_DEIMOS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"
#include "../frame/CoordinateTransformation.h"
#include "../utils/SQLiteDB.h"
#include <vector>
#include <bitset>
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

    /// Default constructor.
    DEIMOS();
    /// Constructor from filename.
    DEIMOS(std::string filename);
    /// Constructor for building moments up to order \p N from \p obj.
    DEIMOS (Object& obj, int N, int C, data_t scale, bool flexed = false);
    /// Save to a file.
    void save(std::string filename) const;
    /// Deconvolve \p obj from \p psf.
    void deconvolve(const DEIMOS& psf);
    /// Get std::complex ellipticity from mo.
    std::complex<data_t> epsilon() const;
    /// Get std::complex ellipticity from mo.
    std::complex<data_t> chi() const;
    /// Get first flexion distortion from mo.
    std::complex<data_t> zeta() const;
    /// Get second flexion distortion from mo.
    std::complex<data_t> delta() const;
    /// The ordered set of multipole moments.
    Moments mo;
    /// The ordered set of multipole moment errors.
    Moments mo_noise;
    /// Deweighting order.
    int C;
    /// Object id.
    unsigned long id;
    /// Width of the weighting function.
    data_t scale;
    /// Ellipticity of weighting function
    std::complex<data_t> eps;
    /// 2nd flexion of weighting function
    std::complex<data_t> G;
    
    friend class DEIMOSList;
  private:
    void match(Object& obj, int N);
    void deweight(const Moments& mo_w, int N);
    std::complex<data_t> epsilon_limited();
    void estimateErrors(const Object& obj, int N);
    bool flexed;
    NumMatrix<data_t> D;
  };

  /// Class for collections of DEIMOS instances.
  class DEIMOSList : public std::vector<boost::shared_ptr<DEIMOS> > {
  public:
    /// Default constructor.
    DEIMOSList ();
    /// Constructor from SQLite
    DEIMOSList (SQLiteDB& sql, std::string table, std::string where);
    /// Save list to table in sql.
    void save(SQLiteDB& sql, std::string table) const;
  };

} // end namespace
#endif
