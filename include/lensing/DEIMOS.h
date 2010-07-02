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
      DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid, const complex<data_t>& eps);
      virtual data_t operator()(const Point<data_t>& P) const;
      //void setEpsilon(const complex<data_t>& eps);
      //const complex<data_t>& getEpsilon() const;
    private:
      LensingTransformation T;
      //GaussianWeightFunction G;
    };

    /// Default constructor.
    DEIMOS();
    /// Constructor from filename.
    DEIMOS(std::string filename);
    /// Constructor for building moments up to order \p N from \p obj.
    DEIMOS (Object& obj, int N, int C, data_t scale);
    /// Save to a file.
    void save(std::string filename) const;
    /// Correct the moments for the application of \p obj.w.
    void deweight();
    /// Deconvolve \p obj from \p psf.
    void deconvolve(const DEIMOS& psf);
    /// Get complex ellipticity from mo.
    complex<data_t> epsilon() const;
    /// Get complex ellipticity from mo.
    complex<data_t> chi() const;
    /// Get first flexion distortion from mo.
    complex<data_t> zeta() const;
    /// Get second flexion distortion from mo.
    complex<data_t> delta() const;
    /// The ordered set of multipole moments.
    MomentsOrdered mo;
    /// Deweighting order.
    int C;
    /// The transformation flags.
    /// If the zeroth bit is set, the moments are deweighted, if the first
    /// one is set, they are deconvolved.
    std::bitset<2> flags;
    /// Object id.
    unsigned long id;
    /// Width of the weighting function.
    data_t scale;
    /// Ellipticity of weighting function
    complex<data_t> eps;

    friend class DEIMOSList;
  private:
    void focus(Object& obj, int N);

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
