#ifndef SHAPELENS_SOURCEMODEL_H
#define SHAPELENS_SOURCEMODEL_H

#include <boost/shared_ptr.hpp>
#include <vector>
#include "../frame/Shapes.h"
#include "../frame/Object.h"
#include "../frame/Catalog.h"
#include "../shapelets/ShapeletObject.h"

namespace shapelens {
/// Base class for idealized source models.
/// A galaxy model is a idealized representation of a two-dimensional shape.
/// It's main advantage: It can be sampled at any resolution.

class SourceModel {
public:
  /// Destructor.
  virtual ~SourceModel();
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const = 0;
  /// Get rectangluar support area of the model.
  const Rectangle<data_t>& getSupport() const;
  /// Get centroid of model.
  const Point<data_t>& getCentroid() const;
  /// Get total flux of model.
  virtual data_t getFlux() const = 0;
  /// Get type of model.
  virtual char getModelType() const = 0;
  /// Get reference ID of model.
  unsigned long getID() const;
 protected:
  /// Rectangluar support area.
  Rectangle<data_t> support;
  /// Centroid position.
  Point<data_t> centroid;
  /// Reference id.
  unsigned long id;
  /// Compute rectangular SourceModel::support for elliptical sources.
  /// Considers position of SourceModel::centroid.
  void setEllipticalSupport(data_t radius, const complex<data_t>& eps);
};
 
/// Collection of SourceModel entities.
 class SourceModelList : public std::vector<boost::shared_ptr<SourceModel> > {
 public:
  /// Create catalog from all entries of SourceModelList.
  Catalog getCatalog() const;
 };

/// Sersic model class.
/// The model has the form
/// \f[I_S\bigl((x,y)\bigl) = exp\Bigl\lbrace -b_n\Bigl[\Bigl(\frac{r}{R_e}\Bigr)^{1/n} -1\Bigr]\Bigr\rbrace\f]
/// with \f$r=\sqrt{x^2 + y^2}\f$ and \f$b_n\f$ defined by the relation between the complete and the incompete Gamma function,
/// \f[\Gamma(2n)=2\gamma(2n,b_n)\ \Rightarrow\ b_n\approx 1.9992 n - 0.3271.\f]
/// The ensure vanishing flux at large radii, the profile is truncated at 
/// \f$5R_e\f$ and the appropriate value at that position is subtracted from 
/// \f$I_S\f$.\n\n
/// Details can be seen in Graham & Driver, 2005, PASA, 22, 118-127.
class SersicModel : public SourceModel {
 public:
  /// Constructor with Sersic index \p n, effective radius \p Re, \p flux, and
  /// intrinsic ellipticity \p eps.
  SersicModel(data_t n, data_t Re, data_t flux, complex<data_t> eps, const Point<data_t>& centroid, unsigned long id=0);
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const;
  /// Get total flux of model.
  virtual data_t getFlux() const;
  /// Get type of model.
  virtual char getModelType() const;
private:
  data_t n, Re, b,limit,flux,flux_limit,shear_norm,flux_scale;
  complex<data_t> eps;
};

/// Moffat model class.
/// The model has the form
///\f[I_M\bigl((x,y)\bigl) = \bigl(1+\alpha r^2\bigr)^{-\beta}\ \text{with}\ r=\sqrt{x^2 + y^2}\ \text{and}\ \alpha = \frac{2^{1/\beta}-1}{(FWHM/2)^2}.\f]
/// The ensure vanishing flux at large radii, the profile is truncated at
/// \f$5FWHM\f$ and the appropriate value at that position is subtracted from 
/// \f$I_M\f$.
class MoffatModel : public SourceModel {
 public:
  /// Constructor with Moffat index \p beta, width \p FWHM, \p flux, and
  /// intrinsic ellipticity \p eps. 
  MoffatModel(data_t beta, data_t FWHM, data_t flux, complex<data_t> eps, const Point<data_t>& centroid, unsigned long id=0);
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const;
  /// Get total flux of model.
  virtual data_t getFlux() const;
  /// Get type of model.
  virtual char getModelType() const;
 private:
  data_t beta, alpha, limit, flux_limit, flux,shear_norm,flux_scale;
  complex<data_t> eps;
};

/// Model from interpolated pixel data.
/// The class provides several interpolation types for the Image given 
/// at construction time.
class InterpolatedModel : public SourceModel {
public:
  /// Constructor.
  /// The image \p im is assumed to have a Grid starting at <tt>(0,0)</tt>,
  /// by specifying \p reference, its reference point is moved to there without
  /// altering \p im. Similarly, \p flux automatically rescales \p im
  /// to the desired total flux.\n
  /// \p order defines order of interpolation.
  /// - <tt>1</tt>: bi-linear
  /// - <tt>n > 1</tt>: polynomial
  /// - <tt>-3</tt>: bi-cubic
  ///
  /// For more details, see Interpolation.
  InterpolatedModel(const Image<data_t>& im, data_t flux, const Point<data_t>& reference, int order = 1, unsigned long id=0);
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const;
  /// Get total flux of model.
  virtual data_t getFlux() const;
  /// Get type of model.
  virtual char getModelType() const;
private:
  const Image<data_t>& im;
  int order;
  data_t flux,flux_scale;
  Point<data_t> reference;
};

/// Model from ShapeletObject.
/// The class provides convenient sampling from a ShapeletObject.
class ShapeletModel : public SourceModel {
public:
  /// Constructor.
  /// By giving \p centroid, one moves \p sobj to this centroid and adjusts
  /// the support without altering \p sobj itself. Similarly, \p flux 
  /// automatically rescales \p sobj to the desired total flux.\n
  ShapeletModel(const ShapeletObject& sobj, data_t flux, const Point<data_t>& centroid);
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const;
  /// Get total flux of model.
  virtual data_t getFlux() const;
  /// Get type of model.
  virtual char getModelType() const;
private:
  const ShapeletObject& sobj;
  const Point<data_t>& scentroid;
  data_t flux_scale;
};
} // end namespace
#endif
