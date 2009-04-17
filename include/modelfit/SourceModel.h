#ifndef SHAPELENS_SOURCEMODEL_H
#define SHAPELENS_SOURCEMODEL_H

#include <frame/Shapes2D.h>
#include <frame/Object.h>
#include <shapelets/ShapeletObject.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace shapelens {
/// Base class for idealized source models.
/// A galaxy model is a idealized representation of a two-dimensional shape.
/// It's main advantage: It can be sampled at any resolution.\n\n
/// For classes which derive from SourceModel, it is assumed that
/// - they provide a well-defined value at any position \f$(x,y)\f$,
/// - their centroids are set to \f$(0,0)\f$,
/// - their total integrated flux is finite,
/// - they define an finite area of support.

class SourceModel {
public:
  /// Destructor.
  virtual ~SourceModel();
  /// Sample model at \p P.
  virtual data_t getValue(const Point2D<data_t>& P) const = 0;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const = 0;
  /// Get rectangluar support area of the model.
  const Rectangle<data_t>& getSupport() const;
  /// Populate \p obj by sampling from the model.
  /// The model is placed at the centroid of \p obj.\n
  /// The sampled model points are divided by \p normalization.\n
  /// If <tt>add == true</tt>, the model values are added to the entries of \p obj.
  void setObject(Object& obj, data_t normalization, bool add = false) const;
  /// Populate \p obj by sampling from a sheared model.
  void setObjectSheared(Object& obj, complex<data_t> gamma, data_t normalization, bool add = false) const;
 protected:
  /// Rectangluar support area.
  Rectangle<data_t> support;
  /// Centroid position.
  Point2D<data_t> centroid;
};
 
/// Collection of SourceModel entities.
 class SourceModelList : public std::vector<boost::shared_ptr<SourceModel> > {
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
  /// Constructor with Sersic index \p n, effective radius \p Re and
  /// intrinsic ellipticity \p eps.
  SersicModel(data_t n, data_t Re, complex<data_t> eps, const Point2D<data_t>& centroid);
  /// Sample model at \p P.
  virtual data_t getValue(const Point2D<data_t>& P) const;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const;
private:
  data_t n, Re, b,limit,flux,flux_limit;
  complex<data_t> eps;
};

/// Sersic model class.
/// The model has the form
///\f[I_M\bigl((x,y)\bigl) = \bigl(1+\alpha r^2\bigr)^{-\beta}\ \text{with}\ r=\sqrt{x^2 + y^2}\ \text{and}\ \alpha = \frac{2^{1/\beta}-1}{(FWHM/2)^2}.\f]
/// The ensure vanishing flux at large radii, the profile is truncated at
/// \f$5FWHM\f$ and the appropriate value at that position is subtracted from 
/// \f$I_M\f$.
class MoffatModel : public SourceModel {
 public:
  /// Constructor with Moffat index \p beta, width \p FWHM and
  /// intrinsic ellipticity \p eps. 
  MoffatModel(data_t beta, data_t FWHM, complex<data_t> eps, const Point2D<data_t>& centroid);
  /// Sample model at \p P.
  virtual data_t getValue(const Point2D<data_t>& P) const;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const;
 private:
  data_t beta, alpha, limit, flux_limit, flux;
  complex<data_t> eps;
};

/// Model from interpolated pixel data.
/// The class provides several interpolation types for the Object given 
/// at construction time.
class InterpolatedModel : public SourceModel {
public:
  /// Constructor.
  /// \p order defines order of interpolation.
  /// - <tt>1</tt>: bi-linear
  /// - <tt>n > 1</tt>: polynomial
  /// - <tt>-3</tt>: bi-cubic
  ///
  /// For more details, see Interpolation.
  InterpolatedModel(const Object& obj, int order = 1);
  /// Sample model at \p P.
  virtual data_t getValue(const Point2D<data_t>& P) const;
  /// Get Object::flux.
  virtual data_t getFlux() const;
private:
  const Object& obj;
  int order;
};

/// Model from ShapeletObject.
/// The class provides convenient sampling from a ShapeletObject.
class ShapeletModel : public SourceModel {
public:
  /// Constructor.
  ShapeletModel(const ShapeletObject& sobj);
  /// Sample model at \p P.
  virtual data_t getValue(const Point2D<data_t>& P) const;
  /// Return ShapeletObject::getShapeletFlux().
  virtual data_t getFlux() const;
private:
  const ShapeletObject& sobj;
};
} // end namespace
#endif
