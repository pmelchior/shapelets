#ifndef GALAXYMODEL_H
#define GALAXYMODEL_H

#include <frame/Object.h>
#include <shapelets/ShapeletObject.h>

/// Base class for idealized galaxy models.
/// A galaxy model is a idealized representation of a two-dimensional shape.
/// It's bigget advantage: It can be sampled at any resolution.\n\n
/// For classes which derive from GalaxyModel, it is assumed that
/// - they provide a well-defined value at any position \f$(x,y)\f$
/// - their centroids are set to \f$(0,0)\f$
/// - their total integrated flux is finite.
class GalaxyModel {
public:
  /// Destructor.
  virtual ~GalaxyModel();
  /// Sample model at \f$(x,y)\f$.
  /// When <tt>normed == true</tt>, all value are normalized to unit flux.
  virtual data_t getValue(data_t x, data_t y,bool normed = false) const = 0;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const = 0;
  /// Populate \p obj by sampling from the model.
  /// The model will be placed at the centroid of \p obj.\n
  /// If <tt>normed == true</tt>, the resulting object has unit flux.\n
  /// If <tt>add == true</tt>, the model values are added to the entries of \p obj.
  void setObject(Object& obj, bool normed = false, bool add = false) const;
  /// Populate \p obj by sampling from a sheared model.
  void setObjectSheared(Object& obj, complex<data_t> gamma, bool normed = false, bool add = false) const;
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
class SersicModel : public GalaxyModel {
 public:
  /// Constructor with Sersic index \p n and effective radius \p Re.
  SersicModel(data_t n, data_t Re);
  /// Sample model at \f$(x,y)\f$.
  /// When <tt>normed == true</tt>, all value are normalized to unit flux.
  virtual data_t getValue(data_t x, data_t y, bool normed = false) const;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const;
private:
  data_t n, Re, b,limit,flux,flux_limit;
};

/// Sersic model class.
/// The model has the form
///\f[I_M\bigl((x,y)\bigl) = \bigl(1+\alpha r^2\bigr)^{-\beta}\ \text{with}\ r=\sqrt{x^2 + y^2}\ \text{and}\ \alpha = \frac{2^{1/\beta}-1}{(FWHM/2)^2}.\f]
/// The ensure vanishing flux at large radii, the profile is truncated at
/// \f$5FWHM\f$ and the appropriate value at that position is subtracted from 
/// \f$I_M\f$.
class MoffatModel : public GalaxyModel {
 public:
  /// Constructor with Moffat index \p beta and width \p FWHM.
  MoffatModel(data_t beta, data_t FWHM);
  /// Sample model at \f$(x,y)\f$.
  /// When <tt>normed == true</tt>, all value are normalized to unit flux.
  virtual data_t getValue(data_t x, data_t y, bool normed = false) const;
  /// Get total integrated flux of model.
  virtual data_t getFlux() const;
 private:
  data_t beta, alpha, limit, flux_limit, flux;
};

/// Model from interpolated pixel data.
/// The class provides a bilinear interpolation of the Object given 
/// at construction time.
class InterpolatedModel : public GalaxyModel {
public:
  /// Constructor.
  InterpolatedModel(Object& obj);
  /// Sample model at \f$(x,y)\f$.
  /// When <tt>normed == true</tt>, all value are normalized to unit 
  /// Object::flux.
  virtual data_t getValue(data_t x, data_t y,bool normed = false) const;
  /// Get Object::flux.
  virtual data_t getFlux() const;
private:
  const Object& obj;
};

/// Model from ShapeletObject.
/// The class provides convenient sampling from a ShapeletObject.
class ShapeletModel : public GalaxyModel {
public:
  /// Constructor.
  ShapeletModel(ShapeletObject& sobj);
  /// Sample model at \f$(x,y)\f$.
  /// When <tt>normed == true</tt>, all value are normalized to unit flux.
  virtual data_t getValue(data_t x, data_t y,bool normed = false) const;
  /// Return ShapeletObject::getShapeletFlux().
  virtual data_t getFlux() const;
private:
  ShapeletObject& sobj;
};


#endif
