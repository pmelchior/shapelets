#ifndef SHAPELENS_COMPOSITE2D_H
#define SHAPELENS_COMPOSITE2D_H

#include <numla/NumMatrix.h>
#include <numla/NumVector.h>
#include "../Typedef.h"
#include "../frame/Image.h"
#include "../frame/Moments.h"
#include "CoefficientVector.h"

namespace shapelens {

/// Shapelet composition class.
/// Provides methods associated with a 2D shapelet model, 
/// \f[f(x_1,x_2) = \sum_{n_1,n_2}^{n_1 + n_2 = n_{max}} f_{n_1,n_2}\ B_{n_1,n_2}(x_1,x_2;\beta)\f]
///
/// Among them are methods to
/// - evaluate a model an a given Grid or at a single point,
/// - integrate the model over the whole 2D plane or within a finite region,
/// - measure lowest order moments of the model (flux, centroid, quadrupole).
///
/// At each evaluation, a covariance matrix of the result can be obtained:
/// \code
/// Composite model ...;
/// data_t integral = model.integrate();
/// NumMatrix<data_t> cov_quad;
/// NumMatrix<data_t> Q = model.getShapelet2ndMoments(&cov_quad);
/// \endcode
/// The example calculates the integral \f$\int dx\ f(x)\f$ without and the quadrupole 
/// moment \f$Q\f$ with the respective covariance matrix.

class Composite {
 public:
  /// Default constructor.
  /// Don't use it for explicit composition.
  Composite();
  /// Argumented constructor.
  Composite(const CoefficientVector<data_t>& Coeffs, data_t beta, Point<data_t>& xcentroid);
  /// Get cartesian shapelet coefficients.
  const CoefficientVector<data_t>& getCoeffs() const;
  /// Set cartesian shapelet coefficients.
  void setCoeffs(const CoefficientVector<data_t>& newCoeffs);
  /// Get the coviariance matrix of the shapelet coefficients.
  const NumMatrix<data_t>& getCovarianceMatrix() const;
  /// Set coefficient covariance matrix.
  void setCovarianceMatrix(const NumMatrix<data_t>& cov);
  /// Get \f$n_{max}\f$, the maximum order of the shapelet model.
  unsigned int getNMax() const;
  /// Set \f$n_{max}\f$.
  /// The underlying CoefficientVector is changed accordingly.
  void setNMax(unsigned int nmax);
  /// Get \f$\beta\f$ from basis function.
  data_t getBeta() const;
  /// Set new \f$\beta\f$ for basis functions.
  void setBeta(data_t beta);
  /// Get centroid position \f$x_c\f$
  const Point<data_t>& getCentroid() const;
  /// Set centroid position \f$x_c\f$
  void setCentroid(const Point<data_t>& inxcentroid);
  /// Get current Grid.
  const Grid& getGrid() const;
  /// Set new Grid.
  /// This changes the default Grid and therefore also the behavior of evalGrid(),
  /// which depends on the stepsize between the grid points.
  void setGrid(const Grid& ingrid);
  /// Get the shapelet model.
  /// This evaluates \f$f\f$ on the whole grid.
  const Image<data_t>& getModel();
  /// Evaluate \f$f(x)\f$.
  /// When given, \p cov will be the (1,1) covariance matrix (= error squared) of \f$f(x)\f$.
  data_t eval(const Point<data_t>& x, NumMatrix<data_t>* cov = NULL) const;
  /// Integrate \f$f(x)\f$.
  /// When given, \p cov will be the (1,1) covariance matrix (= error squared) of 
  /// \f$\int\ dx\ f(x)\f$.
  data_t integrate(NumMatrix<data_t>* cov = NULL) const;
  /// Integrate \f$f(x)\f$ in the area bounded by \p P1 and \p P2.
  /// When given, \p cov will be the (1,1) covariance matrix (= error squared) of 
  /// the integral.
  data_t integrate(const Point<data_t>& P1,const Point<data_t>& P2,  NumMatrix<data_t>* cov = NULL) const;
  /// Calculate the object flux \f$F\f$ from the coefficients.
  /// cf. Paper I, eq. 26.\n
  /// When given, \p cov will be the (1,1) covariance matrix (= error squared) of \f$F\f$
  data_t getShapeletFlux(NumMatrix<data_t>* cov = NULL) const;
  /// Get the object centroid \f$\vec{x}_c\f$ from the coefficients.
  /// cf. Paper I, eq. 27\n
  /// When given, \p cov will be the (2,2) covariance matrix of \f$\vec{x}_c\f$.
  Point<data_t> getShapeletCentroid(NumMatrix<data_t>* cov = NULL) const;
  /// Get 2nd brightness moments \f$Q_{ij}\f$.
  /// When given, \p cov will be the (3,3) covariance matrix of \f$Q_{ij}\f$, ordered as
  /// \f$Q_{11},\ Q_{12},\ Q_{22}\f$.
  Quadrupole getShapelet2ndMoments(NumMatrix<data_t>* cov = NULL) const;
  /// Get the object RMS radius \f$r\f$ from the coefficients.
  /// cf. Paper I, eq. 28\n
  /// When given, \p cov will be the (1,1) covariance matrix (= error squared) of \f$r\f$.
  data_t getShapeletRMSRadius(NumMatrix<data_t>* cov = NULL) const;

  friend class SIFFile;
  friend class Decomposite;

 protected:
  /// The scale size
  data_t beta;
  /// The shapelet coefficients.
  CoefficientVector<data_t> coeffs;
  /// The shapelet design matrix.
  NumMatrix<data_t> M, MInt;
  /// The coefficient covariance matrix.
  NumMatrix<data_t> cov;
  /// The shapelet model (which contains the Grid).
  Image<data_t> model;
  /// The centroid position.
  Point<data_t> xcentroid;
  /// Whether M must be updated.
  bool changeM;
  /// Wheter model must be updated.
  bool changeModel;

 private:
  void evalGrid();
  void makeShapeletMatrix();
  void updateOrders();
};
} // end namespace
#endif
