#ifndef COMPOSITE2D_H
#define COMPOSITE2D_H

/// 2D Composition class.
/// Provides values of a 2D composite shapelet function 
/// \f$f(x0,x1) = \sum_{n0,n1}^{n_{max}} f_{n0,n1} \cdot B_{n0,n1}(x0,x1;\beta)\f$.\n

#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <shapelets/Shapelets2D.h>
#include <shapelets/CoefficientVector.h>

class Composite2D : private Shapelets2D {
 public:
  /// Default constructor.
  /// Don't use it for explicit composition.
  Composite2D();
  /// Argumented constructor.
  /// Different shapelet orders are allowed by giving matrix of appropriate dimensions.
  Composite2D(data_t beta, Point2D& xcentroid, const NumMatrix<data_t>& Coeffs);
  /// Copy constructor.
  Composite2D(const Composite2D &source);
  /// Copy operator
  Composite2D & operator = (const Composite2D &source);
  /// Get cartesian shapelet coefficients.
  const CoefficientVector<data_t>& getCoeffs() const;
  /// Set cartesian shapelet coefficients.
  void setCoeffs(const CoefficientVector<data_t>& newCoeffs);
  /// Get \f$n_{max}\f$, the maximum order of the shapelet model.
  unsigned int getNMax() const;
  /// Get \f$\beta\f$ from basis function.
  data_t getBeta() const;
  /// Set new \f$\beta\f$ for basis functions.
  void setBeta(data_t beta);
  /// Get centroid position \f$x_c\f$
  const Point2D& getCentroid() const;
  /// Set centroid position \f$x_c\f$
  void setCentroid(const Point2D& inxcentroid);
  /// Get current Grid.
  const Grid& getGrid() const;
  /// Set new Grid.
  /// This changes the default Grid and therefore also the behavior of evalGrid(),
  /// which depends on the stepsize between the grid points.
  void setGrid(const Grid& ingrid);
  /// Evaluate \f$f(x)\f$.
  data_t eval(const Point2D& x);
  /// Get the shapelet model.
  /// This evaluates \f$f(x)\f$ on the whole grid.
  const NumVector<data_t>& getModel();
  /// Integrate \f$f(x)\f$.
  data_t integrate();
  /// Integrate \f$f(x)\f$ in the range (x0min,x1min) .. (x0max,x1max).
  data_t integrate(data_t x0min, data_t x0max, data_t x1min,data_t x1max);
  /// Calculate the object flux from the coefficients.
  /// see Paper I, eq. 26
  data_t getShapeletFlux() const;
  /// Get the object centroid from the coefficients.
  /// see Paper I, eq. 27
  Point2D getShapeletCentroid() const;
  /// Get 2nd brightness moments \f$Q_{ij}\f$.
  NumMatrix<data_t> getShapelet2ndMoments() const;
  /// Get the object RMS radius from the coefficients.
  /// see Paper I, eq. 28
  data_t getShapeletRMSRadius() const;

  friend class SIFFile;

 protected:
  Grid grid;
  CoefficientVector<data_t> coeffs;
  NumMatrix<data_t> M;
  NumVector<data_t> model;
  Point2D xcentroid;
  bool change;

 private:
  void evalGrid();
  void makeShapeletMatrix(NumMatrix<data_t>& M);
  NumVector<data_t>& accessModel();
  void updateOrders();
};

#endif
