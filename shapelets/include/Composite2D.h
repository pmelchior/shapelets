#ifndef COMPOSITE2D_H
#define COMPOSITE2D_H

/// 2D Composition class.
/// Provides values of a 2D composite shapelet function 
/// \f$f(x0,x1) = \sum f_{n0,n1} \cdot B_{n0,n1}(x0,x1)\f$.\n
/// A default Grid will be defined on construction, that will be appropriately sized
/// to plot the whole image and all included details.

#include <Shapelets2D.h>
#include <NumMatrix.h>
#include <NumVector.h>
#include <Point2D.h>
#include <Grid.h>

class Composite2D : private Shapelets2D {
 public:
  /// Default constructor.
  /// Don't use it for explicit composition.
  Composite2D();
  /// Argumented constructor.
  /// Different shapelet orders are allowed by giving matrix of appropriate dimensions.
  Composite2D(double beta, Point2D& xcentroid, const NumMatrix<double>& startCoeffs);
  /// Copy operator
  Composite2D & operator = (const Composite2D &source);
  /// Get maximal order for composition in direction (0/1).
  int getOrder(bool direction) const;
  /// Get \f$n_{max}\f$, the maximum order of the shapelet model.
  /// This assumes that the orders in direction 0 and 1 are identical.
  int getNMax() const;
  /// Set the maximum composition order without affecting the shapelet coefficients.
  /// This lowers the composition order to see effect of truncation. */
  void setOrderLimit(bool direction, int orderlimit);
  /// Set the shapelet coefficients to new values.
  /// Enlarge basis set when needed
  void setCoeffs(const NumMatrix<double>& newCoeffs);
  /// Get \f$\beta\f$ from basis function.
  double getBeta() const;
  /// Access \f$\beta\f$ from basis function.
  double& accessBeta();
  /// Set new \f$\beta\f$ for basis functions.
  void setBeta(double beta);
  /// Get centroid position \f$x_c\f$
  const Point2D& getCentroid() const;
  /// Access centroid position \f$x_c\f$
  Point2D& accessCentroid();
  /// Set centroid position \f$x_c\f$
  void setCentroid(const Point2D& inxcentroid);
  /// Get active Grid.
  /// This is especially useful for retrieving the default grid.
  const Grid& getGrid();
  /// Set new Grid.
  /// This changes the default Grid and therefore also the behavior of evalGrid(),
  /// which depends on the stepsize between the grid points.
  void setGrid(const Grid& ingrid);
  /// Evaluate \f$f(x)\f$.
  double eval(const Point2D& x);
  /// Get the shapelet model.
  /// This evaluates \f$f(x)\f$ on the whole grid.
  const NumVector<double>& getModel();
  /// Access the shapelet model directly.
  NumVector<double>& accessModel();
  /// Integrate \f$f(x)\f$.
  double integrate();
  /// Integrate \f$f(x)\f$ in the range (x0min,x1min) .. (x0max,x1max).
  double integrate(double x0min, double x0max, double x1min,double x1max);
  /// Calculate the object flux from the coefficients.
  /// see Paper I, eq. 26
  double getShapeletFlux() const;
  /// Calculate the object centroid from the coefficients.
  /// see Paper I, eq. 27
  void getShapeletCentroid(Point2D& xc) const;
  /// Calculate 2nd brightness moments \f$Q_{ij}\f$.
  void getShapelet2ndMoments(NumMatrix<double>& Q) const;
  /// Calculate the object RMS radius from the coefficients.
  /// see Paper I, eq. 28
  double getShapeletsRMSRadius() const;
 private:
  Grid grid;
  NumMatrix<double> shapeletCoeffs, M;
  NumVector<double> model;
  Point2D xcentroid;
  int order0, order1, orderlimit0, orderlimit1;
  double beta, stepsize0, stepsize1;
  bool change, changeGrid, lockGrid;
  // defines dimensions of the grid depending on the shapelet orders
  void defineGrid();
  double evalGridPoint(const Point2D& x);
  void evalGrid();
  void makeShapeletMatrix(NumMatrix<double>& M, NumMatrix<int>& nVector);
};

#endif
