#ifndef COMPOSITE1D_H
#define COMPOSITE1D_H

#include <Grid.h>
#include <NumVector.h>
#include <Shapelets1D.h>
#include <gsl/gsl_math.h>

/// 1D Composition class.
///  Provides values of a 1D composite shapelet function 
/// \f$f(x) = \sum f_n \cdot B_n(x)\f$.
/// A default Grid will be defined on construction, that will be appropriatly sized
/// to plot the whole image in all included details.

class Composite1D : private Shapelets1D {
 public:
  /// Argumented constructor.
  /// Shapelet coefficients determine the maximum shapelet order.
  Composite1D(double beta, double xcentroid, const NumVector<double>& coeffs);
  /// Get maximal order for composition.
  int getOrder();
  /// Set maximal order for composition.
  /// This introduces a cutoff for the shapelet series.
  /// It's usefull to see the effect of truncation.
  void setOrderLimit(int orderlimit);
  /// Get \f$\beta\f$ from basis function.
  double getBeta();
  /// Set new \f$\beta\f$ for basis functions.
  void setBeta(double beta);
   /// Get smallest reproducible object size \f$\theta_{min}\f$.
  double getThetaMin();
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  double getThetaMax();
  /// Set new Grid for evaluation.
  /// Set the shapelet coeffs.
  /// Enlarge basis set when needed.
  void setCoeffs(const NumVector<double>& newCoeffs);
    /// Set new Grid for evaluation.
  void setGrid(const Grid& grid);
  /// Get active Grid.
  /// This is especially useful when retrieving the default grid.
  const Grid& getGrid();
  /// Evaluate \f$f(x)\f$.
  double eval(double x);
  /// Evaluate \f$f(x)\f$ on a 1D grid.
  void evalGrid(NumVector<double>& values);
  /// Integrate \f$f(x)\f$
  double integrate();
  /// Integrate  \f$f(x)\f$ within the range xmin..xmax.
  double integrate(double xmin, double xmax);
  
 private:
  NumVector<double> shapeletCoeffs;
  int order,orderlimit;
  double xcentroid;
  Grid grid;
};

#endif
