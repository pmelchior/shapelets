#ifndef COMPOSITE1D_H
#define COMPOSITE1D_H

#include <NumVector.h>
#include <Typedef.h>
#include <frame/Grid.h>
#include <shapelets/Shapelets1D.h>
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
  Composite1D(data_t beta, data_t xcentroid, const NumVector<data_t>& coeffs);
  /// Get maximal order for composition.
  int getOrder();
  /// Set maximal order for composition.
  /// This introduces a cutoff for the shapelet series.
  /// It's usefull to see the effect of truncation.
  void setOrderLimit(int orderlimit);
  /// Get \f$\beta\f$ from basis function.
  data_t getBeta();
  /// Set new \f$\beta\f$ for basis functions.
  void setBeta(data_t beta);
   /// Get smallest reproducible object size \f$\theta_{min}\f$.
  data_t getThetaMin();
  /// Get biggest reproducible object size \f$\theta_{max}\f$.
  data_t getThetaMax();
  /// Set new Grid for evaluation.
  /// Set the shapelet coeffs.
  /// Enlarge basis set when needed.
  void setCoeffs(const NumVector<data_t>& newCoeffs);
    /// Set new Grid for evaluation.
  void setGrid(const Grid& grid);
  /// Get active Grid.
  /// This is especially useful when retrieving the default grid.
  const Grid& getGrid();
  /// Evaluate \f$f(x)\f$.
  data_t eval(data_t x);
  /// Evaluate \f$f(x)\f$ on a 1D grid.
  void evalGrid(NumVector<data_t>& values);
  /// Integrate \f$f(x)\f$
  data_t integrate();
  /// Integrate  \f$f(x)\f$ within the range xmin..xmax.
  data_t integrate(data_t xmin, data_t xmax);
  
 private:
  NumVector<data_t> shapeletCoeffs;
  int order,orderlimit;
  data_t xcentroid;
  Grid grid;
};

#endif
