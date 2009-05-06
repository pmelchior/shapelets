#ifndef SHAPELENS_OPTIMALDECOMPOSITE2D_H
#define SHAPELENS_OPTIMALDECOMPOSITE2D_H

#include <gsl/gsl_vector.h>
#include <map>
#include <bitset>
#include <numla/NumMatrix.h>
#include <numla/NumVector.h>
#include "../Typedef.h"
#include "../frame/Grid.h"
#include "../frame/Object.h"
#include "Decomposite2D.h"

namespace shapelens {

/// Class for optimal shapelet decomposition.
/// Provides minimization of the decomposition's \f$\chi^2\f$.
/// This class delivers best fit \f$\beta\f$ and \f$n_{max}\f$, assuming a known position
/// of the centroid \f$x_c\f$.\n\n
/// The procedure:
/// - Create Decomposite2D object with \f$n_{max}=2\f$ as starting point.
/// - Search for \f$\beta\f$ that minimizes \f$\chi^2\f$.
/// - If above minimizations does not converge, increase \f$n_{max}\f$ by 2.\n
///   If this is still not sufficient, there are probably image distortions
///   close to the object and the decomposition will therefore be aborted. 
/// - Increase \f$n_{max}\f$ in steps of 2 until 
///  - \f$\chi^2 < 1 + \sigma(\chi^2)\f$ or 
///  - \f$\chi^2\f$ flattens out (if \p ShapeLensConfig::ALLOW_FLATTENING is set) or 
///  - correlation function of the residuals becomes lower than the one of the noise
///    (if <tt>ShapeLensConfig::NOISEMODEL == COVARIANCE</tt>)
/// - If \f$n_{max} = 6\f$ or \f$n_{max} \ mod\ 12 = 0\f$, do an intermediate search 
///   for \f$\beta\f$ at current value of \f$n_{max}\f$
/// - Search for \f$\beta\f$ that minimizes \f$\chi^2\f$ at the best-fit \f$n_{max}\f$.
/// - Reset to a somewhat lower \f$n_{max}\f$ given the current best-fit value of \f$\beta\f$
///   and search for a potentially lower \f$n_{max}\f$ (now in steps of 1).
/// - If \f$n_{max}\f$ has decreased during last step, seach again for  
///   \f$\beta\f$. If \f$n_{max}\f$ has increased now, ignore it, since
///   we have had already better values.
///
/// Apart from that, the optimization is constrained by these configuration parameters
/// - \p ShapeLensConfig::ALLOW_FLATTENING
/// - \p ShapeLensConfig::BETA_LOW 
/// - \p ShapeLensConfig::BETA_HIGH
/// - \p ShapeLensConfig::DELTA_BETA
/// - \p ShapeLensConfig::NMAX_LOW 
/// - \p ShapeLensConfig::NMAX_HIGH
///
/// and additionally bound by \f$2 \theta_{min}> 1\f$ and \f$\theta_{max}<\frac{D}{2}\f$ with \f$D\f$
/// being the minimum sidelength of the image frame.
///
/// See Paper III, sect. 3.4 and 7.5, and Paper IV, sect. 4.3, for details.

class OptimalDecomposite2D : private Decomposite2D {
 public:
  /// Constructor for decomposing an Object.
  OptimalDecomposite2D(const Object& obj, Composite2D& model);
  /// Return best fit  \f$\beta\f$
  data_t getOptimalBeta();
  /// Return best fit \f$\chi^2\f$.
  data_t getOptimalChiSquare();
  /// Return the decomposition flags.
  /// if \p decompositionFlags(i) is set, it indicates problems during the optimization:
  /// - <tt>i = 0</tt>: \f$\chi^2 > 1\f$
  /// - <tt>i = 1</tt>: optimization stopped due to either 
  ///   <tt>ShapeLensConfig::ALLOW_FLATTENING = 1</tt> or the correlation function of
  ///   residuals becoming lower than the one of the noise (only if 
  ///   <tt>ShapeLensConfig::NOISEMODEL = COVARIANCE</tt>).
  /// - <tt>i = 2</tt>: \f$n_{max}\f$ limited by <tt>ShapeLensConfig::NMAX_HIGH</tt>
  /// - <tt>i = 3</tt>: \f$n_{max}\f$ limited due to \f$2\ \theta_{min}<1\f$
  /// - <tt>i = 4</tt>: \f$n_{max}\f$ limited due to \f$n_{pixels} \leq n_{coeffs}\f$
  /// - <tt>i = 5</tt>: \f$\chi^2\f$ was not optimized because constraints from 
  ///   \p ShapeLensConfig::BETA_LOW or \p ShapeLensConfig::BETA_HIGH were too severe.
  /// - <tt>i = 6</tt>: \f$\chi^2\f$ became worse during optimization and did not improve
  ///    anymore.
  /// - <tt>i = 7</tt>: \f$\chi^2(\beta) \bigl|_{n_{max}=2}\f$ does not have a useful minimum.
  const std::bitset<8>& getDecompositionFlags();
  /// Get the decomposition History.
  std::string getHistory();

private:
  int npixels, optimalNMax;
  data_t optimalBeta,optimalChiSquare,bestChiSquare,image_dimension;
  bool nmaxTrouble, noise_correlated;
  char comp_corr;
  std::bitset<8> flags;
  History history;
  std::string comp_corr_string;
  std::map<int, data_t> bestBeta, bestChi2;
  int findOptimalBeta(unsigned char step);
  void findOptimalNMax(unsigned char step), optimize(), checkBeta();
  void getBetaTrafoMatrix(NumMatrix<data_t>& betaTrafo, data_t beta1, data_t beta2);
  void checkCorrelationFunctionFromResiduals();
  data_t getBetaLimit(bool upper);

};
} // end namespace
#endif
