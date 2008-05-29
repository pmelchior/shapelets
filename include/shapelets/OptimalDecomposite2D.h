#ifndef OPTIMAL_DECOMPOSITE2D_H
#define OPTIMAL_DECOMPOSITE2D_H

#include <gsl/gsl_vector.h>
#include <map>
#include <bitset>
#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <History.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <frame/Object.h>
#include <shapelets/Decomposite2D.h>

/// Class for optimal 2D decomposition.
/// Provides minimization of decomposition's \f$\chi^2\f$.\n
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
///  - \f$\chi^2\f$ flattens out (if ShapeLensConfig::ALLOW_FLATTENING is set) or 
///  - correlation function of the residuals becomes lower than the one of the noise
///    (if <tt>ShapeLensConfig::NOISEMODEL = COVARIANCE</tt>)
/// - If \f$n_{max}\ mod\ 6 = 0\f$, do an intermediate search for \f$\beta\f$ at current
///   value of \f$n_{max}\f$
/// - Search for \f$\beta\f$ that minimizes \f$\chi^2\f$ at the best-fit \f$n_{max}\f$.
/// - Reset to a somewhat lower \f$n_{max}\f$ given the current best-fit value of \f$\beta\f$
///   and search for a potentially lower \f$n_{max}\f$ (now in steps of 1).
/// - If \f$n_{max}\f$ has decreased during last step, seach again for  
///   \f$\beta\f$. If \f$n_{max}\f$ has increased now, ignore it, since
///   we have had already better values.
///
/// See Paper III, sec. 3.4 and 7.5 for details.

class OptimalDecomposite2D : private Decomposite2D {
 public:
  /// Constructor for decomposing an Object.
  /// The optimization is constrained by the parameters
  /// - ShapeLensConfig::ALLOW_FLATTENING
  /// - ShapeLensConfig::BETA_LOW 
  /// - ShapeLensConfig::BETA_HIGH
  /// - ShapeLensConfig::DELTA_BETA
  /// - ShapeLensConfig::NMAX_LOW 
  /// - ShapeLensConfig::NMAX_HIGH
  OptimalDecomposite2D(const Object& obj);
  /// Employs a regularization to lower the negative flux.
  /// The minimization ends, whenever \f$R=- F^-/F^+ < \text{wanted}R\f$, where \f$F^{\pm}\f$
  /// is the positive or negative flux of the shapelet reconstruction. 
  /// The return value is the actual \f$R\f$.\n
  /// If \p wantedR is chosen too low (e.g. \f$10^{-6}\f$, in some cases already 
  /// \f$10^{-3}\f$), the minimization procedure
  /// must increase \f$n_{max}\f$ and can therefore take a considerable amount of time.
  data_t regularize(data_t wantedR);
  /// Get best fit shapelet coefficients.
  const CoefficientVector<data_t>& getCoeffs();
   /// Get covariance matrix of the best fit shapelet coefficients.
  NumMatrix<data_t> getCovarianceMatrix();
  /// Get best fit residuals of the decomposition.
  const NumVector<data_t>& getResiduals();
  /// Get best fit shapelet model.
  const NumVector<data_t>&  getModel();
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
  /// - <tt>i = 5</tt>: \f$\chi^2\f$ was not optimized because constraints from ShapeLensConfig::BETA_LOW
  ///   or ShapeLensConfig::BETA_HIGH were too severe.
  /// - <tt>i = 6</tt>: \f$\chi^2\f$ became worse during optimization and did not improve
  ///    anymore.
  /// - <tt>i = 7</tt>: \f$\chi^2(\beta) \bigl|_{n_{max}=2}\f$ does not have a useful minimum.
  const std::bitset<8>& getDecompositionFlags();
  /// Get the decomposition History.
  std::string getHistory();

private:
  const Object& obj;
  int npixels, optimalNMax;
  data_t beta,optimalBeta,optimalChiSquare,bestChiSquare,image_dimension;
  bool optimized, nmaxTrouble, noise_correlated;
  char comp_corr;
  std::bitset<8> flags;
  History history;
  std::string comp_corr_string;
  std::map<int, data_t> bestBeta, bestChi2;
  struct regResults {
    int nmax;
    data_t f;
    data_t lambda;
    data_t beta;
    data_t chi2;
    data_t R;
    NumVector<data_t> coeffs;
  };
  int findOptimalBeta(unsigned char step);
  void findOptimalNMax(unsigned char step), optimize(), checkBeta();
  void getBetaTrafoMatrix(NumMatrix<data_t>& betaTrafo, data_t beta1, data_t beta2);
  void appendRegResults(std::vector<regResults>& results, int nmax, data_t f, data_t lambda, data_t beta, data_t chi2, data_t R, const NumVector<data_t>& coeffs);
  int findNMaxofBestF(std::vector<regResults>& results);
  int findSuboptimalResultIndex(std::vector<regResults>& result);
  void checkCorrelationFunctionFromResiduals();
  data_t getBetaLimit(bool upper);

};

#endif
