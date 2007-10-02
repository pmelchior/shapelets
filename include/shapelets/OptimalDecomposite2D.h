#ifndef OPTIMAL_DECOMPOSITE2D_H
#define OPTIMAL_DECOMPOSITE2D_H

#include <gsl/gsl_vector.h>
#include <map>
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
/// This class delivers best fit \f$\beta\f$, 
/// centroid position \f$x_c\f$ and \f$n_{max}\f$.\n\n
/// The procedure:
/// - Create Decomposite2D object of order (2,2) as starting point.
/// - Search for \f$\beta\f$ and \f$x_c\f$ that minimize \f$\chi^2\f$.
/// - If above minimizations doesn't converge with, increase order by 2
///   in each direction.\n
///   If this is still not sufficient, there are probably image distortions
///   close to the object and the decomposotion will therefore be aborted. 
///   Try using an appropriate filter before decomposition.
/// - Increase \f$n_{max}\f$ in steps of 2 until the decomposition has
///   \f$\chi^2 < 1 + \sigma(\chi^2)\f$ or flattens out (for cases with 
///   \f$\chi^2 > 2\f$ only).
/// - Search again for the \f$\beta\f$ and \f$x_c\f$ that minimizes \f$\chi^2\f$.
/// - Reset the order to (2,2) and search for a new \f$n_{max}\f$ (now in steps of 1 
///   for possibly lower \f$n_{max}\f$ then before) with the presumably
///   better values for \f$\beta\f$ and \f$x_c\f$.
/// - If \f$n_{max}\f$ has decreased during last step, seach a last time for  
///   \f$\beta\f$ and \f$x_c\f$. If \f$n_{max}\f$ has increased now, ignore it, since
///   we have had already better values.
///
/// See Paper III, sec. 3.4 and 7.5 for details.

class OptimalDecomposite2D : private Decomposite2D {
 public:
  /// Constructor for decomposing an Object.
  /// nmaxLimit is the Limit for the maximum order of the decompososition.
  /// Estimators for \f$\beta\f$ and image centroid  \f$x_c\f$
  /// are obtained from Object, also the noisemodel.
  OptimalDecomposite2D(const Object& obj, int nmaxLow, int nmaxHigh, data_t betaLow, data_t betaHigh);
  /// Employs a regularization to lower the negative flux.
  /// The minimization ends, whenever \f$R=- F^-/F^+ < wantedR\f$, where \f$F^{\pm}\f$
  /// is the positive or negative flux of the shapelet reconstruction. 
  /// The return value is the actual \f$R\f$.\n
  /// If wantedR is chosen too low (e.g. \f$10^{-6}\f$, in some cases already 
  /// \f$10^{-5}\f$), the minimization procedure
  /// must increase \f$n_{max}\f$ and can therefore take a considerable amount of time.
  data_t regularize(data_t wantedR);
  /// Deliver best fit shapelet coefficients.
  void getShapeletCoeffs(NumMatrix<data_t>& coeffs);
   /// Deliver best fit shapelet coefficient errors.
  void getShapeletErrors(NumMatrix<data_t>& errors);
  /// Get best fit residuals of the decomposition.
  const NumVector<data_t>& getResiduals();
  /// Get best fit shapelet model.
  const NumVector<data_t>&  getModel();
  /// Return best fit  \f$\beta\f$
  data_t getOptimalBeta();
  /// Return best fit \f$\chi^2\f$.
  data_t getOptimalChiSquare();
  /// Return the decomposition flag.
  /// It indicates problems during the optimization.
  /// - 0: OK
  /// - 1: \f$n_{max}\f$ limited by user defined range
  /// - 2: \f$n_{max}\f$ limited due to \f$2\theta_{min}<1\f$
  /// - 3: \f$n_{max}\f$ limited due to \f$n_{pixels} \leq n_{coeffs}\f$
  /// - 4: decomposition not optimal, object probably ill-mathced by shapelets
  /// - 5: optimal \f$\beta\f$ not found.
  char getDecompositionFlag();
  /// Get the decomposition History.
  const History& getHistory();

private:
  const Object& obj;
  int npixels, optimalNMax, nmaxLow, nmaxHigh;
  data_t beta,optimalBeta,optimalChiSquare,bestChiSquare,image_dimension, betaHigh, betaLow;
  bool optimized, nmaxTrouble, noise_correlated;
  char flag, comp_corr;
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
  void getCoeffErrorFromBeta(const NumVector<data_t>& coeffVector, NumVector<data_t>& errorVector);
  void getBetaTrafoMatrix(NumMatrix<data_t>& betaTrafo, data_t beta1, data_t beta2);
  void appendRegResults(std::vector<regResults>& results, int nmax, data_t f, data_t lambda, data_t beta, data_t chi2, data_t R, const NumVector<data_t>& coeffs);
  int findNMaxofBestF(std::vector<regResults>& results);
  int findSuboptimalResultIndex(std::vector<regResults>& result);
  void checkCorrelationFunctionFromResiduals();
  bool testBetaLowerLimit(data_t& beta);
  bool testBetaUpperLimit(data_t& beta);
};

#endif
