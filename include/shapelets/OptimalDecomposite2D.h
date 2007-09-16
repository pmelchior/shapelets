#ifndef OPTIMAL_DECOMPOSITE2D_H
#define OPTIMAL_DECOMPOSITE2D_H

#include <gsl/gsl_vector.h>
#include <map>
#include <NumMatrix.h>
#include <NumVector.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <frame/Object.h>
#include <frame/History.h>
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
  OptimalDecomposite2D(const Object& obj, int nmaxLow, int nmaxHigh, double betaLow, double betaHigh);
  /// Employs a regularization to lower the negative flux.
  /// The minimization ends, whenever \f$R=- F^-/F^+ < wantedR\f$, where \f$F^{\pm}\f$
  /// is the positive or negative flux of the shapelet reconstruction. 
  /// The return value is the actual \f$R\f$.\n
  /// If wantedR is chosen too low (e.g. \f$10^{-6}\f$, in some cases already 
  /// \f$10^{-5}\f$), the minimization procedure
  /// must increase \f$n_{max}\f$ and can therefore take a considerable amount of time.
  double regularize(double wantedR);
  /// Deliver best fit shapelet coefficients.
  void getShapeletCoeffs(NumMatrix<double>& coeffs);
   /// Deliver best fit shapelet coefficient errors.
  void getShapeletErrors(NumMatrix<double>& errors);
  /// Get best fit residuals of the decomposition.
  const NumVector<double>& getResiduals();
  /// Get best fit shapelet model.
  const NumVector<double>&  getModel();
  /// Return best fit  \f$\beta\f$
  double getOptimalBeta();
  /// Return best fit \f$\chi^2\f$.
  double getOptimalChiSquare();
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
  double beta,optimalBeta,optimalChiSquare,bestChiSquare,image_dimension, betaHigh, betaLow;
  bool optimized, nmaxTrouble, noise_correlated;
  char flag, comp_corr;
  History history;
  std::ostringstream text;
  std::string comp_corr_string;
  std::map<int, double> bestBeta, bestChi2;
  struct regResults {
    int nmax;
    double f;
    double lambda;
    double beta;
    double chi2;
    double R;
    NumVector<double> coeffs;
  };
  int findOptimalBeta(unsigned char step);
  void findOptimalNMax(unsigned char step), optimize(), checkBeta();
  void getCoeffErrorFromBeta(const NumVector<double>& coeffVector, NumVector<double>& errorVector);
  void getBetaTrafoMatrix(NumMatrix<double>& betaTrafo, double beta1, double beta2);
  void appendRegResults(std::vector<regResults>& results, int nmax, double f, double lambda, double beta, double chi2, double R, const NumVector<double>& coeffs);
  int findNMaxofBestF(std::vector<regResults>& results);
  int findSuboptimalResultIndex(std::vector<regResults>& result);
  void checkCorrelationFunctionFromResiduals();
  bool testBetaLowerLimit(double& beta);
  bool testBetaUpperLimit(double& beta);
};

#endif
