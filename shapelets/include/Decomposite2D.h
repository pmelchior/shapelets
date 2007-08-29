#ifndef DECOMPOSITE2D_H
#define DECOMPOSITE2D_H

#include <string>
#include <Point2D.h>
#include <Grid.h>
#include <Image.h>
#include <NumMatrix.h>
#include <NumVector.h>
#include <NumMatrixDiagonal.h>
#include <IndexVector.h>
#include <Object.h>

/// 2D Decomposition class.
/// Provides decomposition of a given 2D function
/// into cartesian shapelet coefficients.\n
/// The coefficients are calculated using
/// \f$coeffs = (M^T V^{-1} M)^{-1} M^T V^{-1}\cdot data\f$, where
/// \f$coeffs\f$ is a vector representation of the coefficient matrix,
/// \f$M\f$ is the matrix of the value of the shapelet basis function of every order
/// in every pixel of the image (e.g. \f$M_{2,0}\f$ is the value of the basis function 
/// \f$B_{0,0}\f$ at pixel 2 of the image), \f$V\f$ is the covariance matrix of the pixels
/// and \f$data\f$ is a vector representation
/// of the image values (cf. Paper III, eq. 23).\n
/// The errors of the coefficients are given by \f$(M^T V^{-1} M)^{-1}\f$.\n
///
/// Also the goodness of the fit in term of 
/// \f$\chi^2 = \frac{(data - reco)^T \cdot V^{-1} 
/// \cdot (data - reco)}{n_{pixels} - n_{coeffs}}\f$ is calculated, 
/// where \f$reco = M \cdot coeffs\f$ is the reconstructed image. (Paper III, eq. 18).
///
/// The error measure (which enters \f$V\f$) can be given in different ways, depending on the
/// value of Object::getNoiseModel():
/// - <tt>GAUSSIAN</tt>: \f$V = diag(\sigma_n^2)\f$
/// - <tt>POISSONIAN</tt>: \f$V = diag(\sigma_n^2) + diag(data')\f$, where \f$data'\f$ 
///   is a vector representation of a smoothed version of \f$data\f$
/// - <tt>COVARIANCE</tt>: \f$V\f$ is the actual PixelCovarianceMatrix, measured 
///   directly from \f$data\f$
/// - <tt>WEIGHT</tt>: \f$V = diag(weight)\f$, where \f$weight\f$ is a vector representation
///   of the error in each pixel obtained from Object::getBackgroundRMSMap()
///
/// \todo
/// - construct 2D LS matrix from 2 1D shapelet matrices (if Grid is square: only 1 matrix!)
/// - shift the grid (and centroid) to start at (0/0) -> matrices are similar at most pixel
/// - POISSONIAN: construct V from shapelet model instead of data

class Decomposite2D {
 public:
  /// Contructor for decomposing a Object into shapelets of maximum order \f$n_{max}\f$.
  Decomposite2D(int nmax, double beta, const Object& obj);
  /// Get the decomposition shapelet coefficients as matrix.
  const NumVector<double>& getCoeffs();
  /// Access the decomposition shapelet coefficients as matrix.
  NumVector<double>& accessCoeffs();
  /// Allows to update coeffs to best fit values w.r.t. \f$\chi^2\f$.
  /// Is set to 1 by default, but can be switched of in case the coefficients
  /// are given from outside.
  void updateCoeffs(bool update);
  /// Get error matrix of the shapelet coefficients.
  const NumVector<double>& getErrors();
  /// Get shapelet model.
  /// This is the model reconstructed from the coefficients provided by getCoeffs().
  const NumVector<double>& getModel();
  /// Access shapelet model.
  /// The model will not be updated automatically. If you want to use best-fit values,
  /// call getChiSquare() before.
  NumVector<double>& accessModel();
  /// Get residuals between data and shapelet model.
  /// \f$res = data-model\f$.
  const NumVector<double>& getResiduals();
  /// Access residuals.
  /// The residuals will not be updated automatically. If you want to use best-fit values,
  /// call getChiSquare() before.
  NumVector<double>& accessResiduals();
  /// Set scale size \f$\beta\f$ for decomposition.
  void setBeta(double beta);
  /// Set maximal shapelet orders for decomposition.
  void setNMax(int nmax);
  /// Return \f$\chi^2\f$ of decomposition.
  double getChiSquare();
  /// Return the variance \f$\sigma(\chi^2)\f$.
  double getChiSquareVariance();
  /// Return maximal decomposition order in each direction.
  int getNMax();
  /// Update model and residuals during the next call to
  /// getModel() or accessModel().
  void updateModelResiduals();

 private:
  void computeCoeffs();
  void computeModel();
  void computeResiduals();
  void makeLSMatrix ();
  NumMatrix<double> M, Mt, LS;
  //Grid grid;
  //Point2D xcentroid;
  double beta, background_variance;
  int nmax,nCoeffs,npixels;
  char noise;
  NumVector<double> coeffVector, errorVector, model, residual;
  const Object& obj;
  NumMatrixDiagonal<double> Weight;
  PixelCovarianceMatrix V_;
  IndexVector nVector;
  bool change, updateC, updateModel, updateResiduals, gaussian;
};

#endif
