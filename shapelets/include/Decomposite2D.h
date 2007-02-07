#ifndef DECOMPOSITE2D_H
#define DECOMPOSITE2D_H

#include <string>
#include <Point2D.h>
#include <Grid.h>
#include <FitsImage.h>
#include <NumMatrix.h>
#include <NumVector.h>
#include <NumMatrixDiagonal.h>
#include <Object.h>

/// 2D Decomposition class.
/// Provides decomposition of a given 2D function
/// into cartesian shapelet coefficients.\n
/// The coefficients are calculated using
/// \f$coeffs = (M^T V^{-1} M)^{-1} M^T V^{-1}\cdot data\f$, where
/// \f$coeffs\f$ is a vector representation of the coefficient matrix,
/// \f$M\f$ is the matrix of the value of the shapelet basis function of every order
/// in every pixel of the image (e.g. \f$M_{2,0}\f$ is the value of the basis function 
/// \f$B_{0,0}\f$ at pixel 2 of the image), V is the covariance matrix of the pixels
/// (by now it's made up from pixel noise only) and \f$data\f$ is a vector representation
/// of the image values (see Paper III, eq. 83).\n
/// The errors of the coefficients are given by \f$(M^T V^{-1} M)^{-1}\f$.\n
///
/// Also the goodness of the fit in term of 
/// \f$\chi^2 = \frac{(data - reco)^T \cdot V^{-1} 
/// \cdot (data - reco)}{n_{pixels} - n_{coeffs}}\f$ is calculated, 
/// where \f$reco = M \cdot coeffs\f$ is the reconstructed image. (PaperIII, eq. 35).
///
/// \todo
/// - construct 2D LS matrix from 2 1D shapelet matrices (if Grid is square: only 1 matrix!)
/// - shift the grid (and centroid) to start at (0/0) -> matrices are similar at most pixel
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
  /// Set centroid position for the decomposition
  void setCentroid(const Point2D& xcentroid);
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
  NumMatrix<double> M, Mt, LS;
  Grid grid;
  Point2D xcentroid;
  void computeCoeffs();
  void computeModel();
  void computeResiduals();
  void makeLSMatrix ();
  void makeV_Matrix();
  double beta, background_variance;
  int nmax,nCoeffs,npixels;
  NumVector<double> coeffVector, errorVector, model, residual;
  const NumVector<double>& data, bg_rms;
  NumMatrixDiagonal<double> V_;
  NumMatrix<int> nVector;
  bool change, updateC, updateModel, updateResiduals, gaussian;
};

#endif
