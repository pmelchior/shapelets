#ifndef DECOMPOSITE2D_H
#define DECOMPOSITE2D_H

#include <string>
#include <NumMatrix.h>
#include <NumVector.h>
#include <NumMatrixDiagonal.h>
#include <Typedef.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <frame/Image.h>
#include <frame/Object.h>
#include <shapelets/Composite2D.h>

/// Shapelet decomposition class.
/// Provides decomposition of a given Object entity into a Composite2D entity, aka
/// a shapelet model. All protected data members of Composite2D are populated during the
/// decomposition, such that after the decomposition (which requires at least one call to 
/// getChiSquare()) one obtains a fully functional shapelet model. Thus, Decomposite2D behaves like
/// an image filter.
///
/// The coefficients are calculated as the \em Least-Squares solution of the equation
/// \f$obj = M \cdot coeffs + noise\f$,
/// where \f$M\f$ is the design matrix of the problem (e.g. \f$M_{2,0}(\beta)\f$ is the value of the 
/// basis function \f$B_{0,0}(x_2;\beta)\f$ at pixel \f$x_2\f$ of the image) which is provided by 
/// the Composite2D entity.\n
/// For obtaining the shapelet coefficients, we construct the pseudo-inverse 
/// \f$M^\dagger\equiv\bigl(M^T V^{-1} M\bigr)^{-1} M^T V^{-1}\f$ of \f$M\f$ subject to the pixel
/// covariance matrix \f$V\f$ and solve for the coefficients, \f$coeffs = M^\dagger \cdot data\f$
/// (cf. Paper III, eq. 23).\n
/// The coefficient covariance matrix is given by \f$(M^T V^{-1} M)^{-1}\f$.
///
/// Also the goodness of the fit in terms of 
/// \f[\chi^2 = \frac{(data - reco)^T \cdot V^{-1} \cdot (data - reco)}{n_{pixels} - n_{coeffs}}\f]
/// is calculated, where \f$reco = M \cdot coeffs\f$ is the reconstructed image (Paper III, eq. 18).
///
/// The error measure (which enters \f$V\f$) can be given in different ways, depending on the
/// value of ShapeLensConfig::NOISEMODEL:
/// - \p GAUSSIAN: \f$V = diag(\sigma_n^2)\f$, where \f$\sigma_n^2\f$ is obtained from 
///   Object::getNoiseRMS().
/// - \p POISSONIAN: \f$V = diag(\sigma_n^2) + diag(reco)\f$, where \f$reco\f$ is the active
///   shapelet model at the time updateWeightMap() is called.
/// - \p WEIGHT: \f$V = diag(weight)\f$, where \f$weight\f$ is the total weight (defined as its 
///   inverse variance) of each pixel obtained from Object::getWeightMap().
/// - \p WEIGHT_POISSONIAN: \f$V = diag(weight) + diag(reco)\f$, a combination of \p WEIGHT and
///   \p POISSONIAN.
/// - \p COVARIANCE: \f$V\f$ is the actual PixelCovarianceMatrix, which has to be measured 
///   from \f$data\f$ and stored in Object::getPixelCovarianceMatrix().
///
/// In addition to that, one can also explicitly set arbitrary shapelet coefficients \f$coeff'\f$
/// by calling setCoeffs() and compute the \f$\chi^2\f$ of this shapelet model 
/// \f$reco' = M\cdot coeff'\f$. For this to happen, one must explicitly fix the coefficients (in
/// contrast to computing those coefficients which minimize\f$\chi^2\f$) by calling fixCoeffs().

class Decomposite2D {
 public:
  /// Contructor for decomposing a Object into shapelets of maximum order \f$n_{max}\f$.
  Decomposite2D(const Object& obj, Composite2D& C2D);
  /// Get the current decomposition scale.
  data_t getBeta();
  /// Set scale size \f$\beta\f$ for decomposition.
  void setBeta(data_t beta);
  /// Return maximal decomposition order.
  int getNMax();
  /// Set maximal shapelet orders for decomposition.
  void setNMax(int nmax);
  /// Return \f$\chi^2\f$ of decomposition.
  data_t getChiSquare();
  /// Return the variance \f$\sigma(\chi^2)\equiv\sqrt{\frac{2}{n_{pixels} - n_{coeffs}}}\f$.
  data_t getChiSquareVariance();
  /// Update weight map.
  void updateWeightMap();
  /// Allows to fix the coefficients.
  /// With <tt>fixed == true</tt>, coefficients are not set to best fit values w.r.t. \f$\chi^2\f$.
  /// By default this is set to \p false, but can be switched on in case the coefficients
  /// are given by setCoeffs().
  void fixCoeffs(bool fixed);
  /// Set the shapelet coefficients for constructing a model.
  /// This only makes sense in conjunction with fixCoeffs() as one is then able to compute
  /// the \f$\chi^2\f$ of this particular shapelet model.\n\n
  void setCoeffs(const CoefficientVector<data_t>& coeffs);
  /// Get residuals between data and shapelet model.
  /// \f$res = data-model\f$.
  const NumVector<data_t>& getResiduals();

 protected:
  /// Reference to Composite2D entity.
  Composite2D& C2D;
  /// Reference to Object entity.
  const Object& obj;
  bool computeCoeffs();
  bool computeModel();
  bool computeResiduals();
  void makeLSMatrix ();

 private:
  data_t background_variance, chi2;
  char noise;
  NumMatrix<data_t> Mt, LS;
  NumVector<data_t> residual;
  PixelCovarianceMatrix V_;
  NumMatrixDiagonal<data_t> Weight;
  bool fixedCoeffs;
};

#endif
