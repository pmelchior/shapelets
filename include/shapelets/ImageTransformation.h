#ifndef SHAPELENS_IMAGETRANSFORMATION_H
#define SHAPELENS_IMAGETRANSFORMATION_H

#include <complex.h>
// for 3D convolution tensor
#include <boost/multi_array.hpp>
#include <numla/NumMatrix.h>
#include <numla/NumMatrixDiagonal.h>
#include "../Typedef.h"
#include "../utils/History.h"
#include "../frame/Point2D.h"
#include "CoefficientVector.h"

namespace shapelens {

/// Class for image transformation in shapelet space.
/// See Paper I, sect. 3.3 and Paper III, sect. 5 for details.\n\n
/// In general, the transformations are done by constructing a transformation matrix
/// \f$M\f$ and applying it to the coefficients, \f$I^\prime_\mathbf{n} = M\ I_\mathbf{n}\f$.\n
/// If a covariance matrix \f$C\f$ is given, it will be transformed accordingly,
/// \f$C^\prime = M\ C\ M^T\f$.\n\n
/// \b CAUTION: The infinitesimal transformation converge(), shear(), 
/// flex(), translate() and rescale() are to correct first order only.\n\n
/// If a History is passed to the method, appropriate statements are logged to it.

class ImageTransformation {
  typedef complex<data_t> Complex;
   public:
  /// Default constructor.
  ImageTransformation();

  /// Translate the image by \f$dx_1, dx_2\f$.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void translate(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t dx1, data_t dx2, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Translate the image by applying the given \p translationMatrix.
  /// See getTranslationMatrix() and translate().
  void translate(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& translationMatrix, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute translation matrix.
  /// cf. Paper I, eq. (34).
  NumMatrix<data_t> getTranslationMatrix(data_t beta, data_t dx1, data_t dx2, const IndexVector& nVector);
  
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void rotate(CoefficientVector<Complex>& polarCoeffs, data_t rho, NumMatrix<Complex>* covariance = NULL, History* history = NULL);
  /// Rotate image by applying the given \p rotationMatrix.
  /// See getRotationMatrix() and rotate().
  void rotate(CoefficientVector<Complex>& polarCoeffs, const NumMatrix<Complex>& rotationMatrix, NumMatrix<Complex>* covariance = NULL, History* history = NULL);
  /// Compute rotation matrix.
  /// cd. Paper III, eq. (38).
  NumMatrix<Complex> getRotationMatrix(data_t rho, const IndexVector& nVector);

  /// Circularize the image.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void circularize(CoefficientVector<Complex>& polarCoeffs,  NumMatrix<Complex>* covariance = NULL, History* history = NULL);
  /// Circularize the image by applying the given \p circularizationMatrix.
  /// See getCircularizationMatrix() and circularize().
  void circularize(CoefficientVector<Complex>& polarCoeffs, const NumMatrixDiagonal<Complex>& circularizationMatrix, NumMatrix<Complex>* covariance = NULL, History* history = NULL);  
  /// Compute circularization matrix
  /// cf. Paper III, eq. (44).
  NumMatrixDiagonal<Complex> getCircularizationMatrix(const IndexVector& nVector);

  
  /// Flip the image w.r.t. horizontal axis.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void flipX(CoefficientVector<Complex>& polarCoeffs,  NumMatrix<Complex>* covariance = NULL, History* history = NULL);
  /// Flip the image w.r.t. horizontal axis by applying the given \p flipMatrix.
  void flipX(CoefficientVector<Complex>& polarCoeffs,  const NumMatrixDiagonal<Complex>& flipMatrix, NumMatrix<Complex>* covariance = NULL, History* history = NULL);
  /// Compute the flip matrix.
  /// cf. Paper III, eq. (45).\n
  /// \b CAUTION: As the entries of this matrix depend on the values of 
  /// \p polarCoefficients, this matrix cannot be reused for different coefficient sets.
  NumMatrixDiagonal<Complex> getFlipMatrix(CoefficientVector<Complex>& polarCoeffs);


  /// Brighten the image by the given factor.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.\n
  /// Because of the simplicity of the transformation (transformation matrix is diagonal
  /// with the constant entry \p factor), the matrix will not be computed explicitly.
  void brighten(CoefficientVector<data_t>& cartesianCoeffs, data_t factor, NumMatrix<data_t>* covariance = NULL, History* history = NULL);

  /// Convolve the image with a kernel.
  /// The convolution raises the maximum order,
  /// \f$n_{max} \to n_{max}+n_{max}^\text{kernel}\f$, and employ the 
  /// <em>natural choice</em> \f$\beta^2 \to \beta^2 + \bigl(\beta^\text{kernel}\bigr)^2\f$ 
  /// for the convolved scale.\n\n
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Convolve the image with a kernel by applying the given \p convolutionMatrix.
  /// See getConvolutionMatrix() and convolve().
  void convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& convolutionMatrix, data_t beta_kernel, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute convolution matrix.
  /// cf. Paper I, eq. (52), and Paper II, sect. 3.1.\n
  /// For tweaking issues, maximum orders and scale sizes of all objects are
  /// configurable (but must be handled with care).
  NumMatrix<data_t> getConvolutionMatrix(const CoefficientVector<data_t>& kernelCoeffs, unsigned int nmax_orig, unsigned int nmax_kernel, unsigned int nmax_convolved, data_t beta_orig, data_t beta_kernel, data_t beta_convolved);

  /// Deconvolve the image from a kernel.
  /// cf. Paper IV, sect. 3.1.\n
  /// The deconvolution lowers the maximum order, \f$n_{max} \to n_{max}-n_{max}^\text{kernel}\f$, and the scale size,
  /// \f$\beta^2 \to \beta^2 - \bigl(\beta^\text{kernel}\bigr)^2\f$ 
  /// if \f$\beta>\beta^\text{kernel}\f$ or \f$\beta\to 0.1\f$ else.\n
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Deconvolve the image from a kernel by applying the given \p deconvolutionMatrix.
  /// See getDeconvolutionMatrix() and deconvolve().
  void deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& deconvolutionMatrix, data_t beta_kernel, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute deconvolution matrix.
  /// cf. Paper IV, sect. 3.1 and deconvolve().\n
  /// The computed matrix is the pseudo-inverse 
  /// \f$\bigl(P^t P\bigr)^{-1} P^t\f$ of getConvolutionMatrix().\n
  /// If \p covariance is given, it will be considered in the computation.
  NumMatrix<data_t> getDeconvolutionMatrix(const CoefficientVector<data_t>& kernelCoeffs, unsigned int nmax_orig, unsigned int nmax_kernel, unsigned int nmax_convolved, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, NumMatrix<data_t>* covariance = NULL);

  /// Change the scale size of the image by using rescaling relation.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void rescale(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, data_t newbeta, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Change the scale size of the image by applying the given \p rescalingMatrix.
  /// See getRescalingMatrix() and rescale().
  /// \b CAUTION: In addition, the scale size of the object has to be changed here.
  void rescale(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& rescalingMatrix, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute rescaling matrix.
  /// cf. Paper I, appendix A.
  NumMatrix<data_t> getRescalingMatrix(data_t beta, data_t newbeta, const IndexVector& nVector);


  /// Apply convergence to image.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void converge(CoefficientVector<data_t>& cartesianCoeffs, data_t kappa, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Apply convergence to image by employing the given \p convergenceMatrix.
  /// See getConvergenceMatrix() and converge().
  void converge(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& convergenceMatrix, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute convergence matrix.
  /// cf. Paper I, eq. (32).
  NumMatrix<data_t> getConvergenceMatrix(data_t kappa, const IndexVector& nVector);

  /// Apply shear to the image.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void shear(CoefficientVector<data_t>& cartesianCoeffs, complex<data_t> gamma, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Apply shear to image by employing the given \p shearMatrix.
  /// See getShearMatrix() and shear().
  void shear(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& shearMatrix, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute shear matrix.
  /// cf. Paper I, eq. (32).
  NumMatrix<data_t> getShearMatrix(complex<data_t> gamma, const IndexVector& nVector);

  /// Apply flexion to the image.
  /// If \p covariance is given, the covariance matrix is updated. 
  /// If \p history is given, an appropriate statement is appended to it.
  void flex(CoefficientVector<data_t>& cartesianCoeffs, complex<data_t> F, complex<data_t> G, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Apply flexion to image by employing the given \p flexionMatrix.
  /// See getFlexionMatrix() and flex().
  void flex(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& flexionMatrix, NumMatrix<data_t>* covariance = NULL, History* history = NULL);
  /// Compute flexion matrix.
  /// See PM's diploma thesis, eq. (5.38) - (5.41)
  NumMatrix<data_t> getFlexionMatrix(complex<data_t> F, complex<data_t> G, const IndexVector& nVector);

 private:
  void makeBTensor(boost::multi_array<data_t,3>& bt, data_t alpha_1, data_t beta_1, data_t gamma_1, int nmax);
  NumMatrix<data_t> make1DRescalingMatrix(data_t beta1, data_t beta2, int nmax);
};
} // end namespace
#endif
