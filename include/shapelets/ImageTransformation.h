#ifndef IMAGETRANSFORMATION_H
#define IMAGETRANSFORMATION_H

#include <complex.h>
// for 3D convolution tensor
#include <boost/multi_array.hpp>
#include <NumMatrix.h>
#include <Typedef.h>
#include <History.h>
#include <frame/Point2D.h>
#include <shapelets/CoefficientVector.h>

/// Class for image transformation in shapelet space.
/// See Paper I, sect. 3.3 and Paper III, sect. 5 for details.\n\n
/// \b CAUTION: The infinitesimal transformation converge(), shear(), 
/// flex(), translate() and rescale() are to correct first order only.\n\n
/// All transformations are logged in History.

class ImageTransformation {
  typedef complex<data_t> Complex;
   public:
  /// Default constructor.
  ImageTransformation();
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(CoefficientVector<Complex>& polarCoeffs, data_t rho, History& history);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(CoefficientVector<Complex>& polarCoeffs, data_t kappa, History& history);
  /// Shear the image by complex \f$\gamma\f$.
  void shear(CoefficientVector<Complex>& polarCoeffs, complex<data_t> gamma, History& history);
  /// Apply flexion to the image.
  /// The flexion is specified by giving the derivatives of the shear 
  /// \f$\partial_i \gamma_j\f$.\n
  void flex(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& dGamma, History& history);
  /// Translate the image by \f$dx0, dx1\f$.
  void translate(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t dx0, data_t dx1, History& history);
  /// Circularize the image.
  /// The image will be averaged over the polar angle and then have a radial dependency only.
  void circularize(CoefficientVector<Complex>& polarCoeffs, History& history);
  /// Flip the image along the X axis.
  void flipX(CoefficientVector<Complex>& polarCoeffs, History& history);
  /// Brighten the image by the given factor.
  void brighten(CoefficientVector<data_t>& cartesianCoeffs, data_t factor, History& history);
  void convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_kernel, History& history);  
  void deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_kernel, History& history);
  /// Change the scale size of the image by using rescaling relation.
  /// see Paper I, appendix A.
  void rescale(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t newbeta, History& history);
  NumMatrix<data_t> makeRescalingMatrix(data_t beta1, data_t beta2, const IndexVector& nVector);
  NumMatrix<data_t> makeConvolutionMatrix(const CoefficientVector<data_t>& cartesianCoeffs, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, unsigned int nmax_convolved);
 private:
  void makeBTensor(boost::multi_array<data_t,3>& bt, data_t alpha_1, data_t beta_1, data_t gamma_1, int nmax);
  NumMatrix<data_t> make1DRescalingMatrix(data_t beta1, data_t beta2, int nmax);
};

#endif
