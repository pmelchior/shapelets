#ifndef IMAGETRANSFORMATION_H
#define IMAGETRANSFORMATION_H

#include <complex.h>
// for 3D convolution tensor
#include <boost/multi_array.hpp>
#include <NumMatrix.h>
#include <Typedef.h>
#include <History.h>
#include <frame/Point2D.h>
#include <shapelets/IndexVector.h>

/// Class for image transformation in shapelet space.
/// See Paper I, sect. 3.3 and Paper III, sect. 5 for details.\n
/// All transformations are logged in History.

class ImageTransformation {
  typedef complex<data_t> Complex;
   public:
  /// Default constructor.
  ImageTransformation();
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(NumMatrix<Complex>& polarCoeffs, data_t rho, History& history);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(NumMatrix<Complex>& polarCoeffs, data_t& beta, data_t kappa, History& history);
  /// Shear the image by \f$\gamma_0, \gamma_1\f$.
  /// If the shapelet order is lower than 10, the order will be increased by 2;
  /// otherwise the action of the shear would not be correctly described.
  void shear(NumMatrix<Complex>& polarCoeffs, data_t gamma0, data_t gamma1, History& history);
  /// Apply flexion to the image.
  /// The flexion is specified by giving the derivatives of the shear 
  /// \f$\partial_i \gamma_j\f$.\n
  /// see Bacon et al. (2005), arxiv: astro-ph/0504478 and 
  /// Goldberg/Bacon, ApJ 619, 741, 2005.\n
  /// WARNING: The operators in the paper are wrong, implementation uses correct form.
  void flex(NumMatrix<data_t>& cartesianCoeffs, const NumMatrix<data_t>& dGamma, History& history);
  /// Translate the image by \f$dx0, dx1\f$.
  /// This is only valid for small translations (below 1 pixel).
  void translate(NumMatrix<data_t>& cartesianCoeffs, data_t beta, data_t dx0, data_t dx1, History& history);
  /// Circularize the image.
  /// The image will be averaged over the polar angle and then have a radial dependency only.
  void circularize(NumMatrix<Complex>& polarCoeffs, History& history);
  /// Flip the image along the X axis.
  void flipX(NumMatrix<Complex>& polarCoeffs, History& history);
  /// Brighten the image by the given factor.
  void brighten(NumMatrix<data_t>& cartesianCoeffs, NumMatrix<Complex>& polarCoeffs, data_t factor, History& history);
  void convolve(NumMatrix<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel, History& history);  
  void deconvolve(NumMatrix<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel, History& history);
  /// Change the scale size of the image by using rescaling relation.
  /// see Paper I, appendix A.
  void rescale(NumMatrix<data_t>& cartesianCoeffs, data_t beta, data_t newbeta, History& history);
  void makeRescalingMatrix(NumMatrix<data_t>& betaTrafo, data_t beta1, data_t beta2, const IndexVector& nVector);

 private:
  void makeConvolutionMatrix(NumMatrix<data_t>& P, const NumMatrix<data_t>& KernelCoeffs, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, int nmax_orig, int nmax_convolved);
  void makeBTensor(boost::multi_array<data_t,3>& bt, data_t alpha_1, data_t beta_1, data_t gamma_1, int nmax);
  void make1DRescalingMatrix(NumMatrix<data_t>& M1D, data_t beta1, data_t beta2, int nmax);
};

#endif
