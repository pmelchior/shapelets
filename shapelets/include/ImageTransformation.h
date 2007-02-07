#ifndef IMAGETRANSFORMATION_H
#define IMAGETRANSFORMATION_H

#include <NumMatrix.h>
#include <complex.h>
#include <iostream>
#include <Point2D.h>
// for 3D convolution tensor
#include <boost/multi_array.hpp>

/// Class for image transformation in shapelet space.
/// see Paper III, chapter 5 for details
/// All transformations are logged in history.

class ImageTransformation {
  typedef complex<double> Complex;
   public:
  /// Default constructor.
  ImageTransformation();
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(NumMatrix<Complex>& polarCoeffs, double rho, std::ostringstream& history);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(NumMatrix<Complex>& polarCoeffs, double& beta, double kappa, std::ostringstream& history);
  /// Shear the image by \f$\gamma_0, \gamma_1\f$.
  /// If the shapelet order is lower than 10, the order will be increased by 2;
  /// otherwise the action of the shear would not be correctly described.
  void shear(NumMatrix<Complex>& polarCoeffs, double gamma0, double gamma1, std::ostringstream& history);
  /// Apply flexion to the image.
  /// The flexion is specified by giving the derivatives of the shear 
  /// \f$\partial_i \gamma_j\f$.\n
  /// see Bacon et al. (2005), arxiv: astro-ph/0504478 and 
  /// Goldberg/Bacon, ApJ 619, 741, 2005.\n
  /// WARNING: The operators in the paper are wrong, implementation uses correct form.
  void flex(NumMatrix<double>& cartesianCoeffs, NumMatrix<double>& dGamma, std::ostringstream& history);
  /// Translate the image by \f$dx0, dx1\f$.
  /// This is not done by transformation of the coefficients, 
  /// but by shifting the centroid poisition.
  void translate(Point2D& xcentroid, double dx0, double dx1, std::ostringstream& history);
  /// Circularize the image.
  /// The image will be averaged over the polar angle and then have a radial dependency only.
  void circularize(NumMatrix<Complex>& polarCoeffs, std::ostringstream& history);
  /// Flip the image along the X axis.
  void flipX(NumMatrix<Complex>& polarCoeffs, std::ostringstream& history);
  /// Brighten the image by the given factor.
  void brighten(NumMatrix<double>& cartesianCoeffs, NumMatrix<Complex>& polarCoeffs, double factor, std::ostringstream& history);
  void convolve(NumMatrix<double>& cartesianCoeffs, double& beta, NumMatrix<double>& KernelCoeffs, double beta_kernel, std::ostringstream& history);  
  void deconvolve(NumMatrix<double>& cartesianCoeffs, double& beta, NumMatrix<double>& KernelCoeffs, double beta_kernel, std::ostringstream& history);
  /// Change the scale size of the image by using rescaling relation.
  /// see Paper I, appendix A.
  void rescale(NumMatrix<double>& cartesianCoeffs, double beta, double newbeta, std::ostringstream& history);
  void makeRescalingMatrix(NumMatrix<double>& betaTrafo, double beta1, double beta2, int nCoeffs, NumMatrix<int>& nVector);

 private:
  void makeConvolutionMatrix(NumMatrix<double>& P, NumMatrix<double>& KernelCoeffs, double beta_orig, double beta_kernel, double beta_convolved, int nmax_orig, int nmax_convolved);
  void makeBTensor(boost::multi_array<double,3>& bt, double alpha_1, double beta_1, double gamma_1, int nmax);
};

#endif
