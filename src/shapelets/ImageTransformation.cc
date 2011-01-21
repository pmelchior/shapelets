#include "../../include/shapelets/ImageTransformation.h"
#include "../../include/shapelets/CoefficientVector.h"
#include <numla/NumVector.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace shapelens;
using std::endl;
typedef std::complex<data_t> Complex;
const Complex I = Complex(0,1);


ImageTransformation::ImageTransformation() {
}

NumMatrix<data_t> ImageTransformation::getTranslationMatrix(data_t beta, data_t dx1, data_t dx2, const IndexVector& nVector) {
  // rescale dx1 and dx2 to be in units of beta and multiply with
  // ubiquious sqrt(1/2)
  dx1 *= M_SQRT1_2/beta;
  dx2 *= M_SQRT1_2/beta;
  // set up translation matrix for this IndexVector
  NumMatrix<data_t> M(nVector.getNCoeffs(), nVector.getNCoeffs());
  int n1,n2,j;
  for (unsigned int i=0; i<nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    M(i,i) = 1;
    if (n1>=1) {
      j = nVector.getIndex(n1-1,n2);
      M(i,j) += dx1*sqrt((data_t)n1);
    }
    if (n2>=1) {
      j = nVector.getIndex(n1,n2-1);
      M(i,j) += dx2*sqrt((data_t)n2);
    }
    if (n1+n2 <= -1 + nVector.getNMax()) {
      j = nVector.getIndex(n1+1,n2);
      M(i,j) -= dx1*sqrt((data_t)n1+1);
      j = nVector.getIndex(n1,n2+1);
      M(i,j) -= dx2*sqrt((data_t)n2+1);
    }
  }
  return M;
}

void ImageTransformation::translate(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t dx1, data_t dx2, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Translating image by " << dx1 << "/" << dx2 << endl;
  NumMatrix<data_t> M = getTranslationMatrix(beta,dx1,dx2,cartesianCoeffs.getIndexVector());
  // perform transformation
  cartesianCoeffs = M*cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M * (*cov)) * M.transpose();
}

void ImageTransformation::translate(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& M, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Translating image via given transformation matrix" << endl;
  // perform transformation
  cartesianCoeffs = M*cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M * (*cov)) * M.transpose();
}

NumMatrix<Complex> ImageTransformation::getRotationMatrix(data_t rho, const IndexVector& nVector) {
  data_t rho_scaled = M_PI*rho/180;
  // set up rotation matrix for this IndexVector
  NumMatrix<Complex> M(nVector.getNCoeffs(), nVector.getNCoeffs());
  int m;
  for (unsigned int i=0; i<nVector.getNCoeffs(); i++) {
    m = nVector.getState2(i);
    M(i,i) = exp(m*rho_scaled*I);
  }
  return M;
}

void ImageTransformation::rotate(CoefficientVector<Complex>& polarCoeffs, data_t rho, NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Rotating image by " << rho << " degrees" << endl;
  NumMatrix<Complex> M = getRotationMatrix(rho,polarCoeffs.getIndexVector());
  // perform transformation
  polarCoeffs = M*polarCoeffs;
  if (cov != NULL && cov->getColumns() == polarCoeffs.getNCoeffs())
    *cov =(M * (*cov)) * M; // M is diagonal, no need to transpose
}

void ImageTransformation::rotate(CoefficientVector<Complex>& polarCoeffs, const NumMatrix<Complex>& M, NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Rotating image via given transformation matrix" << endl;
  // perform transformation
  polarCoeffs = M*polarCoeffs;
  if (cov != NULL && cov->getColumns() == polarCoeffs.getNCoeffs())
    *cov = (M * (*cov)) * M;
}

NumMatrixDiagonal<Complex> ImageTransformation::getCircularizationMatrix(const IndexVector& nVector) {
  NumMatrixDiagonal<Complex> M(nVector.getNCoeffs());
   int m;
  for (unsigned int i=0; i<nVector.getNCoeffs(); i++) {
    m = nVector.getState2(i);
    if (m == 0)
      M(i) = Complex(1,0);
  }
  return M;
}

void ImageTransformation::circularize(CoefficientVector<Complex>& polarCoeffs, NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Circularizing image" << endl;
  NumMatrixDiagonal<Complex> M = getCircularizationMatrix(polarCoeffs.getIndexVector());
  // perform transformation
  polarCoeffs = M*polarCoeffs;
  // as M is diagonal use a simpler form of M*cov*M^T
  if (cov != NULL) {
    for (unsigned int i=0; i<cov->getRows(); i++)
      for (unsigned int j=0; j<cov->getColumns(); j++)
	(*cov)(i,j) *= M(i)*M(j);
  }
}

void ImageTransformation::circularize(CoefficientVector<Complex>& polarCoeffs, const NumMatrixDiagonal<Complex>& M, NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Circularizing image" << endl;
  polarCoeffs = M*polarCoeffs;
  // as M is diagonal use a simpler form of M*cov*M^T
  if (cov != NULL) {
    for (unsigned int i=0; i<cov->getRows(); i++)
      for (unsigned int j=0; j<cov->getColumns(); j++)
	(*cov)(i,j) *= M(i)*M(j);
  }
}

NumMatrixDiagonal<Complex> ImageTransformation::getFlipMatrix(CoefficientVector<Complex>& polarCoeffs) {
  NumMatrixDiagonal<Complex> M(polarCoeffs.getNCoeffs());
  for (unsigned int i=0; i<polarCoeffs.getNCoeffs(); i++)
    M(i) = exp(-2*arg(polarCoeffs(i))*I);
  return M;
}

void ImageTransformation::flipX(CoefficientVector<Complex>& polarCoeffs,  NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Flipping image arround its X-axis" << endl;
  NumMatrixDiagonal<Complex> M = getFlipMatrix(polarCoeffs);
  // perform transformation
  polarCoeffs = M*polarCoeffs;
  if (cov != NULL) {
    for (unsigned int i=0; i<cov->getRows(); i++)
      for (unsigned int j=0; j<cov->getColumns(); j++)	
	(*cov)(i,j) *= M(i)*M(j);
  }
}

void ImageTransformation::flipX(CoefficientVector<Complex>& polarCoeffs,  const NumMatrixDiagonal<Complex>& M, NumMatrix<Complex>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Flipping image arround its X-axis" << endl;
  // perform transformation
  polarCoeffs = M*polarCoeffs;
  if (cov != NULL) {
    for (unsigned int i=0; i<cov->getRows(); i++)
      for (unsigned int j=0; j<cov->getColumns(); j++)	
	(*cov)(i,j) *= M(i)*M(j);
  }
}

void ImageTransformation::brighten(CoefficientVector<data_t>& cartesianCoeffs, data_t factor, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Changing image brightness by the factor " << factor << endl;
  for (unsigned int i=0; i < cartesianCoeffs.size(); i++)
    cartesianCoeffs(i) *= factor;
  if (cov != NULL) {
    for (unsigned int i=0; i<cov->getRows(); i++)
      for (unsigned int j=0; j<cov->getColumns(); j++)	
	(*cov)(i,j) *= factor*factor;
  }
}

NumMatrix<data_t> ImageTransformation::getConvolutionMatrix(const CoefficientVector<data_t>& kernelCoeffs, unsigned int nmax_orig, unsigned int nmax_kernel, unsigned int nmax_convolved, data_t beta_orig, data_t beta_kernel, data_t beta_convolved) {
  unsigned int nmax = GSL_MAX_INT(nmax_orig,nmax_convolved);
  nmax = GSL_MAX_INT(nmax_kernel,nmax);
  boost::multi_array<data_t,3> bt(boost::extents[nmax+1][nmax+1][nmax+1]);
  data_t alpha = beta_orig;
  data_t beta = beta_kernel;
  data_t gamma = beta_convolved;
  makeBTensor(bt,1./gamma,1./alpha,1./beta,nmax);
 
  // the 1D convolution tensor C_nml
  boost::multi_array<data_t,3> c1(boost::extents[nmax_convolved+1][nmax_orig+1][nmax_kernel+1]);
  for (int n=0; n <= nmax_convolved; n++)
    for (int m=0; m <= nmax_orig; m++)
      for (int l=0; l <= nmax_kernel; l++)
	if ((n+m+l)%2 == 0)
	  c1[n][m][l] = M_SQRT2 * M_SQRTPI * gsl_pow_int(-1,(3*n + m + l)/2) * bt[n][m][l];
  
  // the 2D convolution tensor C_nml
  // the matrix->vector projection
  IndexVectorCartesian  nVector_convolved(nmax_convolved), nVector_orig(nmax_orig);
  const IndexVector& nVector_kernel = kernelCoeffs.getIndexVector();
  int nCoeffs_orig = nVector_orig.getNCoeffs();
  int nCoeffs_kernel = nVector_kernel.getNCoeffs();
  int nCoeffs_convolved = nVector_convolved.getNCoeffs();
  boost::multi_array<data_t,3> c2(boost::extents[nCoeffs_convolved][nCoeffs_orig][nCoeffs_kernel]);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++) 
      for (int l=0; l < nCoeffs_kernel; l++)
	c2[n][m][l] = 
	  c1[nVector_convolved.getState1(n)][nVector_orig.getState1(m)][nVector_kernel.getState1(l)] * 
	  c1[nVector_convolved.getState2(n)][nVector_orig.getState2(m)][nVector_kernel.getState2(l)];
  
  // now the 2D convolution matrix P_nm
  NumMatrix<data_t> P (nCoeffs_convolved,nCoeffs_orig);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++)
      for (int l=0; l < nCoeffs_kernel; l++)
	P(n,m) += c2[n][m][l] * kernelCoeffs(l);
  return P;
}

void ImageTransformation::convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Convolving image with kernel of order " << kernelCoeffs.getNMax() << ", beta = " << beta_kernel << endl;
  
  // natural choice for convolved beta
  data_t beta_convolved = sqrt(beta*beta + beta_kernel*beta_kernel);
  // theoretical nmax after convolution
  int nmax_convolved = cartesianCoeffs.getNMax() + kernelCoeffs.getNMax();

  // construct convolution matrix
  NumMatrix<data_t> P = getConvolutionMatrix(kernelCoeffs,cartesianCoeffs.getNMax(),kernelCoeffs.getNMax(),nmax_convolved,beta,beta_kernel,beta_convolved);
  // perform the convolution
  cartesianCoeffs = P*cartesianCoeffs;
  if (cov != NULL && cov->getRows() == P.getColumns())
    *cov = (P*(*cov))*P.transpose();
  beta = beta_convolved;
}

void ImageTransformation::convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& P, data_t beta_kernel, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Convolving image with kernel via given convolution matrix, beta = " << beta_kernel << endl;
  // perform the convolution
  cartesianCoeffs = P*cartesianCoeffs;
  if (cov != NULL && cov->getRows() == P.getColumns())
    *cov = (P*(*cov))*P.transpose();
  beta = sqrt(beta*beta + beta_kernel*beta_kernel);
} 

NumMatrix<data_t> ImageTransformation::getDeconvolutionMatrix(const CoefficientVector<data_t>& kernelCoeffs, unsigned int nmax_orig, unsigned int nmax_kernel, unsigned int nmax_convolved, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, NumMatrix<data_t>* cov) {
  // construct pseudoinverse of convolution matrix
  NumMatrix<data_t> P = getConvolutionMatrix(kernelCoeffs,nmax_orig,nmax_kernel,nmax_convolved,beta_orig,beta_kernel,beta_convolved);
  if (cov == NULL)
    return (P.transpose() * P).invert() * P.transpose();
  else {
    NumMatrix<data_t> Pt = P.transpose() * (*cov);
    return (Pt*P).invert()*Pt;
  } 
}

void ImageTransformation::deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Deconvolving image with kernel of order " << kernelCoeffs.getNMax() << ", beta = " << beta_kernel << endl;
  // invert 'natural choice'
  data_t beta_orig;
  if (beta_kernel < beta)
    beta_orig = sqrt(beta*beta - beta_kernel*beta_kernel);
  else
    beta_orig = 0.1;
  int nmax_orig = cartesianCoeffs.getNMax();
  // coefficient reduction here:
  // if (cartesianCoeffs.getNMax() >= kernelCoeffs.getNMax())
  //    nmax_orig = cartesianCoeffs.getNMax() - kernelCoeffs.getNMax();
  //   else
  //    nmax_orig = 0;
  NumMatrix<data_t> P_1 = getDeconvolutionMatrix(kernelCoeffs,nmax_orig,kernelCoeffs.getNMax(),cartesianCoeffs.getNMax(),beta_orig,beta_kernel,beta,cov);
  // perform the deconvolution
  cartesianCoeffs = P_1*cartesianCoeffs;
  if (cov != NULL && cov->getRows() == P_1.getColumns())
    *cov = (P_1*(*cov))*P_1.transpose();
  beta = beta_orig;
}

void ImageTransformation::deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& P_1, data_t beta_kernel, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Deconvolving image with kernel via given deconvolution matrix, beta = " << beta_kernel << endl;
  // invert 'natural choice'
  data_t beta_orig;
  if (beta_kernel < beta)
    beta_orig = sqrt(beta*beta - beta_kernel*beta_kernel);
  else
    beta_orig = 0.1;
  // perform the deconvolution
  cartesianCoeffs = P_1*cartesianCoeffs;
  if (cov != NULL && cov->getRows() == P_1.getColumns())
    *cov = (P_1*(*cov))*P_1.transpose();
  beta = beta_orig;
}

// the 2D rescaling matrix is obtained by a tensor multiplications of
// the 1D rescaling matrix with itself.
NumMatrix<data_t> ImageTransformation::getRescalingMatrix(data_t beta2, data_t beta1, const IndexVector& nVector) {
  int nmax = nVector.getNMax();
  int nCoeffs = nVector.getNCoeffs();
  // compute 1D rescaling matrix
  NumMatrix<data_t> M1D = make1DRescalingMatrix(beta2,beta1,nmax);
  // now build tensor product of 1D matrix to give 2D rescaling matrix
  NumMatrix<data_t> M2D(nCoeffs,nCoeffs);
  int i1,i2,j1,j2;
  for (int l = 0; l < nCoeffs; l++) {
    i1 = nVector.getState1(l);
    i2 = nVector.getState2(l);
    // here we assume that both sets of coefficients have same length
    for (int i=0; i< nCoeffs; i++) {
       j1 = nVector.getState1(i);
       j2 = nVector.getState2(i);
       M2D(l,i) = M1D(i1,j1)*M1D(i2,j2);
    }
  }
  return M2D;
}

void ImageTransformation::rescale(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, data_t newbeta, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Rescaling image from beta = "<< beta << " to new beta = " << newbeta << endl;
  NumMatrix<data_t> R = getRescalingMatrix(newbeta,beta,cartesianCoeffs.getIndexVector());
  cartesianCoeffs = R * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (R*(*cov))*R.transpose();
  beta = newbeta;
}
 
void ImageTransformation::rescale(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& R, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Rescaling image via given rescaling matrix" << endl;
  cartesianCoeffs = R * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (R*(*cov))*R.transpose();
}

// ***************** lensing transformations ***************** //

NumMatrix<data_t> ImageTransformation::getConvergenceMatrix(data_t kappa, const IndexVector& nVector) {
  data_t factor = -0.5*kappa; // minus sign as we want to have the inverse transformation
  // set up matrix
  NumMatrix<data_t> M(nVector.getNCoeffs(), nVector.getNCoeffs());
  int n1,n2,j;
  for (unsigned int i=0; i<nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    M(i,i) = 1+kappa;
    if (n1>=2) {
      j = nVector.getIndex(n1-2,n2);
      M(i,j) -= factor*sqrt((data_t)n1*(n1-1));
    }
    if (n2>=2) {
      j = nVector.getIndex(n1,n2-2);
      M(i,j) -= factor*sqrt((data_t)n2*(n2-1));
    }
    if (n1+n2 <= nVector.getNMax() - 2) {
      j = nVector.getIndex(n1+2,n2);
      M(i,j) += factor*sqrt((data_t)(n1+1)*(n1+2));
      j = nVector.getIndex(n1,n2+2);
      M(i,j) += factor*sqrt((data_t)(n2+1)*(n2+2));
    }
  }
  return M;
}
void ImageTransformation::converge(CoefficientVector<data_t>& cartesianCoeffs, data_t kappa, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Applying convergence kappa = " << kappa << endl;
  NumMatrix<data_t> M = getConvergenceMatrix(kappa,cartesianCoeffs.getIndexVector());
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

void ImageTransformation::converge(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& M, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Applying convergence via given transformation matrix" << endl;
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

NumMatrix<data_t> ImageTransformation::getShearMatrix(Complex gamma, const IndexVector& nVector) {
  data_t gamma_1 = -0.5*real(gamma), gamma_2 = -imag(gamma); // minus signs = inverse trafo!
  // set up matrix
  NumMatrix<data_t> M(nVector.getNCoeffs(), nVector.getNCoeffs());
  int n1, n2,j;
  for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    M(i,i) = 1;
    // gamma_1 * S1
    if (n1 >= 2) {
      j = nVector.getIndex(n1-2,n2);
      M(i,j) -= gamma_1* sqrt((data_t)n1*(n1-1));
    }
    if (n2 >= 2) {
      j = nVector.getIndex(n1,n2-2);
      M(i,j) += gamma_1 * sqrt((data_t)n2*(n2-1));
    }
    if (n1 + n2 <= nVector.getNMax() - 2) {
      j = nVector.getIndex(n1+2,n2);
      M(i,j) += gamma_1 * sqrt((data_t)(n1+1)*(n1+2));
      j = nVector.getIndex(n1,n2+2);
      M(i,j) -= gamma_1 * sqrt((data_t)(n2+1)*(n2+2));
    }
    // gamma_2 * S2
    if (n1 >= 1 && n2 >= 1) {
      j = nVector.getIndex(n1-1,n2-1);
      M(i,j) -= gamma_2 * sqrt((data_t)n1*n2);
    }
    if (n1 + n2 <= nVector.getNMax() - 2) {
      j = nVector.getIndex(n1+1,n2+1);
      M(i,j) += gamma_2 * sqrt((data_t)(n1+1)*(n2+1));
    }
  }
  return M;
}

void ImageTransformation::shear(CoefficientVector<data_t>& cartesianCoeffs, Complex gamma, NumMatrix<data_t>* cov,  History* history) {
  if (history != NULL)
    (*history) << "# Applying shear gamma = " << gamma << endl;
  NumMatrix<data_t> M = getShearMatrix(gamma,cartesianCoeffs.getIndexVector());
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

void ImageTransformation::shear(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& M, NumMatrix<data_t>* cov,  History* history) {
  if (history != NULL)
    (*history) << "# Applying shear via given transformation matrix" << endl;
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

NumMatrix<data_t> ImageTransformation::getFlexionMatrix(Complex F, Complex G, const IndexVector& nVector) {
  // transform 1st and 2nd flexion into shear derivatives
  NumMatrix<data_t> dGamma(2,2);
  dGamma(0,0) = -0.5*(real(F) + real(G)); // minus signs = inverse trafo!
  dGamma(0,1) = -0.5*(imag(G) - imag(F));
  dGamma(1,0) = -0.5*(imag(F) + imag(G));
  dGamma(1,1) = -0.5*(real(F) - real(G));  
  // set up flexion matrix
  NumMatrix<data_t> M(nVector.getNCoeffs(), nVector.getNCoeffs());
  int n1, n2,j;
  for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    M(i,i) = 1;
    // dGamma(0,0) * S11
    if (n1 >= 3) {
      j = nVector.getIndex(n1-3,n2);
      M(i,j) -= 0.5 * M_SQRT1_2 * dGamma(0,0) * sqrt((data_t)n1*(n1-1)*(n1-2));
    }
    if (n1 + n2 <= nVector.getNMax() - 3) {
      j = nVector.getIndex(n1+3,n2);
      M(i,j) += dGamma(0,0) * sqrt((data_t)(n1+1)*(n1+2)*(n1+3));
    }
    if (n1 >= 1) {
      j = nVector.getIndex(n1-1,n2);
      M(i,j) -= 0.5 * M_SQRT1_2 * dGamma(0,0) * (n2-2)*sqrt((data_t)n1);
    }
    if (n1 + n2 <= nVector.getNMax() - 1) {
      j = nVector.getIndex(n1+1,n2);
      M(i,j) += 0.5 * M_SQRT1_2 * dGamma(0,0) * (n1+3)*sqrt((data_t)n1+1);
    }
    // dGamma(0,1) * S12; S12 = - {S11}| 1 <-> 2
    if (n2 >= 3) {
      j = nVector.getIndex(n1,n2-3);
      M(i,j) += 0.5 * M_SQRT1_2 * dGamma(0,1) * sqrt((data_t)n2*(n2-1)*(n2-2));
    }
    if (n1 + n2 <= nVector.getNMax() - 3) {
      j = nVector.getIndex(n1,n2+3);
      M(i,j) -= 0.5 * M_SQRT1_2 * dGamma(0,1) * sqrt((data_t)(n2+1)*(n2+2)*(n2+3));
    }
    if (n2 >= 1) {
      j = nVector.getIndex(n1,n2-1);
      M(i,j) += 0.5 * M_SQRT1_2 * dGamma(0,1) * (n1-2)*sqrt((data_t)n2);
    }
    if (n1 + n2 <= nVector.getNMax() -1) {
      j = nVector.getIndex(n1,n2+1);
      M(i,j) -= 0.5 * M_SQRT1_2 * dGamma(0,1) * (n2+3)*sqrt((data_t)n2+1);
    }
    // dGamma(1,0) * S21
    if (n2 >= 3) {
      j = nVector.getIndex(n1,n2-3);
      M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,0) * sqrt((data_t)n2*(n2-1)*(n2-2));
    }
    if (n1 + n2 <= nVector.getNMax() - 3) {
      j = nVector.getIndex(n1,n2+3);
      M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,0) * sqrt((data_t)(n2+1)*(n2+2)*(n2+3));
    }
    if (n2 >= 1) {
      j = nVector.getIndex(n1,n2-1);
      M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,0) * (2*n1+n2-3)*sqrt((data_t)n2);
      if (n1 >= 2) {
	j = nVector.getIndex(n1-2,n2-1);
	M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,0) * 3*sqrt((data_t)n1*(n1-1)*n2);
      }
      if (n1 + n2 <= nVector.getNMax() - 2) {
	j = nVector.getIndex(n1+2,n2-1);
	M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,0) * sqrt((data_t)(n1+1)*(n1+2)*n2);
      }
    }
    if (n1 + n2 <= nVector.getNMax() -1) {
      j = nVector.getIndex(n1,n2+1);
      M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,0) * (n1+n2+5)*sqrt((data_t)n2+1);
      if (n1 + n2 <= nVector.getNMax() - 3) {
	j = nVector.getIndex(n1+2,n2+1);
	M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,0) * 3*sqrt((data_t)(n1+1)*(n1+2)*(n2+1));
      }
    }
    // dGamma(1,1) * S22; S22 = + {S21} | 1 <-> 2
    if (n1 >= 3) {
      j = nVector.getIndex(n1-3,n2);
      M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,1) * sqrt((data_t)n1*(n1-1)*(n1-2));
    }
    if (n1 + n2 <= nVector.getNMax() - 3) {
      j = nVector.getIndex(n1+3,n2);
      M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,1) * sqrt((data_t)(n1+1)*(n1+2)*(n1+3));
    }
    if (n1 >= 1) {
      j = nVector.getIndex(n1-1,n2);
      M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,1) * (2*n2+n1-3)*sqrt((data_t)n1);
      if (n2 >= 2) {
	j = nVector.getIndex(n1-1,n2-2);
	M(i,j) -= 0.25 * M_SQRT1_2 * dGamma(1,1) * 3*sqrt((data_t)n2*(n2-1)*n1);
      }
      if (n1 + n2 <= nVector.getNMax() - 2) {
	j = nVector.getIndex(n1-1,n2+2);
	M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,1) * sqrt((data_t)(n2+1)*(n2+2)*n1);
      }
    }
    if (n1 + n2 <= nVector.getNMax() - 1) {
      j = nVector.getIndex(n1+1,n2);
      M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,1) * (n1+n2+5)*sqrt((data_t)n1+1);
      if (n1 + n2 <= nVector.getNMax() - 3) {
	j = nVector.getIndex(n1+1,n2+2);
	M(i,j) += 0.25 * M_SQRT1_2 * dGamma(1,1) * 3*sqrt((data_t)(n2+1)*(n2+2)*(n1+1));
      }
    }
  }
  return M;
}


void ImageTransformation::flex(CoefficientVector<data_t>& cartesianCoeffs, Complex F, Complex G, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Applying flexion to the image, F = " << F << ", G = " << G << endl;
  NumMatrix<data_t> M = getFlexionMatrix(F,G,cartesianCoeffs.getIndexVector());
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

void ImageTransformation::flex(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& M, NumMatrix<data_t>* cov, History* history) {
  if (history != NULL)
    (*history) << "# Applying flexion to the image via given transformation matrix" << endl;
  cartesianCoeffs = M * cartesianCoeffs;
  if (cov != NULL && cov->getColumns() == cartesianCoeffs.getNCoeffs())
    *cov = (M*(*cov)*M.transpose());
}

// *********** private member functions ****************** //
// computes tensor of Paper II, eq. (8)
void ImageTransformation::makeBTensor(boost::multi_array<data_t,3>& bt, data_t alpha_1, data_t beta_1, data_t gamma_1, int nmax) {
  data_t nu=1./sqrt(1./(alpha_1*alpha_1) +1./(beta_1*beta_1) +1./(gamma_1*gamma_1));
  data_t a=M_SQRT2*nu/alpha_1;
  data_t b=M_SQRT2*nu/beta_1;
  data_t c=M_SQRT2*nu/gamma_1;
  
  boost::multi_array<data_t,3> prefactor(boost::extents[nmax+1][nmax+1][nmax+1]);
  
  for (int i=0; i<= nmax; i++)
    for (int j=0; j<= nmax; j++)
      for (int k=0; k<= nmax; k++)
	prefactor[i][j][k] = nu * 1./sqrt(gsl_pow_int(2.,i+j+k-1)* M_SQRTPI * 
	  gsl_sf_fact(i)*gsl_sf_fact(j)*gsl_sf_fact(k)*alpha_1*beta_1*gamma_1);

  bt[0][0][0] = 1;
  
  for (int i=0; i<= nmax; i++) {
    for (int l=0; l<= i; l++) {
    for (int m=0; m<= i; m++) {
    for (int n=0; n<= i; n++) {
      data_t c1,c2,c3;
      if (l >= m && l >= n) {
        c1 = c2 = c3 = 0;
        if ((l-1) >= 0 && (l-2) >= 0) c1 = 2*(l-1)*(a*a -1) * bt[l-2][m][n];
        if ((l-1) >= 0 && (m-1) >= 0) c2 = 2*m*a*b* bt[l-1][m-1][n];
        if ((l-1) >= 0 && (n-1) >= 0) c3 = 2*n*a*c* bt[l-1][m][n-1];
      }
      if (m >= l && m >= n) {
        c1 = c2 = c3 = 0;
        if ((m-1) >= 0 && (m-2) >= 0) c1 = 2*(m-1)*(b*b -1) * bt[l][m-2][n];
        if ((m-1) >= 0 && (l-1) >= 0) c2 = 2*l*a*b * bt[l-1][m-1][n];
        if ((m-1) >= 0 && (n-1) >= 0) c3 = 2*n*b*c * bt[l][m-1][n-1];
      }
      if (n >= l && n >= m) {
        c1 = c2 = c3 = 0;
        if ((n-1) >= 0 && (n-2) >= 0) c1 = 2*(n-1)*(c*c - 1) * bt[l][m][n-2];
        if ((n-1) >= 0 && (l-1) >= 0) c2 = 2*l*a*c * bt[l-1][m][n-1];
        if ((n-1) >= 0 && (m-1) >= 0) c3 = 2*m*b*c * bt[l][m-1][n-1];
      }
      bt[l][m][n] = c1+c2+c3;
      if (l == 0 && m == 0 && n == 0) bt[l][m][n]=1;
    }
    }
    }
  }

  for (int l=0; l<= nmax; l++)
    for (int m=0; m<= nmax; m++)
      for (int n=0; n<= nmax; n++)
	bt[l][m][n] *= prefactor[l][m][n];
}



// this does not look alike, but is effectively identical to eq. (A3) in Paper I,
// but obtained by myself. This formulation contains less factorials and powers
// and is therefore somewhat faster for large matrices than the original formulation.
NumMatrix<data_t> ImageTransformation::make1DRescalingMatrix(data_t beta1, data_t beta2, int nmax) {
  data_t b1 = (beta1*beta1 - beta2*beta2)/(beta1*beta1 + beta2*beta2);
  data_t b2 = 2*beta1*beta2/(beta1*beta1 + beta2*beta2);
  // loop over all entries i,j
  NumMatrix<data_t> M1D(nmax+1,nmax+1);
  for (int i=0; i<= nmax; i++) {
    for (int j=0; j<= nmax; j++) {
      M1D(i,j) = 0;
      // the three sums over n,m,l
      for (int n=0; n <= (i+j)/2; n++) {
	for (int m=0; m <= n; m++) {
	  for (int l=0; l<= m; l++) {
	    // the selection function S(i,j,n,m,l)
	    if (i==(n-m+2*l) && j==(n+m-2*l)) {
	      M1D(i,j) += gsl_pow_int(-1,l) / (gsl_sf_fact(n-m)*gsl_sf_fact(m-l)*gsl_sf_fact(l)) *
		gsl_pow_int(2*b2,n-m) * gsl_pow_int(b1,m);
	    }
	  }
	}
      }
      M1D(i,j) *= sqrt(b2*gsl_sf_fact(i)*gsl_sf_fact(j)/(gsl_pow_int(2.,i)*gsl_pow_int(2.,j)));
    }
  }
  return M1D;
}
