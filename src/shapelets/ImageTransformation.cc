#include <shapelets/ImageTransformation.h>
#include <shapelets/CoefficientVector.h>
#include <NumVector.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;
typedef complex<data_t> Complex;
const Complex I = Complex(0,1);


ImageTransformation::ImageTransformation() {
}

void ImageTransformation::rotate(CoefficientVector<Complex>& polarCoeffs, data_t rho, History& history) {
  data_t rho_scaled = 2*M_PI*rho/360;
  history << "# Rotating image by " << rho << " degrees" << endl;
  const IndexVector& nVector = polarCoeffs.getIndexVector();
  int m;
  for (unsigned int i=0; i < polarCoeffs.size(); i++) {
    m = nVector.getState2(i);
    polarCoeffs(i) *= exp(m*rho_scaled*I);
  }
}

void ImageTransformation::converge(CoefficientVector<Complex>& polarCoeffs, data_t kappa, History& history) {
  history << "# Converging image by a factor 1 + kappa = " << 1 + kappa << endl;
  NumMatrix<Complex> tmp = polarCoeffs.getCoeffMatrix();
  // FIXME: which method to use: 
  // change beta, convolve with delta function, or the ladder ops
  // try eq. 39 in Paper III, only for infinitesimal transformations
  // use convolution method for bigger kappas
  const IndexVector& nVector = polarCoeffs.getIndexVector();
  int n,m;
  for (unsigned int i=0; i < polarCoeffs.size(); i++) {
    n = nVector.getState1(i);
    m = nVector.getState2(i);
    polarCoeffs(i) = (1+kappa) * tmp(n,m);
    if (n >= 2) 
      polarCoeffs(i) += kappa/2 * sqrt((data_t)(n-m)*(n+m)) * tmp(n-2,m);
    if (n < nVector.getNMax() - 1)
      polarCoeffs(i) -= kappa/2 * sqrt((data_t)(n-m+2)*(n+m+2)) * tmp(n+2,m);
  }
}

// shear with polar coefficients
// see Paper III, eq. 41
void ImageTransformation::shear(CoefficientVector<Complex>& polarCoeffs, Complex gamma, History& history) {
  history << "# Shearing image by gamma = " << gamma << endl;
  NumMatrix<Complex> tmp = polarCoeffs.getCoeffMatrix();
  Complex factor = data_t(0.25)*gamma;
  const IndexVector& nVector = polarCoeffs.getIndexVector();
  int n,m;
  for (unsigned int i=0; i < polarCoeffs.size(); i++) {
    n = nVector.getState1(i);
    m = nVector.getState2(i);
    if (n >= 2) {
      if(m >= -n+2) {
	polarCoeffs(i) += factor * sqrt((data_t)(n+m)*(n+m-2)) * tmp(n-2,m-2);
      }
      if(m <= n-2) {
	polarCoeffs(i) += conj(factor) * sqrt((data_t)(n-m)*(n-m-2)) * tmp(n-2,m+2);
      }
    }
    if (n < nVector.getNMax() - 1) { 
      if(m >= -n+2) {
	polarCoeffs(i) -= factor * sqrt((data_t)(n-m+2)*(n-m+4)) * tmp(n+2,m-2);
      }
      if(m <= n-2) {
	polarCoeffs(i) -= conj(factor) * sqrt((data_t)(n+m+2)*(n+m+4)) * tmp(n+2,m+2);
      }
    }
  }
}


// void ImageTransformation::converge(CoefficientVector<data_t>& cartesianCoeffs, data_t kappa, History& history) {
//   history << "# Applying convergence kappa = " << kappa << endl;
//   const IndexVector& nVector = cartesianCoeffs.getIndexVector();
//   NumMatrix<data_t> tmp = cartesianCoeffs.getCoeffMatrix();
//   int n1, n2;
//   data_t conv;
//   for (unsigned int i=0; i < cartesianCoeffs.size(); i++) {
//     n1 = nVector.getState1(i);
//     n2 = nVector.getState2(i);
//     conv = 0;
//     if (n1 >= 2)
//       conv -= sqrt((data_t)n1*(n1-1)) * tmp(n1-2,n2);
//     if (n2 >= 2)
//       conv -= sqrt((data_t)n2*(n2-1)) * tmp(n1,n2-2);
//     if (n1 < nVector.getNMax() - 1)
//       conv += sqrt((data_t)(n1+1)*(n1+2)) * tmp(n1+2,n2);
//     if (n2 < nVector.getNMax() - 1)
//       conv += sqrt((data_t)(n2+1)*(n2+2)) * tmp(n1,n2+2);
//     cartesianCoeffs(i) *= 1+kappa;
//     cartesianCoeffs(i) += 0.5*kappa*conv;
//   }
// }

// void ImageTransformation::shear(CoefficientVector<data_t>& cartesianCoeffs, Complex gamma, History& history) {
//   history << "# Applying shear gamma = " << gamma << endl;
//   const IndexVector& nVector = cartesianCoeffs.getIndexVector();
//   NumMatrix<data_t> tmp = cartesianCoeffs.getCoeffMatrix();
//   int n1, n2;
//   data_t sheared1, sheared2;
//   for (unsigned int i=0; i < cartesianCoeffs.size(); i++) {
//     n1 = nVector.getState1(i);
//     n2 = nVector.getState2(i);
//     sheared1 = sheared2 = 0;
//     // gamma(0) * S1
//     if (n1 >= 2)
//       sheared1 -= sqrt((data_t)n1*(n1-1)) * tmp(n1-2,n2);
//     if (n2 >= 2)
//       sheared1 += sqrt((data_t)n2*(n2-1)) * tmp(n1,n2-2);
//     if (n1 < nVector.getNMax() - 1)
//       sheared1 += sqrt((data_t)(n1+1)*(n1+2)) * tmp(n1+2,n2);
//     if (n2 < nVector.getNMax() - 1)
//       sheared1 -= sqrt((data_t)(n2+1)*(n2+2)) * tmp(n1,n2+2);
//     sheared1 *= 0.5*real(gamma);
//     // gamma(1) * S2
//     if (n1 >= 1 && n2 >= 1)
//       sheared2 += sqrt((data_t)n1*n2) * tmp(n1-1,n2-1);
//     if (n1 < nVector.getNMax() && n2 < nVector.getNMax())
//       sheared2 += sqrt((data_t)(n1+1)*(n2+1)) * tmp(n1+1,n2+1);
//     sheared2 *= imag(gamma);
//     // add to cartesianCoeff
//     cartesianCoeffs(i) += sheared1 + sheared2;
//   }
// }

void ImageTransformation::flex(CoefficientVector<data_t>& cartesianCoeffs, const NumMatrix<data_t>& dGamma, History& history) {
  history << "# Applying flexion to the image by " << dGamma << endl;
  const IndexVector& nVector = cartesianCoeffs.getIndexVector();
  NumMatrix<data_t> tmp = cartesianCoeffs.getCoeffMatrix();
  int n1, n2;
  data_t flexed;
  for (unsigned int i=0; i < cartesianCoeffs.size(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    flexed = 0;
    // dGamma(0,0) * S11
    if (n1 >= 3)
      flexed += dGamma(0,0) * 2*sqrt((data_t)n1*(n1-1)*(n1-2)) * tmp(n1-3,n2);
    if (n1 < nVector.getNMax() - 2)
      flexed -= dGamma(0,0) * 2*sqrt((data_t)(n1+1)*(n1+2)*(n1+3))* tmp(n1+3,n2);
    if (n1 >= 1)
      flexed += dGamma(0,0) * 2*(n2-2)*sqrt((data_t)n1) * tmp(n1-1,n2);
    if (n1 < nVector.getNMax())
      flexed -= dGamma(0,0) * 2*(n1+3)*sqrt((data_t)n1+1) * tmp(n1+1,n2);

    // dGamma(0,1) * S12; S12 = - {S11}| 1 <-> 2
    if (n2 >= 3)
      flexed -= dGamma(0,1) * 2*sqrt((data_t)n2*(n2-1)*(n2-2)) * tmp(n1,n2-3);
    if (n2 < nVector.getNMax() - 2)
      flexed += dGamma(0,1) * 2*sqrt((data_t)(n2+1)*(n2+2)*(n2+3))* tmp(n1,n2+3);
    if (n2 >= 1)
      flexed -= dGamma(0,1) * 2*(n1-2)*sqrt((data_t)n2) * tmp(n1,n2-1);
    if (n2 < nVector.getNMax())
      flexed += dGamma(0,1) * 2*(n2+3)*sqrt((data_t)n2+1) * tmp(n1,n2+1);

    // dGamma(1,0) * S21
    if (n2 >= 3)
      flexed += dGamma(1,0) * sqrt((data_t)n2*(n2-1)*(n2-2)) * tmp(n1,n2-3);
    if (n2 < nVector.getNMax() - 2)
      flexed -= dGamma(1,0) * sqrt((data_t)(n2+1)*(n2+2)*(n2+3)) * tmp(n1,n2+3);
    if (n2 >= 1) {
      flexed += dGamma(1,0) * (2*n1+n2-3)*sqrt((data_t)n2) * tmp(n1,n2-1);
      if (n1 >= 2)
	flexed += dGamma(1,0) * 3*sqrt((data_t)n1*(n1-1)*n2) * tmp(n1-2,n2-1);
      if (n1 < nVector.getNMax() - 1)
	flexed -= dGamma(1,0) * sqrt((data_t)(n1+1)*(n1+2)*n2) * tmp(n1+2,n2-1);
    }
    if (n2 < nVector.getNMax()) {
      flexed -= dGamma(1,0) * (n1+n2+5)*sqrt((data_t)n2+1) * tmp(n1,n2+1);
      if (n1 < nVector.getNMax() - 1) 
	flexed -= dGamma(1,0) * 3*sqrt((data_t)(n1+1)*(n1+2)*(n2+1)) * tmp(n1+2,n2+1);
    }
    // dGamma(1,1) * S22; S22 = + {S21} | 1 <-> 2
    if (n1 >= 3)
      flexed += dGamma(1,1) * sqrt((data_t)n1*(n1-1)*(n1-2)) * tmp(n1-3,n2);
    if (n1 < nVector.getNMax() - 2)
      flexed -= dGamma(1,1) * sqrt((data_t)(n1+1)*(n1+2)*(n1+3)) * tmp(n1+3,n2);
    if (n1 >= 1) {
      flexed += dGamma(1,1) * (2*n2+n1-3)*sqrt((data_t)n1) * tmp(n1-1,n2);
      if (n2 >= 2)
	flexed += dGamma(1,1) * 3*sqrt((data_t)n2*(n2-1)*n1) * tmp(n1-1,n2-2);
      if (n2 < nVector.getNMax() - 1)
	flexed -= dGamma(1,1) * sqrt((data_t)(n2+1)*(n2+2)*n1) * tmp(n1-1,n2+2);
    }
    if (n1 < nVector.getNMax()) {
      flexed -= dGamma(1,1) * (n1+n2+5)*sqrt((data_t)n1+1) * tmp(n1+1,n2);
      if (n2 < nVector.getNMax() - 1)
	flexed -= dGamma(1,1) * 3*sqrt((data_t)(n2+1)*(n2+2)*(n1+1)) * tmp(n1+1,n2+2);
    }

    // prefactor
    flexed *= -0.25*M_SQRT1_2;
    cartesianCoeffs(i) += flexed;
  }
}

void ImageTransformation::translate(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t dx1, data_t dx2, History& history) {
  history << "# Translating image by " << dx1 << "/" << dx2 << endl;
  // rescale dx1 and dx2 to be in units of beta, 
  // change sign because eq. in paper gives invers transformation
  dx1 *= -1./beta;
  dx2 *= -1./beta;
  // copy coeffs and calculate new coeffs according to eq. (32) in shapelets I.
  NumMatrix<data_t> tmp = cartesianCoeffs.getCoeffMatrix();
  const IndexVector& nVector = cartesianCoeffs.getIndexVector();
  int n1,n2;
  for (unsigned int i=0; i<cartesianCoeffs.size(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1>0)
      cartesianCoeffs(i) -= dx1*M_SQRT1_2*tmp(n1-1,n2);
    if (n2>0)
      cartesianCoeffs(i) -= dx2*M_SQRT1_2*tmp(n1,n2-1);
    if (n1< nVector.getNMax())
      cartesianCoeffs(i) += dx1*M_SQRT1_2*tmp(n1+1,n2);
    if (n2< nVector.getNMax())
      cartesianCoeffs(i) += dx2*M_SQRT1_2*tmp(n1,n2+1);
  }
}

void ImageTransformation::circularize(CoefficientVector<Complex>& polarCoeffs, History& history) {
  history << "# Circularizing image" << endl; 
  const IndexVector& nVector = polarCoeffs.getIndexVector();
  int m;
  for (unsigned int i=0; i < polarCoeffs.size(); i++) {
    m = nVector.getState2(i);
    if (m != 0)
      polarCoeffs(i) = Complex(0,0);
  }
}

void ImageTransformation::flipX(CoefficientVector<Complex>& polarCoeffs, History& history) {
  history << "# Flipping image arround its X-axis" << endl;
  const IndexVector& nVector = polarCoeffs.getIndexVector();
  for (unsigned int i=0; i < polarCoeffs.size(); i++) {
    polarCoeffs(i) = conj(polarCoeffs(i));
  }  
}

void ImageTransformation::brighten(CoefficientVector<data_t>& cartesianCoeffs, data_t factor, History& history) {
  history << "# Changing image brightness by the factor " << factor << endl;
  for (unsigned int i=0; i < cartesianCoeffs.size(); i++) {
    cartesianCoeffs(i) *= factor;
  }
}

void ImageTransformation::convolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_kernel, History& history) {
  history << "# Convolving image with kernel of order " << KernelCoeffs.getNMax() << ", beta = " << beta_kernel << endl;
  int nmax_orig = cartesianCoeffs.getNMax();

  // natural choice for convolved beta
  data_t beta_convolved = sqrt(beta*beta + beta_kernel*beta_kernel);
  int nmax_convolved = nmax_orig;

  // construct covolution matrix
  NumMatrix<data_t> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // perform the convolution
  cartesianCoeffs = P*cartesianCoeffs;
  beta = beta_convolved;
}

void ImageTransformation::deconvolve(CoefficientVector<data_t>& cartesianCoeffs, data_t& beta, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_kernel, History& history) {
  history << "# Deconvolving image with kernel of order " << KernelCoeffs.getNMax() << ", beta = " << beta_kernel << endl;
  int nmax_convolved = cartesianCoeffs.getNMax();
  data_t beta_convolved  = beta;
  // see Paper II, chapter 3.2 and 5 for details
  data_t beta_orig;
  if (beta_kernel < beta_convolved)
    beta_orig = sqrt(beta_convolved*beta_convolved - beta_kernel*beta_kernel);
  else
    beta_orig = beta_convolved;
  int nmax_orig = nmax_convolved;

  NumMatrix<data_t> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta_orig,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // since matrix could be singular, use SVD for inversion
  NumMatrix<data_t> P_1 = P.svd_invert();
  // deconvolve
  cartesianCoeffs = P_1*cartesianCoeffs;
  // go back the the smaller beta
  beta = beta_orig;
}

void ImageTransformation::makeConvolutionMatrix(NumMatrix<data_t>& P, const CoefficientVector<data_t>& KernelCoeffs, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, unsigned int nmax_orig, unsigned int nmax_convolved) {
  unsigned int nmax = GSL_MAX_INT(nmax_orig,nmax_convolved);
  unsigned int nmax_kernel = KernelCoeffs.getNMax();
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
  IndexVectorCartesian  nVector_orig(nmax_orig), nVector_convolved(nmax_convolved);
  const IndexVector& nVector_kernel = KernelCoeffs.getIndexVector();
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
  P = NumMatrix<data_t>(nCoeffs_convolved,nCoeffs_orig);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++)
      for (int l=0; l < nCoeffs_kernel; l++)
	P(n,m) += c2[n][m][l] * KernelCoeffs(l);
}

// this is a direct conversion from the shapelets IDL code
// see operations/shapelets_convolution_matrix.pro
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

void ImageTransformation::rescale(CoefficientVector<data_t>& cartesianCoeffs, data_t beta, data_t newbeta, History& history) {
  history << "# Rescaling image from beta = "<< beta << " to new beta = " << newbeta << endl;
  const IndexVector& nVector = cartesianCoeffs.getIndexVector();
  int nCoeffs = nVector.getNCoeffs();
  NumMatrix<data_t> R(nCoeffs,nCoeffs);
  makeRescalingMatrix(R,newbeta,beta,nVector);
  cartesianCoeffs = R * cartesianCoeffs;
}
 
// the 2D rescaling matrix is obtained by a tensor multiplications of
// the 1D rescaling matrix with itself.
void ImageTransformation::makeRescalingMatrix(NumMatrix<data_t>& M2D, data_t beta2, data_t beta1, const IndexVector& nVector) {
  int nmax = nVector.getNMax();
  int nCoeffs = nVector.getNCoeffs();
  // compute 1D rescaling matrix according to my own calculations
  NumMatrix<data_t>M1D(nmax+1,nmax+1);
  make1DRescalingMatrix(M1D,beta2,beta1,nmax);
  // now build tensor product of 1D matrix to give 2D rescaling matrix
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
}

// this does not look alike, but is effectively identical to eq. (A3) in Paper I,
// but obtained by myself. This formulation contains less factorials and powers
// and is therefore somewhat faster for large matrices than the original formulation.
void ImageTransformation::make1DRescalingMatrix(NumMatrix<data_t>& M1D, data_t beta2, data_t beta1, int nmax) {
  data_t b1 = (beta1*beta1 - beta2*beta2)/(beta1*beta1 + beta2*beta2);
  data_t b2 = 2*beta1*beta2/(beta1*beta1 + beta2*beta2);
  // loop over all entries i,j
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
}
