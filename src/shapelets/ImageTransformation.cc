#include <shapelets/ImageTransformation.h>
#include <shapelets/CoefficientVector.h>
#include <shapelets/MatrixManipulations.h>
#include <NumVector.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;
typedef complex<data_t> Complex;
const Complex I = Complex(0,1);

ImageTransformation::ImageTransformation() {
}

void ImageTransformation::rotate(NumMatrix<Complex>& polarCoeffs, data_t rho, ostringstream& history) {
  history << "# Rotating image by " << rho << endl;
  for (int n = 0; n < polarCoeffs.getRows(); n++) {
    for (int m = -n; m<=n; m++) {
      if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	polarCoeffs(n,mIndex(n,m)) *= exp(m*rho*I);
      }
    }
  }
}

void ImageTransformation::converge(NumMatrix<Complex>& polarCoeffs, data_t& beta, data_t kappa, ostringstream& history) {
  history << "# Converging image by a factor 1 + kappa = " << 1 + kappa << endl;
  NumMatrix<Complex> tmp = polarCoeffs;
  // FIXME: which method to use: 
  // change beta, convolve with delta function, or the ladder ops
  // try eq. 39 in Paper III, only for infinitesimal transformations
  // use convolution method for bigger kappas
  if (kappa <= 0.5) { 
    for (int n = 0; n < polarCoeffs.getRows(); n++) {
      for (int m = -n; m<=n; m++) {
	if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	  polarCoeffs(n,mIndex(n,m)) = (1+kappa)*tmp(n,mIndex(n,m));
	  if (n >= 2) 
	    polarCoeffs(n,mIndex(n,m)) += kappa/2 * sqrt((data_t)(n-m)*(n+m)) * tmp(n-2,mIndex(n-2,m));
	  if (n < polarCoeffs.getRows() -2)
	    polarCoeffs(n,mIndex(n,m)) -= kappa/2 * sqrt((data_t)(n-m+2)*(n+m+2)) * tmp(n+2,mIndex(n+2,m));
	}
      }
    }
  }
// this can work also for finite dilations
// beta += kappa*beta;
}

// see Paper III, eq. 41
void ImageTransformation::shear(NumMatrix<Complex>& polarCoeffs, data_t gamma0, data_t gamma1, ostringstream& history) {
  history << "# Shearing image by gamma0 = " << gamma0 << ", gamma1 = " << gamma1 << endl;
  NumMatrix<Complex> tmp = polarCoeffs;
  Complex factor = data_t(0.25)*(gamma0 + I*gamma1);
  // FIXME: only infinite trafos
  //if (gamma0 <= 0.5 && gamma1 <= 0.5) { 
    for (int n = 0; n < polarCoeffs.getRows(); n++) {
      for (int m = -n; m<=n; m++) {
	if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	  polarCoeffs(n,mIndex(n,m)) = tmp(n,mIndex(n,m));
	  if (n >= 2) {
	    if(m >= -n+2) {
	      polarCoeffs(n,mIndex(n,m)) += factor * sqrt((data_t)(n+m)*(n+m-2))
		* tmp(n-2,mIndex(n-2,m-2));
	    }
	    if(m <= n-2) {
	      polarCoeffs(n,mIndex(n,m)) += conj(factor) * sqrt((data_t)(n-m)*(n-m-2))
		* tmp(n-2,mIndex(n-2,m+2));
	    }
	  }
	  if (n < polarCoeffs.getRows() -2) {
	    if(m >= -n+2) {
	      polarCoeffs(n,mIndex(n,m)) -= factor * sqrt((data_t)(n-m+2)*(n-m+4))
		* tmp(n+2,mIndex(n+2,m-2));
	    }
	    if(m <= n-2) {
	      polarCoeffs(n,mIndex(n,m)) -= conj(factor) * sqrt((data_t)(n+m+2)*(n+m+4))
		* tmp(n+2,mIndex(n+2,m+2));
	    }
	  }
	}
      }
    }
    //}
}
void ImageTransformation::flex(NumMatrix<data_t>& cartesianCoeffs, const NumMatrix<data_t>& dGamma, ostringstream& history) {
 history << "# Applying flexion to the image by " << dGamma << endl;
 NumMatrix<data_t> flexed(cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
 for (int n1=0; n1 < cartesianCoeffs.getRows(); n1++) {
   for (int n2=0; n2 < cartesianCoeffs.getColumns(); n2++) {
     // dGamma(0,0) * S11
     if (n1 >= 3)
       flexed(n1,n2) += dGamma(0,0) * 2*sqrt((data_t)n1*(n1-1)*(n1-2)) * cartesianCoeffs(n1-3,n2);
     if (n1 < cartesianCoeffs.getRows() -3)
       flexed(n1,n2) -= dGamma(0,0) * 2*sqrt((data_t)(n1+1)*(n1+2)*(n1+3))* cartesianCoeffs(n1+3,n2);
     if (n1 >= 1)
       flexed(n1,n2) += dGamma(0,0) * 2*(n2-2)*sqrt((data_t)n1) * cartesianCoeffs(n1-1,n2);
     if (n1 < cartesianCoeffs.getRows() -1)
       flexed(n1,n2) -= dGamma(0,0) * 2*(n1+3)*sqrt((data_t)n1+1) * cartesianCoeffs(n1+1,n2);

     // dGamma(0,1) * S12; S12 = - {S11}| 1 <-> 2
     if (n2 >= 3)
       flexed(n1,n2) -= dGamma(0,1) * 2*sqrt((data_t)n2*(n2-1)*(n2-2)) * cartesianCoeffs(n1,n2-3);
     if (n2 < cartesianCoeffs.getColumns() -3)
       flexed(n1,n2) += dGamma(0,1) * 2*sqrt((data_t)(n2+1)*(n2+2)*(n2+3))* cartesianCoeffs(n1,n2+3);
     if (n2 >= 1)
       flexed(n1,n2) -= dGamma(0,1) * 2*(n1-2)*sqrt((data_t)n2) * cartesianCoeffs(n1,n2-1);
     if (n2 < cartesianCoeffs.getColumns() -1)
       flexed(n1,n2) += dGamma(0,1) * 2*(n2+3)*sqrt((data_t)n2+1) * cartesianCoeffs(n1,n2+1);

     // dGamma(1,0) * S21
     if (n2 >= 3)
       flexed(n1,n2) += dGamma(1,0) * sqrt((data_t)n2*(n2-1)*(n2-2)) * cartesianCoeffs(n1,n2-3);
     if (n2 < cartesianCoeffs.getColumns() -3)
       flexed(n1,n2) -= dGamma(1,0) * sqrt((data_t)(n2+1)*(n2+2)*(n2+3)) * cartesianCoeffs(n1,n2+3);
     if (n2 >= 1) {
       flexed(n1,n2) += dGamma(1,0) * (2*n1+n2-3)*sqrt((data_t)n2) * cartesianCoeffs(n1,n2-1);
       if (n1 >= 2)
	 flexed(n1,n2) += dGamma(1,0) * 3*sqrt((data_t)n1*(n1-1)*n2) * cartesianCoeffs(n1-2,n2-1);
       if (n1 < cartesianCoeffs.getRows() -2)
	 flexed(n1,n2) -= dGamma(1,0) * sqrt((data_t)(n1+1)*(n1+2)*n2) * cartesianCoeffs(n1+2,n2-1);
     }
     if (n2 < cartesianCoeffs.getColumns() -1) {
       flexed(n1,n2) -= dGamma(1,0) * (n1+n2+5)*sqrt((data_t)n2+1) * cartesianCoeffs(n1,n2+1);
       if (n1 < cartesianCoeffs.getRows() -2) 
	 flexed(n1,n2) -= dGamma(1,0) * 3*sqrt((data_t)(n1+1)*(n1+2)*(n2+1)) * cartesianCoeffs(n1+2,n2+1);
     }
     // dGamma(1,1) * S22; S22 = + {S21} | 1 <-> 2
     if (n1 >= 3)
       flexed(n1,n2) += dGamma(1,1) * sqrt((data_t)n1*(n1-1)*(n1-2)) * cartesianCoeffs(n1-3,n2);
     if (n1 < cartesianCoeffs.getRows() -3)
       flexed(n1,n2) -= dGamma(1,1) * sqrt((data_t)(n1+1)*(n1+2)*(n1+3)) * cartesianCoeffs(n1+3,n2);
     if (n1 >= 1) {
       flexed(n1,n2) += dGamma(1,1) * (2*n2+n1-3)*sqrt((data_t)n1) * cartesianCoeffs(n1-1,n2);
       if (n2 >= 2)
	 flexed(n1,n2) += dGamma(1,1) * 3*sqrt((data_t)n2*(n2-1)*n1) * cartesianCoeffs(n1-1,n2-2);
       if (n2 < cartesianCoeffs.getColumns() -2)
	 flexed(n1,n2) -= dGamma(1,1) * sqrt((data_t)(n2+1)*(n2+2)*n1) * cartesianCoeffs(n1-1,n2+2);
     }
     if (n1 < cartesianCoeffs.getRows() -1) {
       flexed(n1,n2) -= dGamma(1,1) * (n1+n2+5)*sqrt((data_t)n1+1) * cartesianCoeffs(n1+1,n2);
       if (n2 < cartesianCoeffs.getColumns() -2)
	 flexed(n1,n2) -= dGamma(1,1) * 3*sqrt((data_t)(n2+1)*(n2+2)*(n1+1)) * cartesianCoeffs(n1+1,n2+2);
     }

     // prefactor
     flexed(n1,n2) *= -0.25*M_SQRT1_2;
   }
 }
 // add to coeffs
 for (int n1=0; n1 < cartesianCoeffs.getRows(); n1++)
   for (int n2=0; n2 < cartesianCoeffs.getColumns(); n2++)
     cartesianCoeffs(n1,n2) += flexed(n1,n2);
}

void ImageTransformation::translate(NumMatrix<data_t>& cartesianCoeffs, data_t beta, data_t dx1, data_t dx2, ostringstream& history) {
  history << "# Translating image by " << dx1 << "/" << dx2 << endl;
  // rescale dx1 and dx2 to be in units of beta, 
  // change sign because eq. in paper gives invers transformation
  dx1 *= -1./beta;
  dx2 *= -1./beta;
  // copy coeffs and calculate new coeffs according to eq. (32) in shapelets I.
  NumMatrix<data_t> tmp = cartesianCoeffs;
  for (int n1=0; n1 < cartesianCoeffs.getRows(); n1++) {
    for (int n2=0; n2 < cartesianCoeffs.getColumns(); n2++) {
      if (n1>0)
	cartesianCoeffs(n1,n2) -= dx1*M_SQRT1_2*tmp(n1-1,n2);
      if (n2>0)
	cartesianCoeffs(n1,n2) -= dx2*M_SQRT1_2*tmp(n1,n2-1);
      if (n1< cartesianCoeffs.getRows()-1)
	cartesianCoeffs(n1,n2) += dx1*M_SQRT1_2*tmp(n1+1,n2);
      if (n2< cartesianCoeffs.getColumns()-1)
	cartesianCoeffs(n1,n2) += dx2*M_SQRT1_2*tmp(n1,n2+1);
    }
  }
}

void ImageTransformation::circularize(NumMatrix<Complex>& polarCoeffs, ostringstream& history) {
  history << "# Circularizing image" << endl; 
  NumMatrix<Complex> tmp = polarCoeffs;
  for (int n = 0; n < polarCoeffs.getRows(); n++) 
    for (int m = -n; m<=n; m++) 
      if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) 
	if (m != 0)
	  polarCoeffs(n,mIndex(n,m)) = Complex(0,0);
}

void ImageTransformation::flipX(NumMatrix<Complex>& polarCoeffs, ostringstream& history) {
  history << "# Flipping image arround its X-axis" << endl;
  for (int n = 0; n < polarCoeffs.getRows(); n++)
    for (int m = -n; m<=n; m++) 
      if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1))
	polarCoeffs(n,mIndex(n,m)) = conj(polarCoeffs(n,mIndex(n,m)));
}  

void ImageTransformation::brighten(NumMatrix<data_t>& cartesianCoeffs, NumMatrix<Complex>& polarCoeffs, data_t factor, ostringstream& history) {
  history << "# Changing image brightness by the factor " << factor << endl;
  // cartesian coeffs first, then the same for polar coeffs
  for (int i =0; i < cartesianCoeffs.getRows();i++)
    for (int j=0; j<cartesianCoeffs.getColumns();j++)
      cartesianCoeffs(i,j) *= factor;

  for (int n = 0; n < polarCoeffs.getRows(); n++) {
    for (int m = -n; m<=n; m++) { 
      if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	polarCoeffs(n,mIndex(n,m)) *= factor;
      }
    }
  }
}

void ImageTransformation::convolve(NumMatrix<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel, ostringstream& history) {
  history << "# Convolving image with kernel of order " << KernelCoeffs.getRows() - 1;
  history << ", beta = " << beta_kernel << endl;
  int nmax_orig = cartesianCoeffs.getRows() -1;

  // default values for convolving
  data_t beta_convolved = sqrt(beta*beta + beta_kernel*beta_kernel);
  int nmax_convolved = nmax_orig;

  // construct covolution matrix
  NumMatrix<data_t> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // transform coeffs into vector form
  CoefficientVector<data_t> f(cartesianCoeffs);
  // perform the convolution
  CoefficientVector<data_t> h = P*(NumVector<data_t>)f;
  // transform vector back to coefficient matrix
  h.fillCoeffMatrix(cartesianCoeffs);
  beta = beta_convolved;
}

void ImageTransformation::deconvolve(NumMatrix<data_t>& cartesianCoeffs, data_t& beta, const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel, ostringstream& history) {
  history << "# Deconvolving image with kernel of order " << KernelCoeffs.getRows() - 1;
  history << ", beta = " << beta_kernel << endl;
  int nmax_convolved = cartesianCoeffs.getRows() -1;
  data_t beta_convolved  = beta;
  // see Paper II, chapter 3.2 and 5 for details
  data_t beta_orig = sqrt(beta_convolved*beta_convolved - beta_kernel*beta_kernel);
  int nmax_orig = nmax_convolved;
  NumMatrix<data_t> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta_orig,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // since matrix could be singular, use SVD for inversion
  NumMatrix<data_t> P_1 = P.svd_invert();

  // transform cartesianCoeffs into vector
  CoefficientVector<data_t> f(cartesianCoeffs);
  // deconvolve
  CoefficientVector<data_t> h = P_1*(NumVector<data_t>)f;
  // transform back to matrix form
  h.fillCoeffMatrix(cartesianCoeffs);

  // in contrast to the original implementation: let's go back the the smaller beta here
  beta = beta_orig;
}

void ImageTransformation::makeConvolutionMatrix(NumMatrix<data_t>& P, const NumMatrix<data_t>& KernelCoeffs, data_t beta_orig, data_t beta_kernel, data_t beta_convolved, int nmax_orig, int nmax_convolved) {
  int nmax = GSL_MAX_INT(nmax_orig,nmax_convolved);
  int nmax_kernel = KernelCoeffs.getRows() - 1;
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
  IndexVector  nVector_orig(nmax_orig), nVector_kernel(nmax_kernel), nVector_convolved(nmax_convolved);
  int nCoeffs_orig = nVector_orig.getNCoeffs();
  int nCoeffs_kernel = nVector_kernel.getNCoeffs();
  int nCoeffs_convolved = nVector_convolved.getNCoeffs();
  boost::multi_array<data_t,3> c2(boost::extents[nCoeffs_convolved][nCoeffs_orig][nCoeffs_kernel]);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++) 
      for (int l=0; l < nCoeffs_kernel; l++)
	c2[n][m][l] = 
	  c1[nVector_convolved.getN1(n)][nVector_orig.getN1(m)][nVector_kernel.getN1(l)] * 
	  c1[nVector_convolved.getN2(n)][nVector_orig.getN2(m)][nVector_kernel.getN2(l)];
  
  // now the 2D convolution matrix P_nm
  // since (n1,n2) are stored as vector, this can be regarded as matrix
  P = NumMatrix<data_t>(nCoeffs_convolved,nCoeffs_orig);
  // vectorize PSFCoeffs for multiplication
  CoefficientVector<data_t> g(KernelCoeffs);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++)
      for (int l=0; l < nCoeffs_kernel; l++)
	P(n,m) += c2[n][m][l] * g(l);
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

void ImageTransformation::rescale(NumMatrix<data_t>& cartesianCoeffs, data_t beta, data_t newbeta, ostringstream& history) {
  history << "# Rescaling image from beta = "<< beta << " to new beta = " << newbeta << endl;
  CoefficientVector<data_t> coeffVector(cartesianCoeffs);
  const IndexVector& nVector = coeffVector.getIndexVector();
  int nCoeffs = nVector.getNCoeffs();
  NumMatrix<data_t> R(nCoeffs,nCoeffs);
  makeRescalingMatrix(R,newbeta,beta,nVector);
  CoefficientVector<data_t> coeffVector_ = R * (NumVector<data_t>)coeffVector;
  coeffVector_.fillCoeffMatrix(cartesianCoeffs);

}
 
// the 2D rescaling matrix is obtained by a tensor multiplications of
// the 1D rescaling matrix with itself.
void ImageTransformation::makeRescalingMatrix(NumMatrix<data_t>& M2D, data_t beta2, data_t beta1, const IndexVector& nVector) {
  int nmax = nVector.getOrder();
  int nCoeffs = nVector.getNCoeffs();
  // compute 1D rescaling matrix according to my own calculations
  NumMatrix<data_t>M1D(nmax+1,nmax+1);
  make1DRescalingMatrix(M1D,beta2,beta1,nmax);
  // now build tensor product of 1D matrix to give 2D rescaling matrix
  int i1,i2,j1,j2;
  for (int l = 0; l < nCoeffs; l++) {
    i1 = nVector.getN1(l);
    i2 = nVector.getN2(l);
    // here we assume that both sets of coefficients have same length
    for (int i=0; i< nCoeffs; i++) {
       j1 = nVector.getN1(i);
       j2 = nVector.getN2(i);
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
