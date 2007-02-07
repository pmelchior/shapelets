#include <ImageTransformation.h>
#include <MatrixManipulations.h>
#include <NumVector.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;
typedef complex<double> Complex;
const Complex I = Complex(0,1);

ImageTransformation::ImageTransformation() {
}

void ImageTransformation::rotate(NumMatrix<Complex>& polarCoeffs, double rho, ostringstream& history) {
  history << "# Rotating image by " << rho << endl;
  for (int n = 0; n < polarCoeffs.getRows(); n++) {
    for (int m = -n; m<=n; m++) {
      if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	polarCoeffs(n,mIndex(n,m)) *= exp(m*rho*I);
      }
    }
  }
}

void ImageTransformation::converge(NumMatrix<Complex>& polarCoeffs, double& beta, double kappa, ostringstream& history) {
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
	    polarCoeffs(n,mIndex(n,m)) += kappa/2 * sqrt((double)(n-m)*(n+m)) * tmp(n-2,mIndex(n-2,m));
	  if (n < polarCoeffs.getRows() -2)
	    polarCoeffs(n,mIndex(n,m)) -= kappa/2 * sqrt((double)(n-m+2)*(n+m+2)) * tmp(n+2,mIndex(n+2,m));
	}
      }
    }
  }
// this can work also for finite dilations
// beta += kappa*beta;
}

// see Paper III, eq. 41
void ImageTransformation::shear(NumMatrix<Complex>& polarCoeffs, double gamma0, double gamma1, ostringstream& history) {
  history << "# Shearing image by gamma0 = " << gamma0 << ", gamma1 = " << gamma1 << endl;
  NumMatrix<Complex> tmp = polarCoeffs;
  Complex factor = 0.25*(gamma0 + I*gamma1);
  // FIXME: only infinite trafos
  //if (gamma0 <= 0.5 && gamma1 <= 0.5) { 
    for (int n = 0; n < polarCoeffs.getRows(); n++) {
      for (int m = -n; m<=n; m++) {
	if ((n%2==0 && m%2==0) || (n%2==1 && abs(m%2)==1)) {
	  polarCoeffs(n,mIndex(n,m)) = tmp(n,mIndex(n,m));
	  if (n >= 2) {
	    if(m >= -n+2) {
	      polarCoeffs(n,mIndex(n,m)) += factor * sqrt((double)(n+m)*(n+m-2))
		* tmp(n-2,mIndex(n-2,m-2));
	    }
	    if(m <= n-2) {
	      polarCoeffs(n,mIndex(n,m)) += conj(factor) * sqrt((double)(n-m)*(n-m-2))
		* tmp(n-2,mIndex(n-2,m+2));
	    }
	  }
	  if (n < polarCoeffs.getRows() -2) {
	    if(m >= -n+2) {
	      polarCoeffs(n,mIndex(n,m)) -= factor * sqrt((double)(n-m+2)*(n-m+4))
		* tmp(n+2,mIndex(n+2,m-2));
	    }
	    if(m <= n-2) {
	      polarCoeffs(n,mIndex(n,m)) -= conj(factor) * sqrt((double)(n+m+2)*(n+m+4))
		* tmp(n+2,mIndex(n+2,m+2));
	    }
	  }
	}
      }
    }
    //}
}
void ImageTransformation::flex(NumMatrix<double>& cartesianCoeffs, NumMatrix<double>& dGamma, ostringstream& history) {
 history << "# Applying flexion to the image by " << dGamma << endl;
 NumMatrix<double> flexed(cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
 for (int n1=0; n1 < cartesianCoeffs.getRows(); n1++) {
   for (int n2=0; n2 < cartesianCoeffs.getColumns(); n2++) {
     // dGamma(0,0) * S11
     if (n1 >= 3)
       flexed(n1,n2) += dGamma(0,0) * 2*sqrt((double)n1*(n1-1)*(n1-2)) * cartesianCoeffs(n1-3,n2);
     if (n1 < cartesianCoeffs.getRows() -3)
       flexed(n1,n2) -= dGamma(0,0) * 2*sqrt((double)(n1+1)*(n1+2)*(n1+3))* cartesianCoeffs(n1+3,n2);
     if (n1 >= 1)
       flexed(n1,n2) += dGamma(0,0) * 2*(n2-2)*sqrt((double)n1) * cartesianCoeffs(n1-1,n2);
     if (n1 < cartesianCoeffs.getRows() -1)
       flexed(n1,n2) -= dGamma(0,0) * 2*(n1+3)*sqrt((double)n1+1) * cartesianCoeffs(n1+1,n2);

     // dGamma(0,1) * S12; S12 = - {S11}| 1 <-> 2
     if (n2 >= 3)
       flexed(n1,n2) -= dGamma(0,1) * 2*sqrt((double)n2*(n2-1)*(n2-2)) * cartesianCoeffs(n1,n2-3);
     if (n2 < cartesianCoeffs.getColumns() -3)
       flexed(n1,n2) += dGamma(0,1) * 2*sqrt((double)(n2+1)*(n2+2)*(n2+3))* cartesianCoeffs(n1,n2+3);
     if (n2 >= 1)
       flexed(n1,n2) -= dGamma(0,1) * 2*(n1-2)*sqrt((double)n2) * cartesianCoeffs(n1,n2-1);
     if (n2 < cartesianCoeffs.getColumns() -1)
       flexed(n1,n2) += dGamma(0,1) * 2*(n2+3)*sqrt((double)n2+1) * cartesianCoeffs(n1,n2+1);

     // dGamma(1,0) * S21
     if (n2 >= 3)
       flexed(n1,n2) += dGamma(1,0) * sqrt((double)n2*(n2-1)*(n2-2)) * cartesianCoeffs(n1,n2-3);
     if (n2 < cartesianCoeffs.getColumns() -3)
       flexed(n1,n2) -= dGamma(1,0) * sqrt((double)(n2+1)*(n2+2)*(n2+3)) * cartesianCoeffs(n1,n2+3);
     if (n2 >= 1) {
       flexed(n1,n2) += dGamma(1,0) * (2*n1+n2-3)*sqrt((double)n2) * cartesianCoeffs(n1,n2-1);
       if (n1 >= 2)
	 flexed(n1,n2) += dGamma(1,0) * 3*sqrt((double)n1*(n1-1)*n2) * cartesianCoeffs(n1-2,n2-1);
       if (n1 < cartesianCoeffs.getRows() -2)
	 flexed(n1,n2) -= dGamma(1,0) * sqrt((double)(n1+1)*(n1+2)*n2) * cartesianCoeffs(n1+2,n2-1);
     }
     if (n2 < cartesianCoeffs.getColumns() -1) {
       flexed(n1,n2) -= dGamma(1,0) * (n1+n2+5)*sqrt((double)n2+1) * cartesianCoeffs(n1,n2+1);
       if (n1 < cartesianCoeffs.getRows() -2) 
	 flexed(n1,n2) -= dGamma(1,0) * 3*sqrt((double)(n1+1)*(n1+2)*(n2+1)) * cartesianCoeffs(n1+2,n2+1);
     }
     // dGamma(1,1) * S22; S22 = + {S21} | 1 <-> 2
     if (n1 >= 3)
       flexed(n1,n2) += dGamma(1,1) * sqrt((double)n1*(n1-1)*(n1-2)) * cartesianCoeffs(n1-3,n2);
     if (n1 < cartesianCoeffs.getRows() -3)
       flexed(n1,n2) -= dGamma(1,1) * sqrt((double)(n1+1)*(n1+2)*(n1+3)) * cartesianCoeffs(n1+3,n2);
     if (n1 >= 1) {
       flexed(n1,n2) += dGamma(1,1) * (2*n2+n1-3)*sqrt((double)n1) * cartesianCoeffs(n1-1,n2);
       if (n2 >= 2)
	 flexed(n1,n2) += dGamma(1,1) * 3*sqrt((double)n2*(n2-1)*n1) * cartesianCoeffs(n1-1,n2-2);
       if (n2 < cartesianCoeffs.getColumns() -2)
	 flexed(n1,n2) -= dGamma(1,1) * sqrt((double)(n2+1)*(n2+2)*n1) * cartesianCoeffs(n1-1,n2+2);
     }
     if (n1 < cartesianCoeffs.getRows() -1) {
       flexed(n1,n2) -= dGamma(1,1) * (n1+n2+5)*sqrt((double)n1+1) * cartesianCoeffs(n1+1,n2);
       if (n2 < cartesianCoeffs.getColumns() -2)
	 flexed(n1,n2) -= dGamma(1,1) * 3*sqrt((double)(n2+1)*(n2+2)*(n1+1)) * cartesianCoeffs(n1+1,n2+2);
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

void ImageTransformation::translate(Point2D& xcentroid, double dx0, double dx1, ostringstream& history) {
  history << "# Translating image by " << dx0 << "/" << dx1 << endl;
  xcentroid(0) += dx0;
  xcentroid(1) += dx1;
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

void ImageTransformation::brighten(NumMatrix<double>& cartesianCoeffs, NumMatrix<Complex>& polarCoeffs, double factor, ostringstream& history) {
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

void ImageTransformation::convolve(NumMatrix<double>& cartesianCoeffs, double& beta, NumMatrix<double>& KernelCoeffs, double beta_kernel, ostringstream& history) {
  history << "# Convolving image with kernel of order " << KernelCoeffs.getRows() - 1;
  history << ", beta = " << beta_kernel << endl;
  int nmax_orig = cartesianCoeffs.getRows() -1;
  // default values for convolving
  double beta_convolved = sqrt(beta*beta + beta_kernel*beta_kernel);
  int nmax_convolved = nmax_orig;

  NumMatrix<double> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // transform cartesianCoeffs into vector
  int nCoeffs_orig = getNCoeffs(nmax_orig);
  NumMatrix<int> nVector_orig;
  makeNVector(nVector_orig,nCoeffs_orig,nmax_orig);
  // this does only work because nmax_convolved = nmax_orig
  // otherwise generate nVector and nCoeffs for _convolved too;
  NumVector<double> f(nCoeffs_orig), h(nCoeffs_orig);
  matrixMapping(cartesianCoeffs,f,0,nVector_orig,nCoeffs_orig);
  h = P*f;
  vectorMapping(h,cartesianCoeffs,nVector_orig,nCoeffs_orig);
  beta = beta_convolved;
}

void ImageTransformation::deconvolve(NumMatrix<double>& cartesianCoeffs, double& beta, NumMatrix<double>& KernelCoeffs, double beta_kernel, ostringstream& history) {
  history << "# Deconvolving image with kernel of order " << KernelCoeffs.getRows() - 1;
  history << ", beta = " << beta_kernel << endl;
  int nmax_convolved = cartesianCoeffs.getRows() -1;
  double beta_convolved  = beta;
  // see Paper II, chapter 3.2 and 5 for details
  double beta_orig = sqrt(beta_convolved*beta_convolved - beta_kernel*beta_kernel);
  int nmax_orig = nmax_convolved;
  NumMatrix<double> P;
  makeConvolutionMatrix(P,KernelCoeffs,beta_orig,beta_kernel,beta_convolved,nmax_orig,nmax_convolved);
  // since matrix could be singular, use SVD for inversion
  NumMatrix<double> P_1 = P.svd_invert();
  //P.printWithIndices(1);

  // transform cartesianCoeffs into vector
  int nCoeffs_convolved = getNCoeffs(nmax_convolved);
  NumMatrix<int> nVector_convolved;
  makeNVector(nVector_convolved,nCoeffs_convolved,nmax_convolved);
  // this does only work because nmax_convolved = nmax_orig
  // otherwise generate nVector and nCoeffs for _orig too;
  NumVector<double> f(nCoeffs_convolved), h(nCoeffs_convolved);
  matrixMapping(cartesianCoeffs,f,0,nVector_convolved,nCoeffs_convolved);
  h = P_1*f;
  vectorMapping(h,cartesianCoeffs,nVector_convolved,nCoeffs_convolved);
  // in contrast to the original implementation: let's go back the the smaller beta here
  beta = beta_orig;
}

void ImageTransformation::rescale(NumMatrix<double>& cartesianCoeffs, double beta, double newbeta, ostringstream& history) {
  history << "# Rescaling image from beta = "<< beta << " to new beta = " << newbeta << endl;
  int nmax = cartesianCoeffs.getRows() -1;
  int nCoeffs = getNCoeffs(nmax);
  NumMatrix<int> nVector;
  makeNVector(nVector,nCoeffs,nmax);
  NumVector<double> coeffVector(nCoeffs);
  matrixMapping(cartesianCoeffs,coeffVector,0,nVector,nCoeffs);
  NumMatrix<double> R(nCoeffs,nCoeffs);
  makeRescalingMatrix(R,beta,newbeta,nCoeffs,nVector);
  NumVector<double> coeffVector_;
  coeffVector_ = R * coeffVector;
  vectorMapping(coeffVector_,cartesianCoeffs,nVector,nCoeffs);
  
}
 
void ImageTransformation::makeConvolutionMatrix(NumMatrix<double>& P, NumMatrix<double>& KernelCoeffs, double beta_orig, double beta_kernel, double beta_convolved, int nmax_orig, int nmax_convolved) {
  int nmax = GSL_MAX_INT(nmax_orig,nmax_convolved);
  int nmax_kernel = KernelCoeffs.getRows() - 1;
  nmax = GSL_MAX_INT(nmax_kernel,nmax);
  boost::multi_array<double,3> bt(boost::extents[nmax+1][nmax+1][nmax+1]);
  double alpha = beta_orig;
  double beta = beta_kernel;
  double gamma = beta_convolved;
  makeBTensor(bt,1./gamma,1./alpha,1./beta,nmax);
 
  // the 1D convolution tensor C_nml
  boost::multi_array<double,3> c1(boost::extents[nmax_convolved+1][nmax_orig+1][nmax_kernel+1]);
  for (int n=0; n <= nmax_convolved; n++)
    for (int m=0; m <= nmax_orig; m++)
      for (int l=0; l <= nmax_kernel; l++)
	if ((n+m+l)%2 == 0)
	  c1[n][m][l] = M_SQRT2 * M_SQRTPI * gsl_pow_int(-1,(3*n + m + l)/2) * bt[n][m][l];
  
  // the 2D convolution tensor C_nml
  // the matrix->vector projection
  int nCoeffs_orig = getNCoeffs(nmax_orig);
  int nCoeffs_kernel = getNCoeffs(nmax_kernel);
  int nCoeffs_convolved = getNCoeffs(nmax_convolved);
  NumMatrix<int> nVector_orig, nVector_kernel, nVector_convolved;
  makeNVector(nVector_orig,nCoeffs_orig,nmax_orig);
  makeNVector(nVector_kernel,nCoeffs_kernel,nmax_kernel);
  makeNVector(nVector_convolved,nCoeffs_convolved,nmax_convolved);
  boost::multi_array<double,3> c2(boost::extents[nCoeffs_convolved][nCoeffs_orig][nCoeffs_kernel]);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++) 
      for (int l=0; l < nCoeffs_kernel; l++)
	c2[n][m][l] = 
	  c1[getN1(nVector_convolved,n)][getN1(nVector_orig,m)][getN1(nVector_kernel,l)] * 
	  c1[getN2(nVector_convolved,n)][getN2(nVector_orig,m)][getN2(nVector_kernel,l)];
  
  // now the 2D convolution matrix P_nm
  // since (n1,n2) are stored as vector, this can be regarded as matrix
  P = NumMatrix<double>(nCoeffs_convolved,nCoeffs_orig);
  // vectorize PSFCoeffs for multiplication
  NumVector<double> g(nCoeffs_kernel); 
  matrixMapping(KernelCoeffs,g,0,nVector_kernel,nCoeffs_kernel);
  for (int n=0; n < nCoeffs_convolved; n++)
    for (int m=0; m < nCoeffs_orig; m++)
      for (int l=0; l < nCoeffs_kernel; l++)
	P(n,m) += c2[n][m][l] * g(l);
}

// this is a direct conversion from the shapelets IDL code
// see operations/shapelets_convolution_matrix.pro
void ImageTransformation::makeBTensor(boost::multi_array<double,3>& bt, double alpha_1, double beta_1, double gamma_1, int nmax) {
  double nu=1./sqrt(1./(alpha_1*alpha_1) +1./(beta_1*beta_1) +1./(gamma_1*gamma_1));
  double a=M_SQRT2*nu/alpha_1;
  double b=M_SQRT2*nu/beta_1;
  double c=M_SQRT2*nu/gamma_1;
  
  boost::multi_array<double,3> prefactor(boost::extents[nmax+1][nmax+1][nmax+1]);
  
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
      double c1,c2,c3;
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

// see Paper I, appendix A
// modified for 2D shapelets with equal number of coefficients
// beta1 is new beta, beta2 is beta for which coeffs were derived
void ImageTransformation::makeRescalingMatrix(NumMatrix<double>& betaTrafo, double beta1, double beta2, int nCoeffs, NumMatrix<int>& nVector) {
  int n1,n2,n1_,n2_;
  double b1 = (beta1*beta1 - beta2*beta2)/(beta1*beta1 + beta2*beta2);
  double b2 = 2*beta1*beta2/(beta1*beta1 + beta2*beta2);
  for (int coeff=0; coeff < nCoeffs; coeff++) {
    n1 = getN1(nVector,coeff);
    n2 = getN2(nVector,coeff);
    for (int coeff_=0; coeff_ < nCoeffs; coeff_++) {
      n1_ = getN1(nVector,coeff_);
      n2_ = getN2(nVector,coeff_);
      // now we know position in coeff matrix
      // we can then get trafo of the coeff 
      double entry = 0;
      for(int l=0; l <= GSL_MIN_INT(n1,n1_); l++) {
	for (int m=0; m <= GSL_MIN_INT(n2,n2_); m++) {
	  // the Pi function
	  if (((n1%2 == 0 && n1_%2 ==0 && l%2 ==0) || (n1%2 == 1 && n1_%2 == 1 && l%2 == 1)) && (((n2%2 == 0 && n2_%2 ==0 && m%2 ==0) || (n2%2 == 1 && n2_%2 == 1 && m%2 == 1)))) 
	    entry += gsl_pow_int(-1,(n1_-l)/2) * gsl_pow_int(-1,(n2_-m)/2) *
	      sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n1_)*gsl_sf_fact(n2)*gsl_sf_fact(n2_)) /
	      (gsl_sf_fact((n1-l)/2)*gsl_sf_fact((n1_-l)/2)*gsl_sf_fact(l) *
	       gsl_sf_fact((n2-m)/2)*gsl_sf_fact((n2_-m)/2)*gsl_sf_fact(m)) *
	      gsl_pow_int(0.5*b1,(n1+n1_)/2 - l + (n2+n2_)/2 - m) *
	      gsl_pow_int(b2,l+m+1);
	}
      }
      betaTrafo(coeff,coeff_) = entry;
    }
  }
}
