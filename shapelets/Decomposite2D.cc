#include <Decomposite2D.h>
#include <MatrixManipulations.h>
#include <IO.h>
#include <gsl/gsl_math.h>

namespace ublas = boost::numeric::ublas;

Decomposite2D::Decomposite2D(int innmax, double inbeta, const Object& obj) : 
data(obj.getData()), bg_rms(obj.getBackgroundRMS()) {
  npixels = data.size();
  // centroid and beta are zero order estimators coming from object extraction
  xcentroid = obj.getCentroid();
  grid = obj.getGrid();
  nmax = innmax;
  setBeta(inbeta);
  nCoeffs = getNCoeffs(nmax);
  makeNVector(nVector,nCoeffs,nmax);
  Mt = NumMatrix<double>(nCoeffs,npixels);
  background_variance = gsl_pow_2(obj.getNoiseRMS());
  if (obj.getNoiseModel().compare("GAUSSIAN")==0)
    gaussian = 1;
  else {
    gaussian = 0;
    makeV_Matrix();
  }
  change = updateC = updateModel = updateResiduals = 1;
}

// see Paper III, eq. 83 and following explanation
void Decomposite2D::makeLSMatrix () {
  NumMatrix<double> M0(nmax+1,npixels), M1(nmax+1,npixels);
  // start with 0th and 1st order shapelets
  double x0_scaled, x1_scaled;
  double factor0 = 1./sqrt(M_SQRTPI*beta);
  for (int i=0; i<npixels; i++) {
    x0_scaled = (grid(i,0) - xcentroid(0))/beta;
    x1_scaled = (grid(i,1) - xcentroid(1))/beta;
    M0(0,i) = factor0*exp(-x0_scaled*x0_scaled/2);
    M1(0,i) = factor0*exp(-x1_scaled*x1_scaled/2);
    if (nmax > 0) {
      M0(1,i) = M_SQRT2*x0_scaled*M0(0,i);
      M1(1,i) = M_SQRT2*x1_scaled*M1(0,i);
    }
  }
  // use recurrance relation to compute higher orders
  for (int n=2;n<=nmax;n++) {
    double factor1 = sqrt(1./(2*n)), factor2 =sqrt((n-1.)/n); 
    for (int i=0; i<npixels; i++) {
      M0(n,i) = 2*(grid(i,0) - xcentroid(0))/beta*factor1*M0(n-1,i) 
	- factor2*M0(n-2,i);
      M1(n,i) = 2*(grid(i,1) - xcentroid(1))/beta*factor1*M1(n-1,i) 
	- factor2*M1(n-2,i);
    }
  }
  // now build tensor product of M0 and M1
  Mt = NumMatrix<double> (nCoeffs,npixels);
  int n0, n1;
  for (int l = 0; l < nCoeffs; l++) {
    n0 = getN1(nVector,l);
    n1 = getN2(nVector,l);
    for (int i=0; i<npixels; i++) 
      Mt(l,i) = M0(n0,i)*M1(n1,i);
  }
  // as long as the noise in each pixel is constant, the covariance matrix is proportional
  // to the identity matrix, so we don't have to consider it here, but when computing chi2.
  // solution for minimizing chi^2:
  M = Mt.transpose();
  if (!gaussian) 
    Mt = Mt*V_;
  LS = (Mt*M);
  LS = LS.invert();
  LS = LS*Mt;
  
  // SVD method 
  //LS = M.svd_invert(); 
  // overlapping integrals
  //LS = M.transpose();
}

void Decomposite2D::makeV_Matrix() {
  NumVector<double> tmp = data;
  V_.resize(npixels);
  for (int i=0; i < npixels; i++)
    tmp(i) += background_variance;
  convolveGaussian(tmp,V_,grid.getSize(0),grid.getSize(1));
  for (int i=0; i< npixels; i++)
    V_(i) = 1./V_(i);
}

void Decomposite2D::computeCoeffs() {
  makeLSMatrix();
  // this is useful only for the regularization in OptimalDecomposite2D
  if (updateC)
    coeffVector = LS * data;
  change = 0;
}

void Decomposite2D::updateModelResiduals() {
  updateModel = updateResiduals = 1;
}

void Decomposite2D::computeModel() {
  if (updateModel) {
    if (change)
      computeCoeffs();
    model = M * coeffVector;
    updateModel = 0;
  }
}

void Decomposite2D::computeResiduals() {
  if(updateResiduals || updateModel) {
    if (updateModel)
      computeModel();
    residual = data;
    residual -= model;
    updateResiduals = 0;
  }
}

const NumVector<double>& Decomposite2D::getCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

NumVector<double>& Decomposite2D::accessCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

const NumVector<double>& Decomposite2D::getErrors() {
  if (change) computeCoeffs();
  if (errorVector.size() != nCoeffs) errorVector = NumVector<double>(nCoeffs);
  NumMatrix<double> MtNoise_1;
  MtNoise_1 = Mt*M;
  MtNoise_1 /= background_variance;
  MtNoise_1 = MtNoise_1.invert();
  
  //NumVector<double> errorVector(nCoeffs);
  // the errors here should be very close to background_rms
  // otherwise orthogonality is spoiled
  for (int i=0; i<errorVector.size(); i++) {
    errorVector(i) = sqrt(MtNoise_1(i,i));
  }
  return errorVector;
  //errors = NumMatrix<double>(nmax+1,nmax+1);
  //vectorMapping(errorVector,errors,nVector,nCoeffs);
}

const NumVector<double>& Decomposite2D::getModel() {
  computeModel();
  return model;
}

NumVector<double>& Decomposite2D::accessModel() {
  computeModel();
  return model;
}

const NumVector<double>& Decomposite2D::getResiduals() {
  computeResiduals();
  return residual;
}

NumVector<double>& Decomposite2D::accessResiduals() {
  computeResiduals();
  return residual;
}

// chi^2 for reconstruction
// see Paper III, eq. 35
double Decomposite2D::getChiSquare() {
  computeResiduals();
  double result = 0;
  // GAUSSIAN noise: no model of the objects brightness employed
  if (gaussian) {
    // rms map provided: use actual rms value in each pixel
    if (bg_rms.size() == data.size()) {
      for (int i=0; i < model.size(); i++)
	result+= gsl_pow_2(residual(i))/gsl_pow_2(bg_rms(i));
      result /= npixels - nCoeffs;
    }
    // default: rms map not provided, global variance from Object's getBackgroundRMS()
    else {
      result = residual*residual;
      result /= (background_variance*(npixels - nCoeffs));
    }
  }
  // POISSONIAN noise: pixel error already included in V_
  else {
    result = residual*(V_*residual);
    result /= npixels - nCoeffs;
  }
  return result;
}

// see Paper III, eq. 36
double Decomposite2D::getChiSquareVariance() {
  return M_SQRT2 /sqrt(npixels - nCoeffs);
}      
  
int Decomposite2D::getNMax() {
  return nmax;
}
  
void Decomposite2D::setCentroid(const Point2D& inxcentroid) {
  if (inxcentroid(0) != xcentroid(0) || inxcentroid(1) != xcentroid(1)) {
      xcentroid = inxcentroid;
      change = updateModel = updateResiduals = 1;
  }
}

void Decomposite2D::setBeta(double inbeta) {
  if (beta!=inbeta) {
    beta = inbeta;
    change = updateModel = updateResiduals = 1;
  }
}

void Decomposite2D::setNMax(int innmax) {
  if (innmax != nmax) {
    nmax = innmax;
    change = updateModel = updateResiduals = 1;
    nCoeffs = getNCoeffs(nmax);
    makeNVector(nVector,nCoeffs,nmax);
    Mt.resize(nCoeffs,npixels);
  }
}

void Decomposite2D::updateCoeffs(bool update) {
  updateC = update;
}
