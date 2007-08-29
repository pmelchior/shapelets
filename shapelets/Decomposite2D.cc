#include <Decomposite2D.h>
#include <MatrixManipulations.h>
#include <IO.h>
#include <gsl/gsl_math.h>

namespace ublas = boost::numeric::ublas;

Decomposite2D::Decomposite2D(int innmax, double inbeta, const Object& O) : 
obj(O) {
  nmax = innmax;
  nVector = IndexVector(nmax);
  nCoeffs = nVector.getNCoeffs();
  setBeta(inbeta);
  npixels = obj.size();
  Mt.resize(nCoeffs,npixels);

  // ways to give pixel errors (depending on noiseModel):
  // GAUSSIAN:   sigma_n -> background_variance
  // WEIGHT :    weight (inverse variance) map -> weight
  // COVARIANCE: pixel cov. matrix -> full cov. matrix
  // POISSONIAN: sigma_n + data -> weight
  if (obj.getNoiseModel().compare("GAUSSIAN")==0) {
    noise = 0;
    background_variance = gsl_pow_2(obj.getNoiseRMS());
  }
  else if (obj.getNoiseModel().compare("WEIGHT")==0) {
    noise = 1;
    Weight = obj.getWeightMap();
  } else if (obj.getNoiseModel().compare("COVARIANCE")==0) {
    noise = 2;
    V_ = obj.getPixelCovarianceMatrix().invert(); 
  }
  else if (obj.getNoiseModel().compare("POISSONIAN")==0) {
    noise = 3;
    NumVector<double> tmp = obj;
    Weight.resize(npixels);
    for (int i=0; i < npixels; i++)
      tmp(i) += background_variance;
    convolveGaussian(tmp,Weight,obj.getSize(0),obj.getSize(1));
    for (int i=0; i< npixels; i++)
      Weight(i) = 1./Weight(i);
  }
  change = updateC = updateModel = updateResiduals = 1;
}

// see Paper III, eq. 83 and following explanation
void Decomposite2D::makeLSMatrix () {
  const Point2D& xcentroid = obj.getCentroid();
  const Grid& grid = obj.getGrid();
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
  int n0, n1;
  for (int l = 0; l < nCoeffs; l++) {
    n0 = nVector.getN1(l);
    n1 = nVector.getN2(l);
    for (int i=0; i<npixels; i++) 
      Mt(l,i) = M0(n0,i)*M1(n1,i);
  }
  M = Mt.transpose();
  
  if (noise == 0)
    Mt /= background_variance;
  else if (noise == 1 || noise == 3)
    Mt = Mt*Weight;
  else if (noise == 2)
    Mt = Mt*V_;
   
  LS = Mt*M;
  LS = LS.invert();
  LS = LS*Mt;
  

  // SVD method 
  //LS = M.svd_invert(); 
  // overlapping integrals
  //LS = M.transpose();
}

void Decomposite2D::computeCoeffs() {
  makeLSMatrix();
  // this is useful only for the regularization in OptimalDecomposite2D
  if (updateC) {
    coeffVector = LS * (const NumVector<double>)obj;
  }
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
    residual = obj;
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
  // since Mt here is in fact Mt*V_, we can compute the covariance matrix
  // of the coefficients simply by Mt*M
  MtNoise_1 = (Mt*M).invert();
  
  // the errors here should be very close to background_rms
  // otherwise orthogonality is spoiled
  for (int i=0; i<errorVector.size(); i++) {
    errorVector(i) = sqrt(MtNoise_1(i,i));
  }
  return errorVector;
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
// see Paper III, eq. 18
double Decomposite2D::getChiSquare() {
  computeResiduals();
  double result = 0;
  if (noise == 0)
    result = (residual*residual)/background_variance;
  else if (noise == 1 || noise == 3)
    result = residual*(Weight*residual);
  else if (noise == 2)
    result = residual*(V_*residual);
  result /= npixels - nCoeffs;
  return result;
}

// see Paper III, eq. 19
double Decomposite2D::getChiSquareVariance() {
  return M_SQRT2 /sqrt(npixels - nCoeffs);
}      
  
int Decomposite2D::getNMax() {
  return nmax;
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
    nVector = IndexVector(nmax);
    nCoeffs = nVector.getNCoeffs();
    Mt.resize(nCoeffs,npixels);
  }
}

void Decomposite2D::updateCoeffs(bool update) {
  updateC = update;
}
