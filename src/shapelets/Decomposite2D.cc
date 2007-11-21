#include <ShapeLensConfig.h>
#include <shapelets/Decomposite2D.h>
#include <shapelets/MatrixManipulations.h>
#include <gsl/gsl_math.h>

namespace ublas = boost::numeric::ublas;

Decomposite2D::Decomposite2D(int innmax, data_t inbeta, const Object& O) : 
obj(O) {
  nmax = innmax;
  nVector = IndexVector(nmax);
  nCoeffs = nVector.getNCoeffs();
  beta = inbeta;
  npixels = obj.size();
  Mt.resize(nCoeffs,npixels);

  // ways to give pixel errors (depending on noiseModel):
  // GAUSSIAN:   sigma_n -> background_variance
  // WEIGHT :    weight (inverse variance) map -> weight
  // COVARIANCE: pixel cov. matrix -> full cov. matrix
  // POISSONIAN: sigma_n + data -> weight
  if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
    noise = 0;
    background_variance = gsl_pow_2(obj.getNoiseRMS());
  }
  else if (ShapeLensConfig::NOISEMODEL == "WEIGHT") {
    noise = 1;
    Weight = obj.getWeightMap();
  } 
  else if (ShapeLensConfig::NOISEMODEL == "COVARIANCE") {
    noise = 2;
    V_ = obj.getPixelCovarianceMatrix().invert(); 
  }
  else if (ShapeLensConfig::NOISEMODEL == "POISSONIAN") {
    noise = 3;
    background_variance = gsl_pow_2(obj.getNoiseRMS());
    Weight.resize(npixels);
    for (int i=0; i < npixels; i++)
      Weight(i) = 1./background_variance;
  }
  change = updateC = updateModel = updateResiduals = 1;
}

// see Paper III, eq. 83 and following explanation
void Decomposite2D::makeLSMatrix () {
  const Point2D& xcentroid = obj.getCentroid();
  const Grid& grid = obj.getGrid();
  NumMatrix<data_t> M0(nmax+1,npixels), M1(nmax+1,npixels);
  // start with 0th and 1st order shapelets
  data_t x0_scaled, x1_scaled;
  data_t factor0 = 1./sqrt(M_SQRTPI*beta);
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
    data_t factor1 = sqrt(1./(2*n)), factor2 =sqrt((n-1.)/n); 
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
    coeffVector = LS * (const NumVector<data_t>)obj;
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

const NumVector<data_t>& Decomposite2D::getCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

NumVector<data_t>& Decomposite2D::accessCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

const NumVector<data_t>& Decomposite2D::getErrors() {
  if (change) computeCoeffs();
  if (errorVector.size() != nCoeffs) errorVector = NumVector<data_t>(nCoeffs);
  NumMatrix<data_t> MtNoise_1;
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

const NumVector<data_t>& Decomposite2D::getModel() {
  computeModel();
  return model;
}

NumVector<data_t>& Decomposite2D::accessModel() {
  computeModel();
  return model;
}

const NumVector<data_t>& Decomposite2D::getResiduals() {
  computeResiduals();
  return residual;
}

NumVector<data_t>& Decomposite2D::accessResiduals() {
  computeResiduals();
  return residual;
}

// chi^2 for reconstruction
// see Paper III, eq. 18
data_t Decomposite2D::getChiSquare() {
  computeResiduals();
  data_t result = 0;
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
data_t Decomposite2D::getChiSquareVariance() {
  return M_SQRT2 /sqrt(npixels - nCoeffs);
}      
  
int Decomposite2D::getNMax() {
  return nmax;
}
  
void Decomposite2D::setBeta(data_t inbeta) {
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

void Decomposite2D::updateWeightMap() {
  if (ShapeLensConfig::NOISEMODEL == "POISSONIAN") {
    computeModel();
    for (int i=0; i < npixels; i++) {
      if (model(i) > 0)
	Weight(i) = 1./(background_variance + model(i));
      else
	Weight(i) = 1./background_variance;
    }
  }
}
