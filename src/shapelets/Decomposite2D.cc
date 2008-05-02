#include <ShapeLensConfig.h>
#include <shapelets/Decomposite2D.h>
#include <gsl/gsl_math.h>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

Decomposite2D::Decomposite2D(int innmax, data_t inbeta, const Object& O) : 
obj(O) {
  nmax = innmax;
  coeffVector.setNMax(nmax);
  nCoeffs = coeffVector.getNCoeffs();
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
  const IndexVector& nVector = coeffVector.getIndexVector();
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
    n0 = nVector.getState1(l);
    n1 = nVector.getState2(l);
    for (int i=0; i<npixels; i++) 
      Mt(l,i) = M0(n0,i)*M1(n1,i);
  }
  
  // for gaussian noise and within reasonable bounds on nmax and beta
  // the coefficient covariance matrix is VERY close to the identity matrix.
  // We can therefore skip the whole computation of (Mt*M).invert().
  // Even when shapelet basis becomes non-orthonormal, we can restrict ourselves
  // to the following code:
  // if (noise == 0) {
    //if (LS.getRows() != nCoeffs) {
    //  LS.resize(nCoeffs,nCoeffs);
    //  LS.clear();
    //}
    // use explicitly multiplication with tranpose of Mt
    //atlas::gemm(CblasNoTrans,CblasTrans,1.0,Mt.getMatrix(),Mt.getMatrix(),0.0,LS.getMatrix());
  // }
  if (noise==0);
  else {
    M = Mt.transpose();
    if (noise == 1 || noise == 3)
      Mt = Mt*Weight;
    if (noise == 2)
      Mt = Mt*V_;
    LS = Mt*M;
  }
}

void Decomposite2D::computeCoeffs() {
  makeLSMatrix();
  // this is useful only for the regularization in OptimalDecomposite2D
  if (updateC) {
    // for gaussian noise (assuming a orthonormal shapelet basis)
    // we can neglect the coeff covariance matrix and perform a direct
    // projection on shapelet states.
    if (noise == 0)
      coeffVector = Mt * obj.getNumVector();
    // otherwise can invert the Least-Squares matrix LS
    // a Cholesky decomposition might be faster here
    else
      coeffVector = LS.invert() * (Mt * obj.getNumVector());
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
    // for gaussian noise model: M is not computed, so use transpose of Mt here
    if (noise == 0) {
      if (model.size() != Mt.getColumns())
	model.resize(Mt.getColumns());
      atlas::gemv(CblasTrans,1.0,Mt.getMatrix(),coeffVector,0.0,model);
    }
    else 
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

const CoefficientVector<data_t>& Decomposite2D::getCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

CoefficientVector<data_t>& Decomposite2D::accessCoeffs() {
  // are shapelet coeffs still effective, otherwise compute new ones
  if (change) computeCoeffs();
  return coeffVector;
}

NumMatrix<data_t> Decomposite2D::getCovarianceMatrix() {
  if (change) computeCoeffs();
  if (noise==0) {
    NumMatrix<data_t> identity(nCoeffs,nCoeffs);
    ublas::matrix_vector_range<ublas::matrix<data_t> > mvr (identity, ublas::range (0, nCoeffs-1), ublas::range (0, nCoeffs-1)); 
    for (unsigned int i=0; i < nCoeffs; i++)
      mvr(i) = 1;
    return identity;
  }
  else
    return LS.invert();
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
    coeffVector.setNMax(nmax);
    nCoeffs = coeffVector.getNCoeffs();
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
