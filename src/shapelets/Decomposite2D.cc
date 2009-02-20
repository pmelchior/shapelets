#include <ShapeLensConfig.h>
#include <shapelets/Decomposite2D.h>
#include <gsl/gsl_math.h>
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

Decomposite2D::Decomposite2D(const Object& O, Composite2D& C) : 
  C2D(C),
  obj(O)
{
  // set up model to obj layout
  C2D.setCentroid(obj.centroid);
  C2D.setGrid(obj.grid);
  residual.grid = obj.grid;
  residual.basefilename = obj.basefilename;

  // ways to give pixel errors (depending on noiseModel):
  // GAUSSIAN:   sigma_n -> background_variance
  // WEIGHT :    weight (inverse variance) map -> Weight
  // COVARIANCE: pixel cov. matrix -> V_
  // POISSONIAN: sigma_n + model -> Weight
  // WEIGTH_POISSONIAN: weight map + model -> Weight
  if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
    noise = 0;
    background_variance = gsl_pow_2(obj.noise_rms);
  }
  else if (ShapeLensConfig::NOISEMODEL == "WEIGHT") {
    noise = 1;
    Weight = obj.weight;
  } 
  else if (ShapeLensConfig::NOISEMODEL == "COVARIANCE") {
    noise = 2;
    V_.setCovarianceMatrix(obj.xi,obj.grid);
    V_ = V_.invert();
  }
  else if (ShapeLensConfig::NOISEMODEL == "POISSONIAN") {
    noise = 3;
    background_variance = gsl_pow_2(obj.noise_rms);
    Weight.resize(obj.size());
    // as we do not have a valid shapelet model yet, we start with a flat prior
    for (int i=0; i < obj.size(); i++)
      Weight(i) = 1./background_variance;
  }
  else if (ShapeLensConfig::NOISEMODEL == "WEIGHT_POISSONIAN") {
    noise = 4;
    // as we do not have a valid shapelet model yet, we start with the original weight map
    Weight = obj.weight;
  }
  else {
    std::cerr << "Decomposite2D: noise model '" << ShapeLensConfig::NOISEMODEL << "' unknown!" << std::endl;
    std::terminate();
  }
  fixedCoeffs = false;
}

// see Paper III, eq. 83 and following explanation
void Decomposite2D::makeLSMatrix () {
  // use Composite2D::makeShapeletMatrix() for computing M for which we have a reference
  C2D.makeShapeletMatrix();

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

  if (noise==0) {
    // of nmax has changed in C2D.makeShapeletMatrix() we have to update the 
    // covariance matrix C2D.cov
    // as basis functions are orthonormal, C2D.cov is the identity times background_variance 
    if ((C2D.M).getColumns() != (C2D.cov).getRows()) {
      int oldncoeffs = (C2D.cov).getRows();
      int newncoeffs = (C2D.M).getColumns();
      (C2D.cov).resize_clear(newncoeffs,newncoeffs);
      if (newncoeffs > oldncoeffs) {
	ublas::matrix_vector_range<ublas::matrix<data_t> > mvr (C2D.cov, ublas::range (oldncoeffs, newncoeffs), ublas::range (oldncoeffs, newncoeffs)); 
	for (unsigned int i=0; i < mvr.size(); i++)
	  mvr(i) = background_variance;
      }
    }
  } else {
    Mt = (C2D.M).transpose();
    if (noise == 1 || noise == 3 || noise == 4)
      Mt = Mt*Weight;
    if (noise == 2)
      Mt = Mt*V_;
    // as the LS matrix contains the coefficient covariance matrix,
    // we store it alredy in place
    C2D.cov = (Mt*(C2D.M)).invert();
  }
}

bool Decomposite2D::computeCoeffs() {
  if (C2D.changeM) {
    makeLSMatrix();
    // bypassing this is useful only for the regularization in OptimalDecomposite2D
    if (!fixedCoeffs) {
      // for gaussian noise (assuming a orthonormal shapelet basis)
      // we can neglect the coeff covariance matrix and perform a direct
      // projection on shapelet states.
      if (noise == 0)
	atlas::gemv(CblasTrans,1.0,(C2D.M).getMatrix(),obj,0.0,C2D.coeffs);
      // otherwise can invert the Least-Squares matrix LS
      // a Cholesky decomposition might be faster here
      else
	C2D.coeffs = C2D.cov * (Mt * obj.getNumVector());
      C2D.changeModel = 1;
    }
    return true;
  } else
    return false;
}

bool Decomposite2D::computeModel() {
  if (computeCoeffs() || C2D.changeModel) {
    C2D.evalGrid();
    return true;
  } else
    return false;
}

bool Decomposite2D::computeResiduals() {
  if (computeModel()) {
    residual = obj;
    residual -= C2D.model;
    return true;
  } else
    return false;
}

const Image<data_t>& Decomposite2D::getResiduals() {
  computeResiduals();
  return residual;
}

// chi^2 for reconstruction
// see Paper III, eq. 18
data_t Decomposite2D::getChiSquare() {
  if (computeResiduals()) {
    if (noise == 0)
      chi2 = (residual.getNumVector()*residual.getNumVector())/background_variance;
    else if (noise == 1 || noise == 3 || noise == 4)
      chi2 = residual.getNumVector()*(Weight*residual.getNumVector());
    else if (noise == 2)
      chi2 = residual.getNumVector()*(V_*residual.getNumVector());
    chi2 /= obj.size() - C2D.coeffs.size();
  }
  return chi2;
}

// see Paper III, eq. 19
data_t Decomposite2D::getChiSquareVariance() {
return M_SQRT2 /sqrt(obj.size() - C2D.coeffs.size());
}      
  
int Decomposite2D::getNMax() {
  return C2D.getNMax();
}
  
void Decomposite2D::setNMax(int nmax) {
  if (C2D.coeffs.getNMax() != nmax) {
    C2D.coeffs.setNMax(nmax);
    C2D.changeM = C2D.changeModel = 1;
  }
}

data_t Decomposite2D::getBeta() {
  return C2D.getBeta();
}

void Decomposite2D::setBeta(data_t beta) {
  C2D.setBeta(beta);
}

void Decomposite2D::fixCoeffs(bool fixed) {
  fixedCoeffs = fixed;
}

void Decomposite2D::setCoeffs(const CoefficientVector<data_t>& coeffs) {
  C2D.setCoeffs(coeffs);
}

void Decomposite2D::updateWeightMap() {
  if (noise == 3 || noise == 4) {
    computeModel();
    if (noise == 3) {
      for (int i=0; i < obj.size(); i++) {
	if ((C2D.model)(i) > 0)
	  Weight(i) = 1./(background_variance + (C2D.model)(i));
	else
	  Weight(i) = 1./background_variance;
      }
    } else {
      Weight = obj.weight;
      for (int i=0; i < obj.size(); i++)
	if ((C2D.model)(i) > 0)
	  Weight(i) = 1./((1./obj.weight(i)) + (C2D.model)(i));
    }
  }
}
