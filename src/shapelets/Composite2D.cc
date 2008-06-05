#include <shapelets/Composite2D.h>
// for factorial and other math functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;

Composite2D::Composite2D() : Shapelets2D() {
  changeM = changeModel = true;
}

Composite2D::Composite2D(const CoefficientVector<data_t>& Coeffs, data_t beta, Point2D& inxcentroid) : Shapelets2D() {  
  coeffs = Coeffs;
  xcentroid = inxcentroid;
  Shapelets2D::setBeta(beta);
  changeM = changeModel = true;
}

Composite2D::Composite2D(const Composite2D& source) {
  Composite2D::operator=(source);
}

Composite2D & Composite2D::operator= (const Composite2D &source) {
  coeffs = source.coeffs;
  cov = source.cov;
  Shapelets2D::setBeta(source.getBeta());
  grid = source.grid;
  xcentroid(0) = source.xcentroid(0);
  xcentroid(1) = source.xcentroid(1);
  changeM = source.changeM;
  changeModel = source.changeModel;
  return *this;
}
  
unsigned int Composite2D::getNMax() const {
  return coeffs.getNMax();
}

void Composite2D::setCoeffs(const CoefficientVector<data_t>& newCoeffs ) {
  if (coeffs.size() != newCoeffs.size())
    changeM = true;
  coeffs = newCoeffs;
  changeModel = true;
}

const CoefficientVector<data_t>& Composite2D::getCoeffs() const {
  return coeffs;
}

void Composite2D::setCovarianceMatrix(const NumMatrix<data_t>& newCov) {
  if (coeffs.getNCoeffs() == newCov.getRows() && coeffs.getNCoeffs() == newCov.getColumns())
    cov = newCov;
  else {
    std::cerr << "Composite2D: covariance matrix doesn not have correct dimensions!" << std::endl;
    std::terminate();
  }
}

const NumMatrix<data_t>& Composite2D::getCovarianceMatrix() const {
  return cov;
}

data_t Composite2D::getBeta() const {
  return Shapelets2D::getBeta();
}

void Composite2D::setBeta(data_t beta) {
  if (Shapelets2D::getBeta() != beta) {
    Shapelets2D::setBeta(beta);
    changeM = changeModel = true;
  }
}

const Point2D& Composite2D::getCentroid() const {
  return xcentroid;
}

void Composite2D::setCentroid(const Point2D& inxcentroid) {
  xcentroid = inxcentroid;
  changeM = changeModel = true;
}

const Grid& Composite2D::getGrid() const {
  return grid;
}

void Composite2D::setGrid(const Grid& ingrid) {
  grid = ingrid;
  changeM = changeModel = true;
}


void Composite2D::evalGrid() {
  if (changeModel) {
    makeShapeletMatrix();
    model = M * coeffs;
    changeModel = 0;
  }
}  

const NumVector<data_t>& Composite2D::getModel() {
  if (changeModel) evalGrid();
  return model;
}

data_t Composite2D::eval(const Point2D& x, NumMatrix<data_t>* cov_est) {
  Point2D xdiff(x(0) - xcentroid(0), x(1) - xcentroid(1));
  // result = sum over all coeffs * basis function
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Eval(1,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) 
    Eval(0,i) = Shapelets2D::eval(nVector.getState1(i), nVector.getState2(i), xdiff);
  NumVector<data_t> result = Eval * coeffs;
  if (cov_est != NULL)
    cov_est->operator=((Eval * cov) * Eval.transpose());
  return result(0);
}

data_t Composite2D::integrate(NumMatrix<data_t>* cov_est) {
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Int(1,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++)
    Int(0,i) = Shapelets2D::integrate(nVector.getState1(i), nVector.getState2(i));
  NumVector<data_t> result = Int*coeffs;
  if (cov_est != NULL)
    cov_est->operator=((Int * cov) * Int.transpose());
  return result(0);
}

data_t Composite2D::integrate(data_t xmin, data_t xmax, data_t ymin,data_t ymax, NumMatrix<data_t>* cov_est) {
  xmin -= xcentroid(0);
  xmax -= xcentroid(0);
  ymin -= xcentroid(1);
  ymax -= xcentroid(1);
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Int(1,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) 
    Int(0,i) = Shapelets2D::integrate(nVector.getState1(i), nVector.getState2(i), xmin, xmax, ymin, ymax); 
  NumVector<data_t> result = Int*coeffs;
  if (cov_est != NULL)
    cov_est->operator=((Int * cov) * Int.transpose());
  return result(0);
}

// compute flux from shapelet coeffs
// see Paper I, eq. 26
data_t Composite2D::getShapeletFlux(NumMatrix<data_t>* cov_est) const {
  int n1, n2;
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Flux(1,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 == 0)
      Flux(0,i) = 2 * gsl_pow_int(2,-(n1+n2)/2) * sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2));
  }
  NumVector<data_t> result = Flux*coeffs;
  data_t factor = M_SQRTPI*Shapelets2D::getBeta();

  if (cov_est != NULL) {
    cov_est->operator=((Flux * cov) * Flux.transpose());
    cov_est->operator*=(factor*factor);
  }

  return result(0) * factor;
}

// compute centroid position from shapelet coeffs
// see Paper I, eq. 27
Point2D Composite2D::getShapeletCentroid(NumMatrix<data_t>* cov_est) const {
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Centroid(2,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 1 && n2%2 ==0) {
      Centroid(0,i) += sqrt((data_t)n1 + 1)*pow(2,0.5*(2-n1-n2)) *
	sqrt(gsl_sf_fact(n1+1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact((n1+1)/2)*gsl_sf_fact(n2/2));
    }
    if (n1%2 == 0 && n2%2 ==1) {
      Centroid(1,i) += sqrt((data_t)n2 + 1)*pow(2,0.5*(2-n1-n2)) * 
	sqrt(gsl_sf_fact(n2+1)*gsl_sf_fact(n1)) /
	(gsl_sf_fact((n2+1)/2)*gsl_sf_fact(n1/2));
    }
  }
  NumVector<data_t> result = Centroid * coeffs;
  Point2D xc;
  data_t flux = getShapeletFlux(); // ignore covariance of flux here
  data_t beta = Shapelets2D::getBeta();
  data_t factor = M_SQRTPI*beta*beta/flux;
  xc(0) = result(0) * factor;
  xc(1) = result(1) * factor;

  if (cov_est != NULL) {
    cov_est->operator=((Centroid * cov) * Centroid.transpose());
    cov_est->operator*=(factor*factor);
  }

  return xc;
}

// compute 2nd brightness moments from shapelet coeffs
NumMatrix<data_t> Composite2D::getShapelet2ndMoments(NumMatrix<data_t>* cov_est) const {
  data_t factor;
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> Moment(3,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 ==0) {
      factor = 2 * gsl_pow_int(2,-(n1+n2)/2) * 
	sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2));
      Moment(0,i) = factor * (1+2*n1);
      Moment(2,i) = factor * (1+2*n2);
    } else if (n1%2 == 1 && n2%2 == 1) {
      Moment(1,i) = 2 * gsl_pow_int(2,-(n1+n2)/2) * sqrt((data_t)(n1+1)*(n2+1)) *
	sqrt(gsl_sf_fact(n1+1)*gsl_sf_fact(n2+1)) /
	(gsl_sf_fact((n1+1)/2)*gsl_sf_fact((n2+1)/2));
    }
  }
  NumVector<data_t> result = Moment * coeffs;
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  factor = M_SQRTPI * beta*beta*beta / flux;
  NumMatrix<data_t> Q(2,2);
  Q(0,0) = result(0) * factor;
  Q(0,1) = result(1) * factor;
  Q(1,1) = result(2) * factor;
  Q(1,0) = Q(0,1);
  if (cov_est != NULL) {
    cov_est->operator=((Moment * cov) * Moment.transpose());
    cov_est->operator*=(factor*factor);
  }
  return Q;
}

// see Paper I, eq. (28)
data_t Composite2D::getShapeletRMSRadius(NumMatrix<data_t>* cov_est) const {
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  NumMatrix<data_t> RMS2(1,nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 ==0) {
      RMS2(0,i) = 4 * gsl_pow_int(2,-(n1+n2)/2) * (1+n1+n2) *
	sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) / (gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2));
    }
  }
  NumVector<data_t> result = RMS2 * coeffs;
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  data_t factor = M_SQRTPI * beta*beta*beta / flux;
  result(0) *= factor;
  if (cov_est != NULL) {
    cov_est->operator=((RMS2 * cov) * RMS2.transpose());
    cov_est->operator*=(factor*factor);
    cov_est->operator()(0,0)=sqrt(cov_est->operator()(0,0));
  }

  return sqrt(result(0));
}

void Composite2D::makeShapeletMatrix() {
  if (changeM) {
    if (M.getRows() != grid.size() || M.getColumns() != coeffs.size())
      M.resize(grid.size(),coeffs.getNCoeffs());

    const IndexVector& nVector = coeffs.getIndexVector();
    int nmax = nVector.getNMax();
    int nCoeffs = nVector.getNCoeffs();
    int npixels = grid.size();
    NumMatrix<data_t> M0(nmax+1,npixels), M1(nmax+1,npixels);
    // start with 0th and 1st order shapelets
    data_t x0_scaled, x1_scaled;
    data_t beta = Shapelets2D::getBeta();
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
	// this access scheme is probably slow but we come around the matrix Mt
	M(i,l) = M0(n0,i)*M1(n1,i);
    }
  }
  changeM = 0;
}
