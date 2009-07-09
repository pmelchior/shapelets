#include "../../include/shapelets/Composite.h"
#include "../../include/shapelets/Shapelets2D.h"
#include "../../include/ShapeLensConfig.h"
// for factorial and other math functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace shapelens;
using namespace std;

Composite::Composite() {
  changeM = changeModel = true;
}

Composite::Composite(const CoefficientVector<data_t>& Coeffs, data_t Beta, Point<data_t>& Xcentroid) :
  coeffs(Coeffs),
  beta(Beta),
  xcentroid(Xcentroid) {  
  changeM = changeModel = true;
}

unsigned int Composite::getNMax() const {
  return coeffs.getNMax();
}

void Composite::setNMax(unsigned int nmax) {
  if (coeffs.getNMax() != nmax) {
    coeffs.setNMax(nmax);
    changeM = true;
    changeModel = true;
    cov.resize(coeffs.size(),coeffs.size());
  }
}

void Composite::setCoeffs(const CoefficientVector<data_t>& newCoeffs ) {
  if (coeffs.size() != newCoeffs.size())
    changeM = true;
  coeffs = newCoeffs;
  changeModel = true;
}

const CoefficientVector<data_t>& Composite::getCoeffs() const {
  return coeffs;
}

void Composite::setCovarianceMatrix(const NumMatrix<data_t>& newCov) {
  if (coeffs.getNCoeffs() == newCov.getRows() && coeffs.getNCoeffs() == newCov.getColumns())
    cov = newCov;
  else {
    std::cerr << "Composite: covariance matrix doesn not have correct dimensions!" << std::endl;
    std::terminate();
  }
}

const NumMatrix<data_t>& Composite::getCovarianceMatrix() const {
  return cov;
}

data_t Composite::getBeta() const {
  return beta;
}

void Composite::setBeta(data_t beta_) {
  if (beta_ != beta) {
    beta = beta_;
    changeM = changeModel = true;
  }
}

const Point<data_t>& Composite::getCentroid() const {
  return xcentroid;
}

void Composite::setCentroid(const Point<data_t>& inxcentroid) {
  xcentroid = inxcentroid;
  changeM = changeModel = true;
}

const Grid& Composite::getGrid() const {
  return model.grid;
}

void Composite::setGrid(const Grid& grid) {
  model.grid = grid;
  changeM = changeModel = true;
}


void Composite::evalGrid() {
  if (changeModel) {
    makeShapeletMatrix();
    if (ShapeLensConfig::PIXEL_INTEGRATION)
      model = MInt * coeffs;
    else
      model = M * coeffs;
    changeModel = 0;
  }
}  

const Image<data_t>& Composite::getModel() {
  if (changeModel) evalGrid();
  return model;
}

data_t Composite::eval(const Point<data_t>& x, NumMatrix<data_t>* cov_est) const {
  Point<data_t> xdiff(x(0) - xcentroid(0), x(1) - xcentroid(1));
  const IndexVector& nVector = coeffs.getIndexVector();
  int nmax = nVector.getNMax();
  int nCoeffs = nVector.getNCoeffs();
  NumVector<data_t> M0(nmax+1), M1(nmax+1);
  data_t x0_scaled = xdiff(0)/beta, x1_scaled = xdiff(1)/beta;
  data_t factor0 = 1./sqrt(M_SQRTPI*beta);
  data_t factor1, factor2;
    
  // start with 0th and 1st order shapelets in each direction
  M0(0) = factor0*exp(-x0_scaled*x0_scaled/2);
  if (nmax > 0)
    M0(1) = M_SQRT2*x0_scaled*M0(0);
  M1(0) = factor0*exp(-x1_scaled*x1_scaled/2);
  if (nmax > 0) 
    M1(1) = M_SQRT2*x1_scaled*M1(0);

  // use recurrance relation to compute higher orders
  for (int n=2; n<=nmax; n++) {
    factor1 = sqrt(1./(2*n));
    factor2 =sqrt((n-1.)/n); 
    M0(n) = 2*xdiff(0)/beta*factor1*M0(n-1) - factor2*M0(n-2);
    M1(n) = 2*xdiff(1)/beta*factor1*M1(n-1) - factor2*M1(n-2);
  }

  // build tensor product of M0 and M1
  NumVector<data_t> Eval(nCoeffs);
  int n0,n1;
  for (int l = 0; l < nCoeffs; l++) {
    n0 = nVector.getState1(l);
    n1 = nVector.getState2(l);
    Eval(l) = M0(n0)*M1(n1);
  }

  // multiply with coeffs
  data_t result = Eval * coeffs;

  // calculate covariance matrix
  if (cov_est != NULL) {
    cov_est->resize(1,1);
    cov_est->operator()(0,0) = Eval * (cov * Eval);
  }
  return result;
}

data_t Composite::integrate(NumMatrix<data_t>* cov_est) const {
  const IndexVector& nVector = coeffs.getIndexVector();
  NumVector<data_t> Int(nVector.getNCoeffs());
  Shapelets2D s2d(beta);
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++)
    Int(i) = s2d.integrate(nVector.getState1(i), nVector.getState2(i));

  data_t result = Int*coeffs;
  if (cov_est != NULL) {
    cov_est->resize(1,1);
    cov_est->operator()(0,0) = Int * (cov * Int);
  }
  return result;
}

data_t Composite::integrate(const Point<data_t>& Pa, const Point<data_t>& Pb, NumMatrix<data_t>* cov_est) const {
  Point<data_t> xdiffa(Pa(0) - xcentroid(0), Pa(1) - xcentroid(1));
  Point<data_t> xdiffb(Pb(0) - xcentroid(0), Pb(1) - xcentroid(1));
  const IndexVector& nVector = coeffs.getIndexVector();
  int nmax = nVector.getNMax();
  int nCoeffs = nVector.getNCoeffs();
  // stores basis function values at xdiffa and xdiffb
  NumMatrix<data_t> M0(nmax+1,2), M1(nmax+1,2);
  // stores integral between a and b
  NumVector<data_t> I0(nmax+1), I1(nmax+1);
  data_t x0_scaled_a = xdiffa(0)/beta, x1_scaled_a = xdiffa(1)/beta;
  data_t x0_scaled_b = xdiffb(0)/beta, x1_scaled_b = xdiffb(1)/beta;
  data_t factor0 = 1./sqrt(M_SQRTPI*beta);
  data_t factor0_int = sqrt(beta*M_SQRTPI/2);
  data_t factor1, factor2;
  data_t factor1_int = -M_SQRT2*beta, factor2_int;
    
  // start with 0th and 1st order shapelets in each direction
  M0(0,0) = factor0*exp(-x0_scaled_a*x0_scaled_a/2);
  M0(0,1) = factor0*exp(-x0_scaled_b*x0_scaled_b/2);
  I0(0) = factor0_int*(gsl_sf_erf(x0_scaled_b) - gsl_sf_erf(x0_scaled_a));
  if (nmax > 0) {
    M0(1,0) = M_SQRT2*x0_scaled_a*M0(0,0);
    M0(1,1) = M_SQRT2*x0_scaled_b*M0(0,1);
    I0(1) = factor1_int*(M0(0,1) - M0(0,0));
  }
  M1(0,0) = factor0*exp(-x1_scaled_a*x1_scaled_a/2);
  M1(0,1) = factor0*exp(-x1_scaled_b*x1_scaled_b/2);
  I1(0) = factor0_int*(gsl_sf_erf(x1_scaled_b) - gsl_sf_erf(x1_scaled_a));
  if (nmax > 0) {
    M1(1,0) = M_SQRT2*x1_scaled_a*M1(0,0);
    M1(1,1) = M_SQRT2*x1_scaled_b*M1(0,1);
    I1(1) = factor1_int*(M1(0,1) - M1(0,0));
  }

  // use recurrance relation to compute higher orders
  for (int n=2; n<=nmax; n++) {
    factor1 = sqrt(1./(2*n));
    factor1_int = -beta*sqrt(2./n);
    factor2 = factor2_int = sqrt((n-1.)/n); 
    M0(n,0) = 2*xdiffa(0)/beta*factor1*M0(n-1,0) - factor2*M0(n-2,0);
    M0(n,1) = 2*xdiffb(0)/beta*factor1*M0(n-1,1) - factor2*M0(n-2,1);
    M1(n,0) = 2*xdiffa(1)/beta*factor1*M1(n-1,0) - factor2*M1(n-2,0);
    M1(n,1) = 2*xdiffb(1)/beta*factor1*M1(n-1,1) - factor2*M1(n-2,1);
    I0(n) = factor1_int*(M0(n-1,1) - M0(n-1,0)) + factor2_int*I0(n-2);
    I1(n) = factor1_int*(M1(n-1,1) - M1(n-1,0)) + factor2_int*I1(n-2);
  }

  // build tensor product of I0 and I1
  NumVector<data_t> Int(nCoeffs);
  int n0, n1;
  for (int l = 0; l < nCoeffs; l++) {
    n0 = nVector.getState1(l);
    n1 = nVector.getState2(l);
    Int(l) = I0(n0)*I1(n1);
  }
  
  // multiply with coeffs
  data_t result = Int * coeffs;

  // calculate covariance matrix
  if (cov_est != NULL) {
    cov_est->resize(1,1);
    cov_est->operator()(0,0) = Int * (cov * Int);
  }
  return result;
}

// compute flux from shapelet coeffs
// see Paper I, eq. 26
data_t Composite::getShapeletFlux(NumMatrix<data_t>* cov_est) const {
  int n1, n2;
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  NumVector<data_t> Flux(nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 == 0)
      Flux(i) = 2 * gsl_pow_int(2,-(n1+n2)/2) * sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2));
  }
  data_t result = Flux*coeffs;
  data_t factor = M_SQRTPI*beta;

  if (cov_est != NULL) {
    cov_est->resize(1,1);
    cov_est->operator()(0,0)=(factor*factor)*(Flux * (cov * Flux));
  }
  return result * factor;
}

// compute centroid position from shapelet coeffs
// see Paper I, eq. 27
Point<data_t> Composite::getShapeletCentroid(NumMatrix<data_t>* cov_est) const {
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
  Point<data_t> xc;
  data_t flux = getShapeletFlux(); // ignore covariance of flux here
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
Quadrupole Composite::getShapelet2ndMoments(NumMatrix<data_t>* cov_est) const {
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
  factor = M_SQRTPI * beta*beta*beta / flux;
  Quadrupole Q;
  Q(0,0) = result(0) * factor;
  Q(0,1) = result(1) * factor;
  Q(1,1) = result(2) * factor;
  if (cov_est != NULL) {
    cov_est->operator=((Moment * cov) * Moment.transpose());
    cov_est->operator*=(factor*factor);
  }
  return Q;
}

// see Paper I, eq. (28)
data_t Composite::getShapeletRMSRadius(NumMatrix<data_t>* cov_est) const {
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  NumVector<data_t> RMS2(nVector.getNCoeffs());
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 ==0) {
      RMS2(i) = 4 * gsl_pow_int(2,-(n1+n2)/2) * (1+n1+n2) *
	sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) / (gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2));
    }
  }
  data_t result = RMS2 * coeffs;
  data_t flux = getShapeletFlux();
  data_t factor = M_SQRTPI * beta*beta*beta / flux;
  result *= factor;
  if (cov_est != NULL) {
    cov_est->resize(1,1);
    cov_est->operator()(0,0) = RMS2 * (cov * RMS2);
    cov_est->operator()(0,0)*= (factor*factor);
    cov_est->operator()(0,0) = sqrt(cov_est->operator()(0,0));
  }

  return sqrt(result);
}

void Composite::makeShapeletMatrix() {
  if (changeM) {
    const Grid& grid = model.grid;
    if (M.getRows() != grid.size() || M.getColumns() != coeffs.size())
      M.resize(grid.size(),coeffs.getNCoeffs());

    const IndexVector& nVector = coeffs.getIndexVector();
    int nmax = nVector.getNMax();
    int nCoeffs = nVector.getNCoeffs();
    int npixels = grid.size();
    int N0 = grid.getSize(0), N1 = grid.getSize(1);
    NumMatrix<data_t> M0(nmax+1,N0), M1(nmax+1,N1);
    // start with 0th and 1st order shapelets
    data_t x0_scaled, x1_scaled;
    data_t factor0 = 1./sqrt(M_SQRTPI*beta);
    data_t factor1, factor2;
    int x,y,n0,n1;
    
    // simple sampling
    for (x=0; x< N0; x++) {
      x0_scaled = (grid(x,0) - xcentroid(0))/beta;
      M0(0,x) = factor0*exp(-x0_scaled*x0_scaled/2);
      if (nmax > 0)
	M0(1,x) = M_SQRT2*x0_scaled*M0(0,x);
    }
    for (y=0; y< N1; y++) {
      x1_scaled = (grid(y*N0,1) - xcentroid(1))/beta;
      M1(0,y) = factor0*exp(-x1_scaled*x1_scaled/2);
      if (nmax > 0) 
	M1(1,y) = M_SQRT2*x1_scaled*M1(0,y);
    }

    // use recurrance relation to compute higher orders
    for (int n=2;n<=nmax;n++) {
      factor1 = sqrt(1./(2*n));
      factor2 =sqrt((n-1.)/n); 
      for (x=0; x < N0; x++)
	M0(n,x) = 2*(grid(x,0) - xcentroid(0))/beta*factor1*M0(n-1,x) 
	  - factor2*M0(n-2,x);
      for (y=0; y< N1; y++)
	M1(n,y) = 2*(grid(y*N0,1) - xcentroid(1))/beta*factor1*M1(n-1,y) 
	  - factor2*M1(n-2,y);
    }
    // now build tensor product of M0 and M1
    for (int l = 0; l < nCoeffs; l++) {
      n0 = nVector.getState1(l);
      n1 = nVector.getState2(l);
      for (int i=0; i<npixels; i++) {
	x = i%N0;
	y = i/N0;
	M(i,l) = M0(n0,x)*M1(n1,y);
      }
    }

    // in-pixel integration:
    // according to PaperIII, eqs (30) - (34)
    if (ShapeLensConfig::PIXEL_INTEGRATION) {
      // it is possible to store this entire in matrix M
      // but then the basis functions are not orthonormal anymore
      // computation of coefficients in Decomposite must be changed then
      if (MInt.getRows() != grid.size() || MInt.getColumns() != coeffs.size())
	MInt.resize(grid.size(),coeffs.getNCoeffs());

      NumMatrix<data_t> I0(nmax+1,N0), I1(nmax+1,N1);
      data_t x0_scaled_b, x1_scaled_b;
      factor0 = sqrt(beta*M_SQRTPI/2);
      factor1 = -M_SQRT2*beta;
      for (x=0; x<N0 - 1; x++) {
	x0_scaled = (grid(x,0) - xcentroid(0))/(M_SQRT2*beta);
	x0_scaled_b = (grid(x+1,0) - xcentroid(0))/(M_SQRT2*beta);
	I0(0,x) = factor0*(gsl_sf_erf(x0_scaled_b) - gsl_sf_erf(x0_scaled));
	if (nmax > 0)
	  I0(1,x) = factor1*(M0(0,x+1) - M0(0,x));
      }
      for (y=0; y<N1 - 1; y++) {
	x1_scaled = (grid(y*N0,1) - xcentroid(1))/(M_SQRT2*beta);
	x1_scaled_b = (grid((y+1)*N0,1) - xcentroid(1))/(M_SQRT2*beta);
	I1(0,y) = factor0*(gsl_sf_erf(x1_scaled_b) - gsl_sf_erf(x1_scaled));
	if (nmax > 0)
	  I1(1,y) = factor1*(M1(0,y+1) - M1(0,y));
      }
      // use recurrance relation to compute higher orders
      for (int n=2;n<=nmax;n++) {
	factor1 = -beta*sqrt(2./n);
	factor2 = sqrt((n-1.)/n);
	for (x=0; x < N0 - 1; x++)
	  I0(n,x) = factor1*(M0(n-1,x+1) - M0(n-1,x)) + factor2*I0(n-2,x);
	for (y=0; y<N1 - 1; y++) 
	  I1(n,y) = factor1*(M1(n-1,y+1) - M1(n-1,y)) + factor2*I1(n-2,y);
      }
      // now build tensor product of I0 and I1
      for (int l = 0; l < nCoeffs; l++) {
	n0 = nVector.getState1(l);
	n1 = nVector.getState2(l);
	for (int i=0; i<npixels; i++) {
	  x = i%N0;
	  y = i/N0;
	  MInt(i,l) = I0(n0,x)*I1(n1,y);
	}
      }
    }
  }
  changeM = 0;
}
