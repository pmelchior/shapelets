#include <shapelets/Composite2D.h>
// for factorial and other math functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;

Composite2D::Composite2D() : Shapelets2D() {
  change = 1;
}

Composite2D::Composite2D(const CoefficientVector<data_t>& Coeffs, data_t beta, Point2D& inxcentroid) : Shapelets2D() {  
  coeffs = Coeffs;
  xcentroid = inxcentroid;
  Shapelets2D::setBeta(beta);
  // grid has to be redefined before evaluation
  change = 1;
}

Composite2D::Composite2D(const Composite2D& source) {
  Composite2D::operator=(source);
}

Composite2D & Composite2D::operator= (const Composite2D &source) {
  coeffs = source.coeffs;
  Shapelets2D::setBeta(source.getBeta());
  grid = source.grid;
  xcentroid(0) = source.xcentroid(0);
  xcentroid(1) = source.xcentroid(1);
  return *this;
}
  
unsigned int Composite2D::getNMax() const {
  return coeffs.getNMax();
}

void  Composite2D::setCoeffs(const CoefficientVector<data_t>& newCoeffs ) {
  coeffs = newCoeffs;
  change = 1;
}

const CoefficientVector<data_t>& Composite2D::getCoeffs() const {
  return coeffs;
}

data_t Composite2D::getBeta() const {
  return Shapelets2D::getBeta();
}

void Composite2D::setBeta(data_t beta) {
  Shapelets2D::setBeta(beta);
  change = 1;
}

const Point2D& Composite2D::getCentroid() const {
  return xcentroid;
}

void Composite2D::setCentroid(const Point2D& inxcentroid) {
  xcentroid = inxcentroid;
  change = 1;
}

const Grid& Composite2D::getGrid() const {
  return grid;
}

void Composite2D::setGrid(const Grid& ingrid) {
  grid = ingrid;
  change = 1;
}

data_t Composite2D::eval(const Point2D& x) {
  data_t result = 0;
  Point2D xdiff;
  xdiff(0) = x(0) - xcentroid(0);
  xdiff(1) = x(1) - xcentroid(1);
  // result = sum over all coeffs * basis function
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) 
    result += coeffs(i) * Shapelets2D::eval(nVector.getState1(i), nVector.getState2(i), xdiff);
  return result;
}

void Composite2D::evalGrid() {
  NumMatrix<data_t> M(grid.size(),coeffs.getNCoeffs());
  makeShapeletMatrix(M);
  model = M * coeffs;
  change = 0;
}  

const NumVector<data_t>& Composite2D::getModel() {
  if (change) evalGrid();
  return model;
}

NumVector<data_t>& Composite2D::accessModel() {
  change = 0;
  return model;
}

data_t Composite2D::integrate() {
   data_t result = 0;
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) 
    result += coeffs(i) * Shapelets2D::integrate(nVector.getState1(i), nVector.getState2(i)); 
  return result;
}
data_t Composite2D::integrate(data_t xmin, data_t xmax, data_t ymin,data_t ymax) {
  data_t result = 0;
  xmin -= xcentroid(0);
  xmax -= xcentroid(0);
  ymin -= xcentroid(1);
  ymax -= xcentroid(1);
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) 
    result += coeffs(i) * Shapelets2D::integrate(nVector.getState1(i), nVector.getState2(i), xmin, xmax, ymin, ymax); 
  return result;
}

// compute flux from shapelet coeffs
// see Paper I, eq. 26
data_t Composite2D::getShapeletFlux() const {
  data_t result = 0;
  int n1, n2;
  // result = sum over n1,n2 (coeff_(n1,n2) * Integral of B_(n1,n2))
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    std::cout << i << "\t" << n1 << "\t" << n2 << "\t" << coeffs(i) << std::endl;
    if (n1%2 == 0 && n2%2 == 0) {
      result += 2 * gsl_pow_int(2,-(n1+n2)/2) * sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2)) * coeffs(i);
    }
  }
  std::cout << "FLUX = " << result  << std::endl;
  return M_SQRTPI*Shapelets2D::getBeta()*result;
}

// compute centroid position from shapelet coeffs
// see Paper I, eq. 27
Point2D Composite2D::getShapeletCentroid() const {
  Point2D xc(0,0);
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 1 && n2%2 ==0) {
      xc(0) += sqrt((data_t)n1 + 1)*pow(2,0.5*(2-n1-n2)) *
	sqrt(gsl_sf_fact(n1+1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact((n1+1)/2)*gsl_sf_fact(n2/2)) *
	coeffs(i);
    }
    if (n1%2 == 0 && n2%2 ==1) {
      xc(1) += sqrt((data_t)n2 + 1)*pow(2,0.5*(2-n1-n2)) * 
	sqrt(gsl_sf_fact(n2+1)*gsl_sf_fact(n1)) /
	(gsl_sf_fact((n2+1)/2)*gsl_sf_fact(n1/2)) *
	coeffs(i);
    }
  }
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  xc(0) = M_SQRTPI*beta*beta*xc(0)/flux;
  xc(1) = M_SQRTPI*beta*beta*xc(1)/flux;
  return xc;
}

// compute 2nd brightness moments from shapelet coeffs
NumMatrix<data_t> Composite2D::getShapelet2ndMoments() const {
  NumMatrix<data_t> Q(2,2);
  data_t factor;
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    if (n1%2 == 0 && n2%2 ==0) {
      factor = 2 * gsl_pow_int(2,-(n1+n2)/2) * 
	sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) /
	(gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2)) * coeffs(i);
      Q(0,0) += factor * (1+2*n1);
      Q(1,1) += factor * (1+2*n2);
    } else if (n1%2 == 1 && n2%2 == 1) {
      Q(0,1) += 2 * gsl_pow_int(2,-(n1+n2)/2) * sqrt((data_t)(n1+1)*(n2+1)) *
	sqrt(gsl_sf_fact(n1+1)*gsl_sf_fact(n2+1)) /
	(gsl_sf_fact((n1+1)/2)*gsl_sf_fact((n2+1)/2)) * coeffs(i);
    }
  }
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  Q(0,0) *= M_SQRTPI * beta*beta*beta / flux;
  Q(0,1) *= M_SQRTPI * beta*beta*beta / flux;
  Q(1,1) *= M_SQRTPI * beta*beta*beta / flux;
  Q(1,0) = Q(0,1);
  return Q;
}

// see Paper I, eq. (28)
data_t Composite2D::getShapeletRMSRadius() const {
  data_t rms = 0;
  int n1, n2;
  const IndexVector& nVector = coeffs.getIndexVector();
  for (unsigned int i = 0; i < nVector.getNCoeffs(); i++) {
    n1 = nVector.getState1(i);
    n2 = nVector.getState2(i);
    std::cout << i << "\t" << n1 << "\t" << n2 << "\t" << coeffs(i) << std::endl;
    if (n1%2 == 0 && n2%2 ==0) {
      rms += 4 * gsl_pow_int(2,-(n1+n2)/2) * (1+n1+n2) *
	sqrt(gsl_sf_fact(n1)*gsl_sf_fact(n2)) / (gsl_sf_fact(n1/2)*gsl_sf_fact(n2/2)) * 
	coeffs(i); 
    }
  }
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  return sqrt(rms * M_SQRTPI * beta*beta*beta / flux);
}

void Composite2D::makeShapeletMatrix(NumMatrix<data_t>& M) {
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
