#include <shapelets/Composite2D.h>
#include <shapelets/CoefficientVector.h>
// for factorial and other math functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;

Composite2D::Composite2D() : Shapelets2D() {
  change = 1;
}

Composite2D::Composite2D(data_t beta, Point2D& inxcentroid, const NumMatrix<data_t>& Coeffs) : Shapelets2D() {  
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
  
int Composite2D::getOrder(bool direction) const{
  if (direction==0) return coeffs.getRows()-1;
  else return coeffs.getColumns()-1;
}

int Composite2D::getNMax() const {
  return GSL_MAX_INT(coeffs.getRows()-1,coeffs.getColumns()-1);
}

void  Composite2D::setCoeffs(const NumMatrix<data_t>& newCoeffs ) {
  coeffs = newCoeffs;
  change = 1;
}

const NumMatrix<data_t>& Composite2D::getCoeffs() const {
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
  // result = sum over l0,l1 (coeff_(l0,l1) * B_(l0,l1)(x0,x1))
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0++)
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1++)
      if (coeffs(l0,l1) != 0)
	result += coeffs(l0,l1) * Shapelets2D::eval(l0,l1,xdiff);
  // get rid of tiny values
  if (fabs(result) < 1e-20) result = 0;
  return result;
}

data_t Composite2D::evalGridPoint(const Point2D& x) {
  data_t result = 0;
  Point2D xdiff;
  // shift to center of pixel
  //xdiff(0) = x(0) + 0.5*stepsize0 - xcentroid(0);    
  //xdiff(1) = x(1) + 0.5*stepsize 1 - xcentroid(1);
  xdiff(0) = x(0) - xcentroid(0);    
  xdiff(1) = x(1) - xcentroid(1);
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0++)
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1++)
      if (coeffs(l0,l1) != 0)
	result += coeffs(l0,l1) * Shapelets2D::eval(l0,l1,xdiff);

  // get rid of tiny values
  if (fabs(result) < 1e-20) result = 0;
  return result;
}

void Composite2D::evalGrid() {
  if ((coeffs.getRows() - 1 != coeffs.getColumns() - 1) || (coeffs(coeffs.getRows() - 1-1,coeffs.getColumns() - 1-1)!=0)) {
    if (model.size() != grid.size())
      model.resize(grid.size());
    for (int j = 0; j < grid.size(); j +=1)
      model(j) = evalGridPoint(grid(j));
  } else {
    // this approach only works for square, upper triangular coeff matrices
    // which is assumed for this ansatz here.
    CoefficientVector<data_t> coeffVector(coeffs);
    const IndexVector& nVector = coeffVector.getIndexVector();
    NumMatrix<data_t> M(grid.size(),nVector.getNCoeffs());
    makeShapeletMatrix(M,nVector);
    model = M * (NumVector<data_t>) coeffVector;
    
  }
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
  // result = sum over l0,l1 (coeff_(l0,l1) * Integral of B_(l0,l1))
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0+=1) 
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1+=1)
      if (coeffs(l0,l1) != 0)
	result += coeffs(l0,l1) * Shapelets2D::integrate(l0,l1);
  return result;
}
data_t Composite2D::integrate(data_t x0min, data_t x0max, data_t x1min,data_t x1max) {
   data_t result = 0;
   x0min -= xcentroid(0);
   x0max -= xcentroid(0);
   x1min -= xcentroid(1);
   x1max -= xcentroid(1);
  // result = sum over l0,l1 (coeff_(l0,l1) * Integral of B_(l0,l1))
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0+=1) 
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1+=1)
      if (coeffs(l0,l1) != 0)
	result += coeffs(l0,l1)* Shapelets2D::integrate(l0,l1,x0min,x0max,x1min,x1max);
  // get rid of tiny values
  if (result < 1e-20) result = 0;
  return result;
}

// compute flux from shapelet coeffs
// see Paper I, eq. 26
data_t Composite2D::getShapeletFlux() const {
  data_t result = 0;
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0+=1) {
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1+=1) {
      if (l0%2 == 0 && l1%2 ==0) {
	//result += pow(2,0.5*(2-l0-l1))* sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	result += 2 * gsl_pow_int(2,-(l0+l1)/2) * sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * coeffs(l0,l1);
      }
    }
  }
  return M_SQRTPI*Shapelets2D::getBeta()*result;
}

// compute centroid position from shapelet coeffs
// see Paper I, eq. 27
Point2D Composite2D::getShapeletCentroid() const {
  Point2D xc(0,0);
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0+=1) {
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1+=1) {
      if (l0%2 == 1 && l1%2 ==0) {
	xc(0) += sqrt((data_t)l0 + 1)*pow(2,0.5*(2-l0-l1)) *
	  sqrt(gsl_sf_fact(l0+1)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact((l0+1)/2)*gsl_sf_fact(l1/2)) *
	  coeffs(l0,l1);
      }
      if (l0%2 == 0 && l1%2 ==1) {
	xc(1) += sqrt((data_t)l1 + 1)*pow(2,0.5*(2-l0-l1)) * 
	  sqrt(gsl_sf_fact(l1+1)*gsl_sf_fact(l0)) /
	  (gsl_sf_fact((l1+1)/2)*gsl_sf_fact(l0/2)) *
	  coeffs(l0,l1);
      }
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
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0++) {
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1++) {
      if (l0%2 == 0 && l1%2 ==0) {
	factor = 2 * gsl_pow_int(2,-(l0+l1)/2) * 
	  sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * coeffs(l0,l1);
	Q(0,0) += factor * (1+2*l0);
	Q(1,1) += factor * (1+2*l1);
      } else if (l0%2 == 1 && l1%2 == 1) {
	Q(0,1) += 2 * gsl_pow_int(2,-(l0+l1)/2) * sqrt((data_t)(l0+1)*(l1+1)) *
	  sqrt(gsl_sf_fact(l0+1)*gsl_sf_fact(l1+1)) /
	  (gsl_sf_fact((l0+1)/2)*gsl_sf_fact((l1+1)/2)) * coeffs(l0,l1);
      }
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
  for (int l0 = 0; l0 <= coeffs.getRows() - 1; l0++) {
    for (int l1 = 0; l1 <= coeffs.getColumns() - 1; l1++) {
      if (l0%2 == 0 && l1%2 ==0) {
	rms += 4 * gsl_pow_int(2,-(l0+l1)/2) * (1+l0+l1) *
	  sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) / (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * 
	  coeffs(l0,l1); 
      }
    }
  }
  data_t flux = getShapeletFlux();
  data_t beta = Shapelets2D::getBeta();
  return sqrt(rms * M_SQRTPI * beta*beta*beta / flux);
}

void Composite2D::makeShapeletMatrix(NumMatrix<data_t>& M, const IndexVector& nVector) {
  int nmax = coeffs.getRows() - 1;
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
    n0 = nVector.getN1(l);
    n1 = nVector.getN2(l);
    for (int i=0; i<npixels; i++)
      // this access scheme is probably slow but we come around the matrix Mt
      M(i,l) = M0(n0,i)*M1(n1,i);
  }
}
