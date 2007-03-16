#include <Composite2D.h>
#include <MatrixManipulations.h>
// for factorial and other math functions
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

using namespace std;

Composite2D::Composite2D() : Shapelets2D() {
  lockGrid = 0;
  change = changeGrid = 1;
}

Composite2D::Composite2D(double inbeta, Point2D& inxcentroid, const NumMatrix<double>& startCoeffs) : Shapelets2D() {  
  shapeletCoeffs = startCoeffs;
  xcentroid = inxcentroid;
  order0 = orderlimit0 = shapeletCoeffs.getRows() -1;
  order1 = orderlimit1 = shapeletCoeffs.getColumns() -1;
  //beta = inbeta;
  Shapelets2D::setOrders(order0,order1);
  // grid hast to be redefined before evaluation
  change = changeGrid = 1;
  lockGrid = 0;
  // set beta in the Shapelets2D and define if in-pixel-Integration will be used
  setBeta(inbeta);
}

Composite2D & Composite2D::operator= (const Composite2D &source) {
  Shapelets2D::setOrders(source.order0,source.order1);
  Shapelets2D::setBeta(beta);
  shapeletCoeffs = source.shapeletCoeffs;
  grid = source.grid;
  xcentroid(0) = source.xcentroid(0);
  xcentroid(1) = source.xcentroid(1);
  order0 = source.order0;
  order1 = source.order1;
  orderlimit0 = source.orderlimit0;
  orderlimit1 = source.orderlimit1;
  stepsize0 = source.stepsize0;
  stepsize1 = source.stepsize1;
  return *this;
}
  
int Composite2D::getOrder(bool direction) const{
  // if orderlimit is set, this could be lower than max order
  if (direction==0) return GSL_MIN_INT(order0,orderlimit0);
  else return GSL_MIN_INT(order1,orderlimit1);
}

int Composite2D::getNMax() const {
  if (order0==order1 && orderlimit0 == orderlimit1) 
    return GSL_MIN_INT(order0,orderlimit0);
}

void Composite2D::setOrderLimit(bool direction, int orderlimit) {
  // only set orderlimit, if it's really lower then max order
  if (Shapelets2D::getOrder(direction) >= orderlimit) {
      if (direction==0) orderlimit0 = orderlimit;
      else orderlimit1 = orderlimit;
      change = 1;
  } else 
    cout << "Warning: setOrderLimit attempts to increase number of orders to evaluate." << endl << endl;
}

void  Composite2D::setCoeffs(const NumMatrix<double>& newCoeffs ) {
  shapeletCoeffs = newCoeffs;
  change = 1;
  // set new orders and limits
  order0 = orderlimit0 = shapeletCoeffs.getRows() -1;
  order1 = orderlimit1 = shapeletCoeffs.getColumns() -1;
  // if new size is different, adapt the shapelet orders
  // and set flag for redefining the grid before evaluation
  if (Shapelets2D::getOrder(0) != order0 || Shapelets2D::getOrder(1) != order1) {
    Shapelets2D::setOrders(order0,order1);
    changeGrid = 1;
  }
}

double Composite2D::getBeta() const {
  return beta;
}

double& Composite2D::accessBeta() {
  change = 1;
  return beta;
}

void Composite2D::setBeta(double inbeta) {
  beta = inbeta;
  Shapelets2D::setBeta(beta);
  change = changeGrid = 1;
}

const Point2D& Composite2D::getCentroid() const {
  return xcentroid;
}

Point2D& Composite2D::accessCentroid() {
  change = 1;
  return xcentroid;
}

void Composite2D::setCentroid(const Point2D& inxcentroid) {
  xcentroid = inxcentroid;
  change = changeGrid = 1;
}

const Grid& Composite2D::getGrid() {
  if (changeGrid) defineGrid();
  return grid;
}

void Composite2D::setGrid(const Grid& ingrid) {
  grid = ingrid;
  stepsize0 = grid.getStepsize(0);
  stepsize1 = grid.getStepsize(1);
  lockGrid = 1;
  change = 1;
}

double Composite2D::eval(const Point2D& x) {
  double result = 0;
  Point2D xdiff;
  xdiff(0) = x(0) - xcentroid(0);
  xdiff(1) = x(1) - xcentroid(1);
  // result = sum over l0,l1 (coeff_(l0,l1) * B_(l0,l1)(x0,x1))
  for (int l0 = 0; l0 <= orderlimit0; l0++)
    for (int l1 = 0; l1 <= orderlimit1; l1++)
      if (shapeletCoeffs(l0,l1) != 0)
	result += shapeletCoeffs(l0,l1) * Shapelets2D::eval(l0,l1,xdiff);
  // get rid of tiny values
  if (fabs(result) < 1e-20) result = 0;
  return result;
}

double Composite2D::evalGridPoint(const Point2D& x) {
  double result = 0;
  Point2D xdiff;
  //double range = 2 * Shapelets2D::getThetaMax(orderlimit0,orderlimit1);
  // shift to center of pixel
  //xdiff(0) = x(0) + 0.5*stepsize0 - xcentroid(0);    
  //xdiff(1) = x(1) + 0.5*stepsize1 - xcentroid(1);
  xdiff(0) = x(0) - xcentroid(0);    
  xdiff(1) = x(1) - xcentroid(1);
  for (int l0 = 0; l0 <= orderlimit0; l0++)
    for (int l1 = 0; l1 <= orderlimit1; l1++)
      if (shapeletCoeffs(l0,l1) != 0)
	result += shapeletCoeffs(l0,l1) * Shapelets2D::eval(l0,l1,xdiff);

  // get rid of tiny values
  if (fabs(result) < 1e-20) result = 0;
  return result;
}

void Composite2D::evalGrid() {
  if (changeGrid) defineGrid();
  if ((orderlimit0 != orderlimit1) || (shapeletCoeffs(orderlimit0-1,orderlimit1-1)!=0)) {
    if (model.size() != grid.size())
      model.resize(grid.size());
    for (int j = 0; j < grid.size(); j +=1)
      model(j) = evalGridPoint(grid(j));
  } else {
    // this approach only works for square, upper triangular coeff matrices
    // thus, nnmax = orderlimit0 = orderlimit1
    int nCoeffs = getNCoeffs(orderlimit0);
    NumVector<double> coeffVector(nCoeffs);
    NumMatrix<int> nVector;
    makeNVector(nVector,nCoeffs,orderlimit0);
    NumMatrix<double> M(grid.size(),nCoeffs);
    makeShapeletMatrix(M,nVector);
    matrixMapping(shapeletCoeffs,coeffVector,0,nVector,nCoeffs);
    model = M * coeffVector;
  }
  change = 0;
}  

const NumVector<double>& Composite2D::getModel() {
  if (change) evalGrid();
  return model;
}

NumVector<double>& Composite2D::accessModel() {
  change = 0;
  return model;
}

double Composite2D::integrate() {
   double result = 0;
  // result = sum over l0,l1 (coeff_(l0,l1) * Integral of B_(l0,l1))
  for (int l0 = 0; l0 <= orderlimit0; l0+=1) 
    for (int l1 = 0; l1 <= orderlimit1; l1+=1)
      if (shapeletCoeffs(l0,l1) != 0)
	result += shapeletCoeffs(l0,l1) * Shapelets2D::integrate(l0,l1);
  return result;
}
double Composite2D::integrate(double x0min, double x0max, double x1min,double x1max) {
   double result = 0;
   x0min -= xcentroid(0);
   x0max -= xcentroid(0);
   x1min -= xcentroid(1);
   x1max -= xcentroid(1);
  // result = sum over l0,l1 (coeff_(l0,l1) * Integral of B_(l0,l1))
  for (int l0 = 0; l0 <= orderlimit0; l0+=1) 
    for (int l1 = 0; l1 <= orderlimit1; l1+=1)
      if (shapeletCoeffs(l0,l1) != 0)
	result += shapeletCoeffs(l0,l1)* Shapelets2D::integrate(l0,l1,x0min,x0max,x1min,x1max);
  // get rid of tiny values
  if (result < 1e-20) result = 0;
  return result;
}

// compute flux from shapelet coeffs
// see Paper I, eq. 26
double Composite2D::getShapeletFlux() const {
  double result = 0;
  for (int l0 = 0; l0 <= orderlimit0; l0+=1) {
    for (int l1 = 0; l1 <= orderlimit1; l1+=1) {
      if (l0%2 == 0 && l1%2 ==0) {
	//result += pow(2,0.5*(2-l0-l1))* sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	result += 2 * gsl_pow_int(2,-(l0+l1)/2) * sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * shapeletCoeffs(l0,l1);
      }
    }
  }
  return M_SQRTPI*beta*result;
}

// compute centroid position from shapelet coeffs
// see Paper I, eq. 27
void Composite2D::getShapeletCentroid(Point2D& xc) const {
  xc = Point2D(0,0);
  for (int l0 = 0; l0 <= orderlimit0; l0+=1) {
    for (int l1 = 0; l1 <= orderlimit1; l1+=1) {
      if (l0%2 == 1 && l1%2 ==0) {
	xc(0) += sqrt((double)l0 + 1)*pow(2,0.5*(2-l0-l1)) *
	  sqrt(gsl_sf_fact(l0+1)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact((l0+1)/2)*gsl_sf_fact(l1/2)) *
	  shapeletCoeffs(l0,l1);
      }
      if (l0%2 == 0 && l1%2 ==1) {
	xc(1) += sqrt((double)l1 + 1)*pow(2,0.5*(2-l0-l1)) * 
	  sqrt(gsl_sf_fact(l1+1)*gsl_sf_fact(l0)) /
	  (gsl_sf_fact((l1+1)/2)*gsl_sf_fact(l0/2)) *
	  shapeletCoeffs(l0,l1);
      }
    }
  }
  double flux = getShapeletFlux();
  xc(0) = M_SQRTPI*beta*beta*xc(0)/flux;
  xc(1) = M_SQRTPI*beta*beta*xc(1)/flux;
}

// compute 2nd brightness moments from shapelet coeffs
// see IDL code: properties/shapelets_quadrupole.pro
void Composite2D::getShapelet2ndMoments(NumMatrix<double>& Q) const {
  if (Q.getRows() != 2 || Q.getColumns() != 2) 
    Q = NumMatrix<double>(2,2);
  for (int l0 = 0; l0 <= orderlimit0; l0++) {
    for (int l1 = 0; l1 <= orderlimit1; l1++) {
      if (l0%2 == 0 && l1%2 ==0) {
	Q(0,0) += 2 * gsl_pow_int(2,-(l0+l1)/2) * (1+2*l0) *
	  sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * shapeletCoeffs(l0,l1);
	Q(1,1) += 2 * gsl_pow_int(2,-(l0+l1)/2) * (1+2*l1) *
	  sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) /
	  (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * shapeletCoeffs(l0,l1);
      } else if (l0%2 == 1 && l1%2 == 1) {
	Q(0,1) += 2 * gsl_pow_int(2,-(l0+l1)/2) * sqrt((double)(l0+1)*(l1+1)) *
	  sqrt(gsl_sf_fact(l0+1)*gsl_sf_fact(l1+1)) /
	  (gsl_sf_fact((l0+1)/2)*gsl_sf_fact((l1+1)/2)) * shapeletCoeffs(l0,l1);
      }
    }
  }
  double flux = getShapeletFlux();
  Q(0,0) *= M_SQRTPI * beta*beta*beta / flux;
  Q(0,1) *= M_SQRTPI * beta*beta*beta / flux;
  Q(1,1) *= M_SQRTPI * beta*beta*beta / flux;
  Q(1,0) = Q(0,1);
}

// see Paper I, eq. (28)
double Composite2D::getShapeletsRMSRadius() const {
  double rms = 0;
  for (int l0 = 0; l0 <= orderlimit0; l0++) {
    for (int l1 = 0; l1 <= orderlimit1; l1++) {
      if (l0%2 == 0 && l1%2 ==0) {
	rms += 4 * gsl_pow_int(2,-(l0+l1)/2) * (1+l0+l1) *
	  sqrt(gsl_sf_fact(l0)*gsl_sf_fact(l1)) / (gsl_sf_fact(l0/2)*gsl_sf_fact(l1/2)) * 
	  shapeletCoeffs(l0,l1); 
      }
    }
  }
  double flux = getShapeletFlux();
  return sqrt(rms * M_SQRTPI * beta*beta*beta / flux);
}

// define Grid in such a way that it contains all the large scale structure
// and all details contained in a shapelet image
// therefore use 2*theta_max as image dimension and theta_min/2 as stepsize of the grid points
void Composite2D::defineGrid() {
  changeGrid = 0;
  if (!lockGrid) {
    stepsize0 = stepsize1 = Shapelets2D::getThetaMin(order0,order1)/2;
    double range = Shapelets2D::getThetaMax(order0,order1)*3; 
    grid = Grid(xcentroid(0)-range,xcentroid(0)+range,stepsize0,xcentroid(1)-range,xcentroid(1)+range,stepsize1);
  }
}

void Composite2D::makeShapeletMatrix(NumMatrix<double>& M, NumMatrix<int>& nVector) {
  int nmax = orderlimit0;
  int nCoeffs = getNCoeffs(nmax);
  int npixels = grid.size();
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
    n0 = getN1(nVector,l);
    n1 = getN2(nVector,l);
    for (int i=0; i<npixels; i++)
      // this access scheme is probably slow but we come around the matrix Mt
      M(i,l) = M0(n0,i)*M1(n1,i);
  }
}
