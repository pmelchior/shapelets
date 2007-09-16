#include <shapelets/Composite1D.h>

// starting with a shapelet of order 2, since this is most common to estimate the
// large scale structure.
Composite1D::Composite1D(double beta, double inxcentroid, const NumVector<double>& coeffs) : 
Shapelets1D (2,beta) {
  shapeletCoeffs = coeffs;
  xcentroid = inxcentroid;
  // get order from lenght of given shapeletCoeffs
  order = orderlimit = shapeletCoeffs.size() - 1;
  Shapelets1D::setBeta(beta);
  Shapelets1D::setOrder(order);
  // construct default Grid from beta and order
  double stepsize = Shapelets1D::getThetaMin(order)/2;
  double range = Shapelets1D::getThetaMax(order)*2;
  grid = Grid(xcentroid-range,xcentroid+range,stepsize);
}

int Composite1D::getOrder() {
  // if orderlimit is set, this will be lower than max order
  return GSL_MIN_INT(order,orderlimit);
}

void Composite1D::setOrderLimit(int inorderlimit) {
  // only set limit if it's really lower
  if (orderlimit <= Shapelets1D::getOrder())
    orderlimit = inorderlimit;
  else std::cout << "# Warning (Composite1D): setOrderLimit attempts to increase number of orders to evaluate." << std::endl;
}

double Composite1D::getBeta() {
  return Shapelets1D::getBeta();
}

void Composite1D::setBeta(double beta) {
  Shapelets1D::setBeta(beta);
}

double Composite1D::getThetaMin() {
  return Shapelets1D::getThetaMin(orderlimit);
}

double Composite1D::getThetaMax() {
  return Shapelets1D::getThetaMax(orderlimit);
}

void Composite1D::setCoeffs(const NumVector<double>& newCoeffs) {
  shapeletCoeffs = newCoeffs;
  // reset order and limit
  order = orderlimit = newCoeffs.size() - 1;
  // maybe we need higher orders
  if (order > Shapelets1D::getOrder()) Shapelets1D::setOrder(order);
}

void Composite1D::setGrid(const Grid& ingrid) {
  grid = ingrid;
}

const Grid& Composite1D::getGrid() {
  return grid;
}

double Composite1D::eval(double x) {
  double result = 0;
  // result = sum over l (coeff_l * B_l(x))
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::eval(l, x - xcentroid);
  return result;
}

void Composite1D::evalGrid(NumVector<double>& values) {
  for (int j = 0; j < grid.size(); j +=1)
    values(j) = eval(grid(j,0));
}

double Composite1D::integrate() {
  double result = 0;
  // result = sum over l (coeff_l * Integral of B_l)
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::integrate(l);
  return result;
}

double Composite1D::integrate(double xmin, double xmax) {
  double result = 0;
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::integrate(l,xmin,xmax);
  return result;
}
  
