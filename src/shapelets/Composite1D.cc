#include <shapelets/Composite1D.h>

// starting with a shapelet of order 2, since this is most common to estimate the
// large scale structure.
Composite1D::Composite1D(data_t beta, data_t inxcentroid, const NumVector<data_t>& coeffs) : 
Shapelets1D (2,beta) {
  shapeletCoeffs = coeffs;
  xcentroid = inxcentroid;
  // get order from lenght of given shapeletCoeffs
  order = orderlimit = shapeletCoeffs.size() - 1;
  Shapelets1D::setBeta(beta);
  Shapelets1D::setOrder(order);
  // construct default Grid from beta and order
  data_t stepsize = Shapelets1D::getThetaMin(order)/2;
  data_t range = Shapelets1D::getThetaMax(order)*2;
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

data_t Composite1D::getBeta() {
  return Shapelets1D::getBeta();
}

void Composite1D::setBeta(data_t beta) {
  Shapelets1D::setBeta(beta);
}

data_t Composite1D::getThetaMin() {
  return Shapelets1D::getThetaMin(orderlimit);
}

data_t Composite1D::getThetaMax() {
  return Shapelets1D::getThetaMax(orderlimit);
}

void Composite1D::setCoeffs(const NumVector<data_t>& newCoeffs) {
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

data_t Composite1D::eval(data_t x) {
  data_t result = 0;
  // result = sum over l (coeff_l * B_l(x))
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::eval(l, x - xcentroid);
  return result;
}

void Composite1D::evalGrid(NumVector<data_t>& values) {
  for (int j = 0; j < grid.size(); j +=1)
    values(j) = eval(grid(j,0));
}

data_t Composite1D::integrate() {
  data_t result = 0;
  // result = sum over l (coeff_l * Integral of B_l)
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::integrate(l);
  return result;
}

data_t Composite1D::integrate(data_t xmin, data_t xmax) {
  data_t result = 0;
  for (int l = 0; l <= orderlimit; l+=1)
    result += shapeletCoeffs(l) * Shapelets1D::integrate(l,xmin,xmax);
  return result;
}
  
