#include <shapelets/Shapelets2D.h>
#include <gsl/gsl_math.h>

using namespace shapelens;

Shapelets2D::Shapelets2D () {
}

Shapelets2D::Shapelets2D (data_t beta) {
  S1D.setBeta(beta);
}

data_t Shapelets2D::getBeta() const {
  return S1D.getBeta();
}

void Shapelets2D::setBeta(data_t beta) {
  S1D.setBeta(beta);
}

data_t Shapelets2D::getThetaMin(int order0, int order1) const {
  return S1D.getBeta()/sqrt((data_t) GSL_MAX_INT(order0,order1)+1);
}

data_t Shapelets2D::getThetaMax(int order0, int order1) const {
  return S1D.getBeta()*sqrt((data_t) GSL_MAX_INT(order0,order1)+1);
}

data_t Shapelets2D::integrate(int order0, int order1) const {
  return S1D.integrate(order0)*S1D.integrate(order1);
}

// integrate basis function within range
// see Paper III. eq. (82)
data_t Shapelets2D::integrate(int order0, int order1, data_t x0min, data_t x0max, data_t x1min,data_t x1max) {
  return S1D.integrate(order0,x0min,x0max)*S1D.integrate(order1,x1min,x1max);
}

data_t Shapelets2D::eval (int order0, int order1, Point2D<data_t>& x) {
  return S1D.eval(order0,x(0)) * S1D.eval(order1,x(1));
}
