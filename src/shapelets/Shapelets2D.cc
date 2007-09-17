#include <shapelets/Shapelets2D.h>
#include <gsl/gsl_math.h>

Shapelets2D::Shapelets2D () {
  // just to see that they are not correctly set yet
  beta = order0 = order1 = -1;
}

Shapelets2D::Shapelets2D (int inorder0, int inorder1, double inbeta) {
  order0 = inorder0;
  order1 = inorder1;
  beta = inbeta;
  // define only one 1D shapelet object of highest order, that's enough
  S1D = Shapelets1D(GSL_MAX_INT(order0,order1),beta);
}

int Shapelets2D::getOrder (bool direction) {
  if (direction == 0) return order0;
  else return order1;
}

void Shapelets2D::setOrders (int inorder0, int inorder1) {
  order0 = inorder0;
  order1 = inorder1;
  // because of the empty default constructor we may not have defined
  // the Shapelets1D yet
   S1D = Shapelets1D(GSL_MAX_INT(order0,order1),beta);
}

double Shapelets2D::getBeta() {
  return beta;
}

void Shapelets2D::setBeta(double inbeta) {
  beta = inbeta;
  S1D.setBeta(beta);
}

double Shapelets2D::getThetaMin(int order0, int order1) {
  return beta/sqrt((double) GSL_MAX_INT(order0,order1)+1);
}

double Shapelets2D::getThetaMax(int order0, int order1) {
  return beta*sqrt((double) GSL_MAX_INT(order0,order1)+1);
}

double Shapelets2D::integrate(int order0, int order1) {
  return S1D.integrate(order0)*S1D.integrate(order1);
}

// integrate basis function within range
// see Paper III. eq. (82)
double Shapelets2D::integrate(int order0, int order1, double x0min, double x0max, double x1min,double x1max) {
  return S1D.integrate(order0,x0min,x0max)*S1D.integrate(order1,x1min,x1max);
}

double Shapelets2D::eval (int order0, int order1, Point2D& x) {
  return S1D.eval(order0,x(0)) * S1D.eval(order1,x(1));
}