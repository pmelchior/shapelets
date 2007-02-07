#include <Shapelets1D.h>
// for factorial function
#include <gsl/gsl_sf.h>
// for computing the error function in integrate(order,xmin,xmax)
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>

Shapelets1D::Shapelets1D() {
}

Shapelets1D::Shapelets1D(int order, double inbeta) {
  H = Hermite(order);
  beta = inbeta;
  // defined for fast computation
  sqrt_beta = sqrt(beta);
}

int Shapelets1D::getOrder () {
  return H.getOrder();
}

void Shapelets1D::setOrder (int order) {
  H.setOrder(order);
}

double Shapelets1D::getBeta() {
  return beta;
}

void Shapelets1D::setBeta(double inbeta) {
  beta = inbeta;
  // defined for fast computation
  sqrt_beta = sqrt(beta);
}

double Shapelets1D::getThetaMin(int order) {
  return beta*1./sqrt(order + 0.5);
}

double Shapelets1D::getThetaMax(int order) {
  return beta*sqrt(order + 0.5);
}

double Shapelets1D::eval (int order, double x) {
  double x_scaled = x/beta;
  return 1./(sqrt_beta * sqrt(M_SQRTPI*gsl_pow_int(2,order)*gsl_sf_fact(order))) * 
  H.eval(order,x_scaled) *
  exp(-x_scaled*x_scaled/2); 
}

// integral over basis function
// see Paper I, eq (17) 
double Shapelets1D::integrate(int order) {
  double result;
  if (order%2 != 0) result = 0;
  else result = sqrt_beta * sqrt(gsl_pow_int(2,1-order)* M_SQRTPI * gsl_sf_fact(order))
       / gsl_sf_fact(order/2);
  return result;
}

// integrate basis function within range xmin-xmax
// see Paper III, eq. (78) - (81)
double Shapelets1D::integrate(int order, double xmin, double xmax) {
  double result;
  if (order == 0) 
    result = (gsl_sf_erf(xmax/(M_SQRT2*beta)) - gsl_sf_erf(xmin/(M_SQRT2*beta))) * sqrt_beta *sqrt(M_SQRTPI/2);
  else if (order == 1) 
    result = -sqrt_beta * M_SQRT2 * (eval(0,xmax) - eval(0,xmin));
  else if (order > 1) {
    result = -beta * sqrt((double) 2/order) * 
      (eval(order-1,xmax) - eval(order-1,xmin)) +
      sqrt((double) (order-1)/order) * integrate(order-2,xmin,xmax);
  }
  return result;
}
   
