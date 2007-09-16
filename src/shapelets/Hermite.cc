#include <shapelets/Hermite.h>

using namespace boost::numeric::ublas;

// constructors
Hermite::Hermite () {
}

Hermite::Hermite (unsigned int order) {
  computed = 0;
  HermiteCoeffs = boost::numeric::ublas::triangular_matrix<double,lower>(order+1,order+1);
  HermiteCoeffs.clear();
  computeHermiteCoeffs(order);
}

void Hermite::computeHermiteCoeffs(unsigned int order) {
// returns vector of vectors of Hermite polynomial coefficients
// H[order][power]
// Hermite polynomial coefficients obey recurrence relation
// H_(l+1)(x) = 2*x*H_l(x) - 2*l*H_(l-1)(x)
  
  // define starting values for iteration
  // nothing computed yet
  if (computed == 0) {
    if (order >= 0) HermiteCoeffs(0,0) = 1;
    if (order >= 1) HermiteCoeffs(1,1) = 2;
    computed = 1;
  }
  while (computed<order) {
    for (int i=0; i<= computed; i+=1) {
      HermiteCoeffs(computed+1,i+1) = 2*HermiteCoeffs(computed,i);
      if (i<= computed-1) HermiteCoeffs(computed+1,i) -= 2*computed*HermiteCoeffs(computed-1,i);
      }
    computed++;
  }
}

int Hermite::getOrder () {
  return HermiteCoeffs.size1() -1;
}

void Hermite::printCoeffs (unsigned int order) {
  if (order > computed) setOrder(order);
  for (int i = 0; i < order+1; i+=1 ) {
    std::cout << HermiteCoeffs(order,i) << "\t";
  }  
  std::cout << std::endl;
}

void Hermite::setOrder (unsigned int order) {
  // only do something if, we need higher orders
  if (order > computed) {
    HermiteCoeffs.resize(order+1,order+1);
    computeHermiteCoeffs(order);
  }
}

double Hermite::eval (unsigned int order, double x) {
  if (order > computed) setOrder(order);
  double result = 0;
  for (int i = 0; i  < order+1; i+=1 )
    result +=  HermiteCoeffs(order,i) * gsl_pow_int(x,i);
  return result;
}

double Hermite::getCoefficient(unsigned int order, unsigned int power) {
  if (order > computed) setOrder(order);
  if (power > order) return 0;
  else return HermiteCoeffs(order,power);
}
