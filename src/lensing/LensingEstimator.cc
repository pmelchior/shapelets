#include <lensing/LensingEstimator.h>
#include <math.h>
#include <gsl/gsl_math.h>

typedef complex<data_t> Complex;

LensingEstimator::LensingEstimator() {
}

complex<data_t> LensingEstimator::getShearMoments(ShapeletObject& so) {
  NumMatrix<data_t> Q = so.getShapelet2ndMoments();
  return getEllipticity(Q);
}

complex<data_t> LensingEstimator::getEllipticity(NumMatrix<data_t>& Q) {
  complex<data_t> I(0,1);
  complex<data_t> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<data_t> denom = Q11+Q22 + data_t(2)*sqrt(Q11*Q22-Q12*Q12);
  return (Q11 - Q22 + data_t(2)*I*Q12)/denom;
}

complex<data_t> LensingEstimator::getShearOrder(ShapeletObject& so, unsigned int n, complex<data_t> norm) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  if (n%2 ==0 && so.getNMax() >= n)
    return data_t(4./sqrt(data_t(n*(n+2)))) * f(n,2) / norm;
  else 
    return Complex(0,0);
}

complex<data_t> LensingEstimator::getNormShearOrder(ShapeletObjectList& ensemble, unsigned int n){
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const CoefficientVector<Complex>& f = (*iter)->getPolarCoeffs();
    if (n%2==0 && n >= 2 && (*iter)->getNMax() >= n+2) {
      norm += f(n-2,0) - f(n+2,0);
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getShearUnweighted(ShapeletObject& so, complex<data_t> norm) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex num(0,0),denom(0,0);
  for (int n = 0; n <= so.getNMax(); n+=2) {
    denom += data_t(sqrt((data_t) n+1)) * f(n,0);
    if (n>=2)
      num += data_t(sqrt(data_t(n*(n+2)))) * f(n,2);
  }
  return num / (norm * denom);
}

complex<data_t> LensingEstimator::getNormShearUnweighted(ShapeletObjectList& ensemble) {
  Complex norm(0,0), epsilon;
  int counter = 0;
  NumMatrix<data_t> Q;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    if ((*iter)->getNMax() >= 2) {
      Q = (*iter)->getShapelet2ndMoments();
      epsilon = getEllipticity(Q);
      norm += epsilon*epsilon;
      counter++;
    }
  }
  // ellipticity dispersion
  norm /= counter;
  // forming shear responsivity
  norm = Complex(2,0) - norm;
  return norm;
}

complex<data_t> LensingEstimator::getShearProfile(ShapeletObject& so, complex<data_t> norm) {
  Complex estimate(0,0);
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  for (int n = 2; n <= (int) so.getNMax() - 2; n+=2)
    estimate += data_t(sqrt((data_t)n*(n+2))) * (f(n-2,0) - f(n+2,0))  * f(n,2);
  return data_t(4) * estimate / norm;
}

complex<data_t> LensingEstimator::getNormShearProfile(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const CoefficientVector<Complex>& f = (*iter)->getPolarCoeffs();
    for (int n = 2; n <= (*iter)->getNMax() - 2; n+=2) {
      Complex summand = f(n-2,0) - f(n+2,0);
      norm += (data_t)n*(n+2) * (summand * summand);
    }
    // only with nmax > 4 there is a contribution to the norm
    if ((*iter)->getNMax() >= 4)
      counter++;
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getShearInvariant(ShapeletObject& so, NumMatrix<complex<data_t> >& norm) {
}

NumMatrix<complex<data_t> > LensingEstimator::getNormShearInvariant(ShapeletObjectList& ensemble) {
}

complex<data_t> LensingEstimator::getFlexionFOrder11(ShapeletObject& so,complex<data_t> norm) {
  return data_t(4./3) * so.getBeta() * (so.getPolarCoeffs())(1,1) / norm;
}

complex<data_t> LensingEstimator::getNormFlexionFOrder11(ShapeletObjectList& ensemble){
  Complex norm(0,0);
  int counter = 0;
  data_t beta2, R2;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const CoefficientVector<Complex>& f = (*iter)->getPolarCoeffs();
    if ((*iter)->getNMax() >= 4) {
      R2 = gsl_pow_2((*iter)->getShapeletRMSRadius());
      beta2 = gsl_pow_2((*iter)->getBeta());
      norm += (beta2-R2)*f(0,0) + R2*f(2,0) - beta2*f(4,0);
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getFlexionGOrder33(ShapeletObject& so, complex<data_t> norm) {
  return data_t((4./3)*sqrt(6.)/so.getBeta()) * (so.getPolarCoeffs())(3,3) / norm;
}

complex<data_t> LensingEstimator::getNormFlexionGOrder33(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const CoefficientVector<Complex>& f = (*iter)->getPolarCoeffs();
    if ((*iter)->getNMax() >= 6) {
      norm += f(0,0) + f(2,0) - f(4,0) - f(6,0);
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getWeightFlexionFProfile(ShapeletObject& so, int n) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex w(0,0);
  if (n >= 3 && n <= so.getNMax() - 3)
    w = data_t(sqrt(1.+n)*(n-1)) * (f(n-3,0) + f(n+1,0)) +
      data_t(sqrt(1.+n)*(n+3)) * (f(n-1,0) + f(n+3,0)) +
      4*data_t(gsl_pow_2(so.getShapeletRMSRadius())/gsl_pow_2(so.getBeta()) * sqrt(1.+n))*
      (f(n+1,0) + f(n-1,0));
  return w;
}

complex<data_t> LensingEstimator::getWeightFlexionGProfile(ShapeletObject& so, int n) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex w(0,0);
  if (n >= 3 && n <= so.getNMax() - 3)
    w = data_t(sqrt((data_t)(n-3)*(n-1)*(n+1))) * 
      (f(n-3,0) + f(n-1,0) - f(n+1,0) - f(n+3,0));
  return w;
}

complex<data_t> LensingEstimator::getFlexionFProfile(ShapeletObject& so, complex<data_t> norm){
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=1; n <= so.getNMax() - 1; n+=2)
    estimate += getWeightFlexionFProfile(so,n) * f(n,1);
  estimate *= (16 * M_SQRT2) / (3 * so.getBeta());
  return estimate / norm;
}

complex<data_t> LensingEstimator::getNormFlexionFProfile(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    for (int n=3; n <= (*iter)->getNMax() - 3; n+=2) {
      Complex weight = getWeightFlexionFProfile(*(*iter),n);
      norm += weight*weight;
    }
    // only with nmax >= 6 there is a contribution to the norm
    if ((*iter)->getNMax() >= 6) 
      counter++;
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getFlexionGProfile(ShapeletObject& so, complex<data_t> norm) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=3; n <= so.getNMax() - 3; n+=2)
    estimate += getWeightFlexionGProfile(so,n) * f(n,3);
  estimate *= (16 * M_SQRT2) / so.getBeta();
  return estimate / norm;
}

complex<data_t> LensingEstimator::getNormFlexionGProfile(ShapeletObjectList& ensemble){
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    for (int n=3; n <= (*iter)->getNMax() - 3; n+=2) {
      Complex weight = getWeightFlexionGProfile(*(*iter),n);
      norm += weight*weight;
    }
    // only with nmax >= 6 there is a contribution to the norm
    if ((*iter)->getNMax() >= 6) 
      counter++;
  }
  norm /= counter;
  return norm;
}

complex<data_t> LensingEstimator::getFlexionGDiagonal(ShapeletObject& so, complex<data_t> norm) {
  const CoefficientVector<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=3; n <= so.getNMax() - 3; n+=2)
    estimate += data_t(sqrt((data_t)(n-1)*(n+1)*(n+3))) * f(n,3);
  return data_t(2*M_SQRT2/so.getBeta()) * estimate / norm;
}

complex<data_t> LensingEstimator::getNormFlexionGDiagonal(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const CoefficientVector<Complex>& f = (*iter)->getPolarCoeffs();
    for (int n=0; n <= (*iter)->getNMax(); n+=2)
      norm += data_t(3*n*n - 4*n +15)*f(n,0);
  }
  // since every object in ensemble has at least f(0,0), all of them contribute to norm
  norm /= ensemble.size();
}
