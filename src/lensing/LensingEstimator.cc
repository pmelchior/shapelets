#include <lensing/LensingEstimator.h>
#include <shapelets/MatrixManipulations.h>
#include <math.h>
#include <gsl/gsl_math.h>

typedef complex<double> Complex;

LensingEstimator::LensingEstimator() {
}

complex<double> LensingEstimator::getShearMoments(ShapeletObject& so) {
  NumMatrix<double> Q;
  so.getShapelet2ndMoments(Q);
  return 0.5*getEllipticity(Q);
}

complex<double> LensingEstimator::getEllipticity(NumMatrix<double>& Q) {
  complex<double> I(0,1);
  complex<double> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<double> denom = Q11+Q22;// + 2.*sqrt(Q11*Q22-Q12*Q12);
  return (Q11 - Q22 + 2.*I*Q12)/denom;
}

complex<double> LensingEstimator::getShearOrder(ShapeletObject& so, unsigned int n, complex<double> norm) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  if (so.getNMax() >= n + 2)
    return 4./double(sqrt((double) n*(n+2))) * f(n,mIndex(n,2)) / norm;
  else 
    return Complex(0,0);
}

complex<double> LensingEstimator::getNormShearOrder(ShapeletObjectList& ensemble, unsigned int n){
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    if (n >= 2 && (*iter)->getNMax() >= n+2) {
      norm += f(n-2,mIndex(n-2,0)) - f(n+2,mIndex(n+2,0));
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<double> LensingEstimator::getShearUnweighted(ShapeletObject& so, complex<double> norm) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n = 2; n <= (int) so.getNMax() - 2; n+=2)
    estimate += double(sqrt((double) n*(n+2))) * f(n,mIndex(n,2));
  return estimate / norm;
}

complex<double> LensingEstimator::getNormShearUnweighted(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  double ellipticity2 = 0;
  int counter = 0;
  NumMatrix<double> Q(2,2);
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    //std::cout << (*iter)->getHistory() << std::endl;
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    for (int n = 0; n <= (*iter)->getNMax(); n+=2)
      norm += double(sqrt((double) n+1)) * f(n,mIndex(n,0));
    if ((*iter)->getNMax() >= 2) {
      (*iter)->getShapelet2ndMoments(Q);
      counter++;
      // TODO: is this the complex ellipticity or its modulus?
      ellipticity2 += gsl_pow_2(abs(getEllipticity(Q)));
    }
  }
  // ellipticity dispersion
  ellipticity2 /= counter;
  // mean of norm
  norm /= counter;
  // including factor 2 and shear responsivity
  norm *= (2-ellipticity2);
  return norm;
}

complex<double> LensingEstimator::getShearProfile(ShapeletObject& so, complex<double> norm) {
  Complex estimate(0,0);
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  for (int n = 2; n <= (int) so.getNMax() - 2; n+=2) {
    Complex summand = f(n-2,mIndex(n-2,0)) - f(n+2,mIndex(n+2,0));
    estimate += double(sqrt((double)n*(n+2))) * summand * f(n,mIndex(n,2));
  }
  return 4. * estimate / norm;
}

complex<double> LensingEstimator::getNormShearProfile(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    for (int n = 2; n <= (*iter)->getNMax() - 2; n+=2) {
      Complex summand = f(n-2,mIndex(n-2,0)) - f(n+2,mIndex(n+2,0));
      norm += (double)n*(n+2) * (summand * summand);
    }
    // only with nmax > 4 there is a contribution to the norm
    if ((*iter)->getNMax() >= 4)
      counter++;
  }
  norm /= counter;
  return norm;
}

complex<double> LensingEstimator::getShearInvariant(ShapeletObject& so, NumMatrix<complex<double> >& norm) {
}

NumMatrix<complex<double> > LensingEstimator::getNormShearInvariant(ShapeletObjectList& ensemble) {
}

complex<double> LensingEstimator::getFlexionFOrder11(ShapeletObject& so,complex<double> norm) {
  return (4./3) * so.getBeta() * (so.getPolarCoeffs())(1,mIndex(1,1)) / norm;
}

complex<double> LensingEstimator::getNormFlexionFOrder11(ShapeletObjectList& ensemble){
  Complex norm(0,0);
  int counter = 0;
  double beta2, R2;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    if ((*iter)->getNMax() >= 4) {
      R2 = gsl_pow_2((*iter)->getShapeletRMSRadius());
      beta2 = gsl_pow_2((*iter)->getBeta());
      norm += (beta2-R2)*f(0,mIndex(0,0)) + R2*f(2,mIndex(2,0)) - beta2*f(4,mIndex(4,0));
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<double> LensingEstimator::getFlexionGOrder33(ShapeletObject& so, complex<double> norm) {
  return ((4./3)*double(sqrt(6.))/so.getBeta()) * (so.getPolarCoeffs())(3,mIndex(3,3)) / norm;
}

complex<double> LensingEstimator::getNormFlexionGOrder33(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  int counter = 0;
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    if ((*iter)->getNMax() >= 6) {
      norm += f(0,mIndex(0,0)) + f(2,mIndex(2,0)) - f(4,mIndex(4,0)) - f(6,mIndex(6,0));
      counter++;
    }
  }
  norm /= counter;
  return norm;
}

complex<double> LensingEstimator::getWeightFlexionFProfile(ShapeletObject& so, int n) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex w(0,0);
  if (n >= 3 && n <= so.getNMax() - 3)
    w = double(sqrt(1.+n)*(n-1)) * (f(n-3,mIndex(n-3,0)) + f(n+1,mIndex(n+1,0))) +
      double(sqrt(1.+n)*(n+3)) * (f(n-1,mIndex(n-1,0)) + f(n+3,mIndex(n+3,0))) +
      4*double(gsl_pow_2(so.getShapeletRMSRadius())/gsl_pow_2(so.getBeta()) * sqrt(1.+n))*
      (f(n+1,mIndex(n+1,0)) + f(n-1,mIndex(n-1,0)));
  return w;
}

complex<double> LensingEstimator::getWeightFlexionGProfile(ShapeletObject& so, int n) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex w(0,0);
  if (n >= 3 && n <= so.getNMax() - 3)
    w = double(sqrt((double)(n-3)*(n-1)*(n+1))) * 
      (f(n-3,mIndex(n-3,0)) + f(n-1,mIndex(n-1,0)) - f(n+1,mIndex(n+1,0)) - f(n+3,mIndex(n+3,0)));
  return w;
}

complex<double> LensingEstimator::getFlexionFProfile(ShapeletObject& so, complex<double> norm){
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=1; n <= so.getNMax() - 1; n+=2)
    estimate += getWeightFlexionFProfile(so,n) * f(n,mIndex(n,1));
  estimate *= (16 * M_SQRT2) / (3 * so.getBeta());
  return estimate / norm;
}

complex<double> LensingEstimator::getNormFlexionFProfile(ShapeletObjectList& ensemble) {
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

complex<double> LensingEstimator::getFlexionGProfile(ShapeletObject& so, complex<double> norm) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=3; n <= so.getNMax() - 3; n+=2)
    estimate += getWeightFlexionGProfile(so,n) * f(n,mIndex(n,3));
  estimate *= (16 * M_SQRT2) / so.getBeta();
  return estimate / norm;
}

complex<double> LensingEstimator::getNormFlexionGProfile(ShapeletObjectList& ensemble){
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

complex<double> LensingEstimator::getFlexionGDiagonal(ShapeletObject& so, complex<double> norm) {
  const NumMatrix<Complex>& f = so.getPolarCoeffs();
  Complex estimate(0,0);
  for (int n=3; n <= so.getNMax() - 3; n+=2)
    estimate += double(sqrt((double)(n-1)*(n+1)*(n+3))) * f(n,mIndex(n,3));
  return (2*M_SQRT2/so.getBeta()) * estimate / norm;
}

complex<double> LensingEstimator::getNormFlexionGDiagonal(ShapeletObjectList& ensemble) {
  Complex norm(0,0);
  for (iter = ensemble.begin(); iter != ensemble.end(); iter++) {
    const NumMatrix<Complex>& f = (*iter)->getPolarCoeffs();
    for (int n=0; n <= (*iter)->getNMax(); n+=2)
      norm += double(3*n*n - 4*n +15)*f(n,mIndex(n,0));
  }
  // since every object in ensemble has at least f(0,0), all of them contribute to norm
  norm /= ensemble.size();
}
