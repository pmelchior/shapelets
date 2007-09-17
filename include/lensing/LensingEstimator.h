#ifndef LENSINGESTIAMTOR_H
#define LENSINGESTIAMTOR_H

#include <shapelets/ShapeletObject.h>
#include <shapelets/ShapeletObjectList.h>
#include <complex>
#include <vector>

class LensingEstimator {
 public:
  LensingEstimator();

  // various passive shear and flexion estimators
  complex<double> getShearMoments(ShapeletObject& so);

  complex<double> getShearOrder(ShapeletObject& so, unsigned int order, complex<double> norm);
  complex<double> getNormShearOrder(ShapeletObjectList& ensemble, unsigned int order);

  complex<double> getShearUnweighted(ShapeletObject& so, complex<double> norm);
  complex<double> getNormShearUnweighted(ShapeletObjectList& ensemble);
  
  complex<double> getShearProfile(ShapeletObject& so, complex<double> norm);
  complex<double> getNormShearProfile(ShapeletObjectList& ensemble);

  complex<double> getShearInvariant(ShapeletObject& so, NumMatrix<complex<double> >& norm);
  NumMatrix<complex<double> > getNormShearInvariant(ShapeletObjectList& ensemble);

  complex<double> getFlexionFOrder11(ShapeletObject& so, complex<double> norm);
  complex<double> getNormFlexionFOrder11(ShapeletObjectList& ensemble);

  complex<double> getFlexionGOrder33(ShapeletObject& so, complex<double> norm);
  complex<double> getNormFlexionGOrder33(ShapeletObjectList& ensemble);

  complex<double> getFlexionFHigherOrders();
  complex<double> getFlexionGHigherOrders();

  complex<double> getFlexionFProfile(ShapeletObject& so, complex<double> norm);
  complex<double> getNormFlexionFProfile(ShapeletObjectList& ensemble);

  complex<double> getFlexionGProfile(ShapeletObject& so, complex<double> norm);
  complex<double> getNormFlexionGProfile(ShapeletObjectList& ensemble);

  complex<double> getFlexionGDiagonal(ShapeletObject& so, complex<double> norm);
  complex<double> getNormFlexionGDiagonal(ShapeletObjectList& ensemble);

  // active shear and flexion estimation
  complex<double> getShearActive();
  void computeShearFlexionActive(complex<double>& gamma, complex<double>& F, complex<float_t>& G);
 private:
  ShapeletObjectList::iterator iter;
  complex<double> getEllipticity(NumMatrix<double>& Q);
  complex<double> getWeightFlexionFProfile(ShapeletObject& so, int n);
  complex<double> getWeightFlexionGProfile(ShapeletObject& so, int n);
};

#endif
