#ifndef LENSINGESTIAMTOR_H
#define LENSINGESTIAMTOR_H

#include <Typedef.h>
#include <shapelets/ShapeletObject.h>
#include <shapelets/ShapeletObjectList.h>
#include <complex>
#include <vector>

class LensingEstimator {
 public:
  LensingEstimator();

  // various passive shear and flexion estimators
  complex<data_t> getShearMoments(ShapeletObject& so);

  complex<data_t> getShearOrder(ShapeletObject& so, unsigned int order, complex<data_t> norm);
  complex<data_t> getNormShearOrder(ShapeletObjectList& ensemble, unsigned int order);

  complex<data_t> getShearUnweighted(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormShearUnweighted(ShapeletObjectList& ensemble);
  
  complex<data_t> getShearProfile(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormShearProfile(ShapeletObjectList& ensemble);

  complex<data_t> getShearInvariant(ShapeletObject& so, NumMatrix<complex<data_t> >& norm);
  NumMatrix<complex<data_t> > getNormShearInvariant(ShapeletObjectList& ensemble);

  complex<data_t> getFlexionFOrder11(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormFlexionFOrder11(ShapeletObjectList& ensemble);

  complex<data_t> getFlexionGOrder33(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormFlexionGOrder33(ShapeletObjectList& ensemble);

  complex<data_t> getFlexionFHigherOrders();
  complex<data_t> getFlexionGHigherOrders();

  complex<data_t> getFlexionFProfile(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormFlexionFProfile(ShapeletObjectList& ensemble);

  complex<data_t> getFlexionGProfile(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormFlexionGProfile(ShapeletObjectList& ensemble);

  complex<data_t> getFlexionGDiagonal(ShapeletObject& so, complex<data_t> norm);
  complex<data_t> getNormFlexionGDiagonal(ShapeletObjectList& ensemble);

  // active shear and flexion estimation
  complex<data_t> getShearActive();
  void computeShearFlexionActive(complex<data_t>& gamma, complex<data_t>& F, complex<data_t>& G);
 private:
  ShapeletObjectList::iterator iter;
  complex<data_t> getEllipticity(NumMatrix<data_t>& Q);
  complex<data_t> getWeightFlexionFProfile(ShapeletObject& so, int n);
  complex<data_t> getWeightFlexionGProfile(ShapeletObject& so, int n);
};

#endif
