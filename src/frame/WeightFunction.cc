#include "../../include/frame/WeightFunction.h"
#include <stdexcept>
#include <gsl/gsl_math.h>

namespace shapelens {

  FlatWeightFunction::FlatWeightFunction() {
  }
  
  data_t FlatWeightFunction::operator()(const Point<data_t>& P) const {
    return 1;
  }

  GaussianWeightFunction::GaussianWeightFunction(data_t scale_, const Point<data_t>& centroid_):
    scale(scale_), n(0), sigma2(scale_*scale_) {
    C = centroid_;
    fptr = &GaussianWeightFunction::Gauss;
  }

  void GaussianWeightFunction::setCentroid(const Point<data_t>& centroid) {
    C = centroid;
  }
  const Point<data_t>& GaussianWeightFunction::getCentroid() const {
    return C;
  }  

  data_t  GaussianWeightFunction::Gauss(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r);
  }
  data_t  GaussianWeightFunction::Gauss_(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(-r/sigma2);
  }
  data_t  GaussianWeightFunction::Gauss__(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(gsl_pow_2(r/sigma2) - 1./sigma2);
  }
  data_t  GaussianWeightFunction::Gauss___(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(gsl_pow_3(-r/sigma2) + 3*r/(sigma2*sigma2));
  }
  data_t  GaussianWeightFunction::Gauss_2(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return -Gauss(r)/(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss__2(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)/gsl_pow_2(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss___2(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return -Gauss(r)/gsl_pow_3(2*sigma2);
  }
  data_t  GaussianWeightFunction::Gauss(data_t r) const {
    return exp(-r*r/(2*sigma2));
  }

  void GaussianWeightFunction::setDerivative(int n_) {
    n = n_;
    switch (n) {
      case 0: fptr = &GaussianWeightFunction::Gauss; break;
      case 1: fptr = &GaussianWeightFunction::Gauss_; break;
      case -1: fptr = &GaussianWeightFunction::Gauss_2; break;
      case 2: fptr = &GaussianWeightFunction::Gauss__; break;
      case -2: fptr = &GaussianWeightFunction::Gauss__2; break;
      case 3: fptr = &GaussianWeightFunction::Gauss___; break;
      case -3: fptr = &GaussianWeightFunction::Gauss___2; break;
      default: throw std::invalid_argument("WeightFunction: derivative of Gaussian invalid");
    }
  }
  int GaussianWeightFunction::getDerivative() const {
    return n;
  }

  data_t GaussianWeightFunction::getScale() const {
    return scale;
  }
  void GaussianWeightFunction::setScale(data_t scale_) {
    scale = scale_;
    sigma2 = scale*scale;
  }
}
