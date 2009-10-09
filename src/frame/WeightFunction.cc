#include "../../include/frame/WeightFunction.h"
#include <stdexcept>
#include <gsl/gsl_math.h>

namespace shapelens {

  WeightFunction::WeightFunction() : n(0), type(0) {
    fptr = &WeightFunction::Flat;
  }
  WeightFunction::WeightFunction(data_t scale_, const Point<data_t>& centroid_):
    scale(scale_), C(centroid_), n(0), type(1) {
    sigma2 = scale*scale;
    fptr = &WeightFunction::Gauss;
  }
  
  data_t  WeightFunction::Gauss(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r);
  }
  data_t  WeightFunction::Gauss_(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(-r/sigma2);
  }
  data_t  WeightFunction::Gauss__(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(gsl_pow_2(r/sigma2) - 1./sigma2);
  }
  data_t  WeightFunction::Gauss___(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)*(gsl_pow_3(-r/sigma2) + 3*r/(sigma2*sigma2));
  }
  data_t  WeightFunction::Gauss_2(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return -Gauss(r)/sigma2;
  }
  data_t  WeightFunction::Gauss__2(const Point<data_t>& P) const {
    data_t r = sqrt(gsl_pow_2(P(0)-C(0)) + gsl_pow_2(P(1)-C(1)));
    return Gauss(r)/gsl_pow_2(2*sigma2);
  }
  data_t  WeightFunction::Gauss(data_t r) const {
    return exp(-r*r/(2*sigma2));
  }
  data_t  WeightFunction::Flat(const Point<data_t>& P) const {
    return 1;
  }
  data_t  WeightFunction::Flat_(const Point<data_t>& P) const {
    return 0;
  }

  void WeightFunction::setType(int type_) {
    if (type_ != type) {
      type = type_;
      setDerivative(n);
    }
  }
  int WeightFunction::getType() const {
    return type;
  }

  void WeightFunction::setDerivative(int n_) {
    n = n_;
    switch (type) {
    case 1: 
      switch (n) {
      case 0: fptr = &WeightFunction::Gauss; break;
      case 1: fptr = &WeightFunction::Gauss_; break;
      case -1: fptr = &WeightFunction::Gauss_2; break;
      case 2: fptr = &WeightFunction::Gauss__; break;
      case -2: fptr = &WeightFunction::Gauss__2; break;
      case 3: fptr = &WeightFunction::Gauss___; break;
      default: throw std::invalid_argument("WeightFunction: derivative of Gaussian invalid");
      }
    default:
      if (n==0)
	fptr = &WeightFunction::Flat;
      else
	fptr = &WeightFunction::Flat_;
    }
  }
  int WeightFunction::getDerivative() const {
    return n;
  }

  data_t WeightFunction::getScale() const {
    return scale;
  }
  void WeightFunction::setScale(data_t scale_) {
    scale = scale_;
    sigma2 = scale*scale;
  }
}
