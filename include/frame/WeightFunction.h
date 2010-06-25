#ifndef SHAPELENS_WEIGHTFUNCTION_H
#define SHAPELENS_WEIGHTFUNCTION_H

#include "../Typedef.h"
#include "Point.h"

namespace shapelens {
  class WeightFunction {
  public:
    virtual data_t operator()(const Point<data_t>& P) const = 0;
  };

  class FlatWeightFunction : public WeightFunction {
  public:
    FlatWeightFunction();
    virtual data_t operator()(const Point<data_t>& P) const;
  };

  class GaussianWeightFunction : public WeightFunction {
  public:
    GaussianWeightFunction(data_t scale, const Point<data_t>& centroid);
    virtual data_t operator()(const Point<data_t>& P) const;
    data_t getScale() const;
    void setScale(data_t scale);
    void setDerivative(int n);
    int getDerivative() const;
    void setCentroid(const Point<data_t>& centroid);
    const Point<data_t>& getCentroid() const;
  protected:
    Point<data_t> C;
    int n;
    data_t scale, sigma2;
    data_t (GaussianWeightFunction::*fptr) (const Point<data_t>&) const;
    data_t Gauss(data_t r) const;
    data_t Gauss(const Point<data_t>& P) const;
    data_t Gauss_(const Point<data_t>& P) const;
    data_t Gauss__(const Point<data_t>& P) const;
    data_t Gauss___(const Point<data_t>& P) const;
    data_t Gauss_2(const Point<data_t>& P) const;
    data_t Gauss__2(const Point<data_t>& P) const;
    data_t Gauss___2(const Point<data_t>& P) const;
  };

  inline data_t GaussianWeightFunction::operator()(const Point<data_t>& P) const {
    return (*this.*fptr)(P);
  }

} // end namespace
#endif
