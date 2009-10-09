#ifndef SHAPELENS_WEIGHTFUNCTION_H
#define SHAPELENS_WEIGHTFUNCTION_H

#include "../Typedef.h"
#include "Point.h"

namespace shapelens {
  class WeightFunction {
  public:
    WeightFunction();
    WeightFunction(data_t scale, const Point<data_t>& centroid);
    void setType(int type);
    int getType() const;
    void setCentroid(const Point<data_t>& centroid);
    const Point<data_t>& getCentroid() const;
    inline data_t operator()(const Point<data_t>& P) const;
    data_t getScale() const;
    void setScale(data_t scale);
    void setDerivative(int n);
    int getDerivative() const;
  private:
    int n, type;
    data_t scale, sigma2;
    Point<data_t> C;
    data_t (WeightFunction::*fptr) (const Point<data_t>&) const;
    data_t Gauss(data_t r) const;
    data_t Gauss(const Point<data_t>& P) const;
    data_t Gauss_(const Point<data_t>& P) const;
    data_t Gauss__(const Point<data_t>& P) const;
    data_t Gauss___(const Point<data_t>& P) const;
    data_t Gauss_2(const Point<data_t>& P) const;
    data_t Gauss__2(const Point<data_t>& P) const;
    data_t Flat(const Point<data_t>& P) const;
    data_t Flat_(const Point<data_t>& P) const;
  };
} // end namespace
#endif
