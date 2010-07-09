#ifndef SHAPELENS_MATHHELPER_H
#define SHAPELENS_MATHHELPER_H

#include "../Typedef.h"

namespace shapelens {

  /// Templated replacement of \p gsl_pow_int.
  template <class T> 
    T pow_int(const T& x, int n) {
    T x2,x3,x4;
    switch (n) {
    case 0: return T(1);
    case 1: return x;
    case 2: return x*x;
    case 3: return x*x*x;
    case 4: x2 = x*x; return x2*x2;
    case 5: x2 = x*x; return x2*x2*x;
    case 6: x3 = x*x*x; return x3*x3;
    case 7: x3 = x*x*x; return x3*x3*x;
    case 8: x2 = x*x; x4 = x2*x2; return x4*x4;
    case 9: x3 = x*x*x; return x3*x3*x3;
    case 10:x2 = x*x; x4 = x2*x2; return x4*x4*x2;
    default: x2 = T(1); // dummy variable for recursive calls
      while (n > 10) {
	x2 *= pow_int(x,10);
	n-=10;
      }
      return x2*pow_int(x,n); 
    }
  }

  /// Computes the mean and variance of a series of \p x_i up to
  /// the index \p i.
  void mean_variance(const data_t& x_i, int i, data_t& mean, data_t& var) {
    if (i==0) {
      mean = x_i;
      var = 0;
    } else {
      data_t mean_j = mean;
      mean = ((mean*i) + x_i)/(i+1);
      var = (1-1./i)*var + (i+1)*pow_int(mean - mean_j,2);
    }
  }
}


#endif
