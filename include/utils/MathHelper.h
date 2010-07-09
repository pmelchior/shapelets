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

  /// Computes sample mean and variance of a series of \p x_i up to
  /// the index \p i.
  /// \code
  /// NumVector<data_t> v(N);
  /// data_t mean, var;
  /// for (int i=0; i < v.size(); i++) {
  ///   v(i) = some_function();
  ///   mean_variance(v(i), i, mean, var);
  /// }
  /// \endcode
  void mean_variance(const data_t& x_i, int i, data_t& mean, data_t& var);

}


#endif
