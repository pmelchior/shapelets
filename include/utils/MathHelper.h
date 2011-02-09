#ifndef SHAPELENS_MATHHELPER_H
#define SHAPELENS_MATHHELPER_H

#include "../Typedef.h"
#include <numla/NumVector.h>
#include <set>

namespace shapelens {

  /// Templated replacement of \p gsl_pow_int.
  template <class T> 
  inline  T pow_int(const T& x, int n) {
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

  /// Perform \f$\kappa-\sigma\f$-clipping to calculate mean and
  /// and standard deviation of \p data.
  /// The algorithm iteratively determines median \f$m\f$ and 
  /// standard deviation \f$\sigma\f$ of all pixels within 
  /// \f$[m-3\sigma,m+3\sigma]\f$. If a new iteration works on the same
  /// set of pixels than the previous one, the algorithm terminates.\n
  std::pair<data_t, data_t> kappa_sigma_clip(const NumVector<data_t>& data, unsigned long maxLength = 1000000);

  /// Compute \f$\frac{1}{2}^n\f$ quantiles of \p data.
  std::set<double> quantiles(const NumVector<data_t>& data, unsigned int n);
}


#endif
