#include "../../include/utils/MathHelper.h"

namespace shapelens {
  /// see Wolfram MathWorld "Sample Variance"
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
