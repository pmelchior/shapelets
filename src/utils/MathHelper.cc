#include "../../include/utils/MathHelper.h"
#include <numla/NumVectorMasked.h>
#include <list>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>

namespace shapelens {
  
  std::pair<data_t, data_t> kappa_sigma_clip(const NumVector<data_t>& data, unsigned long maxLength) {
    data_t mean, sigma;
    // since sorting is performed in place, we need a copy of the vector
    // if the vector is too large we take a subset of 1 million entries
    // randomly chosen
    NumVector<data_t> D;
    unsigned long datamax = data.size(), jmax;
    if (datamax > maxLength) {
      jmax = maxLength;
      gsl_rng_env_setup();
      const gsl_rng_type * R = gsl_rng_default;
      gsl_rng * r = gsl_rng_alloc (R);
      D.resize(jmax);
      for (unsigned long i=0; i < jmax; i++) {
	unsigned long j = (unsigned long) floor(datamax * gsl_rng_uniform (r));
	D(i) = data(j);
      }
      gsl_rng_free (r);
    } else {
      jmax = datamax;
      D = data;
    }

    gsl_sort(D.c_array(),1,jmax);
    sigma = gsl_stats_sd(D.c_array(),1,jmax);
    mean = gsl_stats_median_from_sorted_data (D.c_array(),1,jmax);

    // sigma clipping here
    unsigned int j = jmax;
    while (j >= 1) {
      j=0;
      for (int i = 0; i < jmax; i++ ) {
	// if data is 3 sigma arround iterative median, keep it
	if (D(i) > mean-3*sigma && D(i) < mean+3*sigma)  {
	  D(j) = D(i);
	  j++;
	}
      }
      // next time only work on the first jmax = j, all others are not sorted
      if (j >= 1) {
	gsl_sort(D.c_array(), 1, j);
 	mean = gsl_stats_median_from_sorted_data (D.c_array(),1,j);
 	sigma = gsl_stats_sd(D.c_array(),1,j);
	// no change in the selected pixels: converged
	if (jmax == j)
	  break;
	else
	  jmax = j;
      }
      else
	mean = sigma = 0;
    }
    return std::pair<data_t, data_t> (mean, sigma);
  }

  double maskedMedian(const NumVector<data_t>& data, NumVector<bool>& mask, double min, double max) {
    int size = 0;
    for (int i=0; i < data.size(); i++) {
      if (data(i) >= min && data(i) <= max)
	mask(i) = false;
      else
	mask(i) = true;
    }
    NumVectorMasked<data_t> data_masked(data,mask);
    return data_masked.median();
  }

  std::set<double> quantiles(const NumVector<data_t>& data, unsigned int n) {
    // use set because it's sorted
    std::set<double> q;
    std::set<double>::iterator siter;
    // for saving current median list
    std::list<double> tmp;
    std::list<double>::iterator iter;
    // add min and max for easier handling later
    q.insert(data.min());
    q.insert(data.max());
    if (n >= 1) {
      // for n = 1 we only need the median of the entire data
      q.insert(data.median());
      // otherwise: we compute the median of data with values between
      // min..median and median..max
      // -> iteratively split the data set into series of median
      // we need mask to exclude entries outside the considered bounds
      if (n > 1) {
	NumVector<bool> mask(data.size());
	for (int l=1; l < n; l++) {
	  tmp.clear();
	  // copy current list of quantiles
	  // because the insert below would change the iterator
	  for(siter = q.begin(); siter != q.end(); siter++)
	    tmp.push_front(*siter);
	  iter = tmp.begin();
	  while (iter != --(tmp.end()))
	    q.insert(maskedMedian(data, mask,*iter,*(iter++))); // increments iter
	}
      }
    }
    return q;
  }
}
