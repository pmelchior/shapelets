#include "../../include/lensing/DEIMOSForward.h"
#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/MathHelper.h"

namespace shapelens {

  DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const MultiExposureMoments& mepsf_, int N, int C, data_t flux, data_t width_) :
    meo(meo_), mepsf(mepsf_), width(width_), K(meo_.size()) {
    DEIMOS::N = N;
    DEIMOS::C = C;
    DEIMOS::D = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+C+1)*(N+C+2))/2);
    DEIMOS::S = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
    DEIMOS::mo.setOrder(DEIMOS::N);

    // set up containers
    {
      Moments tmp(N);
      NumMatrix<data_t> P(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
      for (int k = 0; k < K; k++) {
	mem.push_back(tmp);
	meP.push_back(P);
	meD.push_back(*this);
      }
    }
    
    // set noise image from each exposure
    for (int k = 0; k < K; k++)
      meD[k].setNoiseImage(meo[k]);

    // initialize 0th and 2nd moments 
    // with circular Gaussian with given flux & scale
    mo(0,0) = flux;
    // s = sqrt(trQ/F) = sqrt(2*Q_ii/F)
    mo(0,2) = mo(2,0) = width*width*mo(0,0)/2;

    // Minimize chi^2
    // FIXME: need convergence criterium
    for (int t = 0; t < 10; t++) {
      //std::cout << mo << std::endl;
      computeMomentsFromGuess();

      // compute chi^2 and best-fit moments
      data_t chi2 = 0;
      NumVector<data_t> diff(mo.size());
      mo.clear();
      S.clear();
      for (int k = 0; k < K; k++) {
	DEIMOS& d = meD[k];
	diff = mem[k];
	diff -= d.mo;
	
	NumMatrix<data_t> S_1 = d.S.invert();
	chi2 += diff * (S_1 * (NumVector<data_t>) diff);
	NumMatrix<data_t> X = meP[k].transpose() * S_1;
	mo += X * (NumVector<data_t>) d.mo;
	S += X*meP[k];
      }
      S = S.invert();
      mo = S * (NumVector<data_t>)mo;
     
      // non-sensical ellipticity check
      data_t tiny = 1e-4;
      mo(0,0) = std::max(tiny, mo(0,0));
      mo(0,2) = std::max(tiny*flux, mo(0,2));
      mo(2,0) = std::max(tiny*flux, mo(2,0));
      if (mo(1,1) > 0)
	mo(1,1) = std::min(mo(1,1), sqrt(mo(0,2)*mo(2,0)));
      else
	mo(1,1) = std::max(mo(1,1), -sqrt(mo(0,2)*mo(2,0)));

      // FIXME: should we update the centroid?
      // globally or individually?
      // pre-convolved or convolved?
      // Should be tied in the size determination of the weight function?

//       // update centroid
//       for (int i = 0; i < meo.size(); i++) {
// 	Object& obj = const_cast<Object&>(meo[i]);
// 	obj.centroid(0) += mo0(1,0)/mem[i](0,0);
// 	obj.centroid(1) += diff(0,1)/mem[i](0,0);
//       }

    }
    // set DEIMOS parameters to best fit
    // FIXME: what to do with SN, scale, centroid etc.
  }

  void DEIMOSForward::computeMomentsFromGuess() {
    // 1) convolve guess with psf moments
    // 2) measure deweighted moments from each exposure
    //    with weight function shape based on guess of convolved moments
    // 3) get the convolved moment errors
    // 4) compute the contribution to chi^2
    for (int k = 0; k < K; k++) {
      convolveExposure(k);
      Moments& mo_c = mem[k];
      DEIMOS& d = meD[k];
      // set ellipticities and sizes for weight functions in each exposure
      // FIXME: how to set the width (think varying PSF FWHM in exposures) 
      d.scale = width;//sqrt((mo0c(0,2) + mo0c(2,0))/mo0c(0,0));
      d.eps = shapelens::epsilon(mo_c);
      data_t abs_eps = abs(d.eps);
      // FIXME: need individual WCS scale_factors here!
      d.scale *= d.scale_factor/sqrt(1 + abs_eps*abs_eps - 2*abs_eps);
      DEIMOS::DEIMOSWeightFunction w(d.scale, meo[k].centroid, d.eps);
      Moments mo_w(meo[k], w, N+C);
      d.deweight(mo_w);
      d.computeCovariances(mo_w);
    }
  }

  void DEIMOSForward::convolveExposure(unsigned int k) {
    const Moments& p = mepsf[k];
    NumMatrix<data_t>& P = meP[k];

    // matrix representation of convolution eq. 9
    for (int n = 0; n <= mo.getOrder(); n++) {
      for (int i = 0; i <= n; i++) {
	int  j = n-i;
	int n = mo.getIndex(i,j);  // convolved moment index
	for (int k = 0; k <= i; k++) {
	  for (int l = 0; l <= j; l++) {
	    int m = mo.getIndex(k,l); // unconvolved moment index
	    P(n,m) = binomial(i,k) * binomial(j,l) * p(i-k,j-l);
	  }
	}
      }
    }
    P.gemv(mo, mem[k]);
  }

}
