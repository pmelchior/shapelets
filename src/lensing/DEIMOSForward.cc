#include "../../include/lensing/DEIMOSForward.h"
#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/MathHelper.h"

namespace shapelens {

  inline data_t sign(data_t x) {
    if (x < 0)
      return -1;
    else
      return 1;
  }

  DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const MultiExposureMoments& mepsf_, int N, int C, data_t flux, data_t width_) :
    meo(meo_), mepsf(mepsf_), width(width_) {
    DEIMOS::N = N;
    DEIMOS::C = C;
    DEIMOS::D = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+C+1)*(N+C+2))/2);
    DEIMOS::S = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
    
    // set up containers
    {
      Moments tmp(N);
      NumMatrix<data_t> P(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
      for (int i = 0; i < meo.size(); i++) {
	mem.push_back(tmp);
	mem_.push_back(tmp);
	meP.push_back(P);
	meS.push_back(P);
      }
    }

    // initialize 0th and 2nd moments 
    // with circular Gaussian with given flux & scale
    mo0 = Moments(N);
    mo0(0,0) = flux;
    // s = sqrt(trQ/F) = sqrt(2*Q_ii/F)
    mo0(0,2) = mo0(2,0) = width*width*mo0(0,0)/2;

    // Minimize chi^2
    // FIXME: need convergence criterium
    for (int t = 0; t < 10; t++) {
      computeMomentsFromGuess();
      // compute chi^2 and best-fit moments
      data_t chi2 = 0;
      NumVector<data_t> diff(mo.size());
      mo0.clear();
      DEIMOS::S.clear();
      for (int k = 0; k < meo.size(); k++) {
	diff = mem[k];
	diff -= mem_[k];
	NumMatrix<data_t> S_1 = meS[k].invert();
	chi2 += diff * (S_1 * (NumVector<data_t>) diff);
	NumMatrix<data_t> X = meP[k].transpose() * S_1;
	NumMatrix<data_t> C = (X*meP[k]).invert();
	DEIMOS::S += C;
	mo0 += C * X * (NumVector<data_t>) mem_[k];
	//std::cout << k << "\t" << diff << "\t" << chi2 << std::endl;
      }
      //std::cout << std::endl;
      mo0 /= meo.size();
      DEIMOS::S /= meo.size()*meo.size();
      //std::cout << DEIMOS::S << std::endl;
     
      // non-sensical ellipticity check
      data_t tiny = 1e-4;
      mo0(0,0) = std::max(tiny, mo0(0,0));
      mo0(0,2) = std::max(tiny*flux, mo0(0,2));
      mo0(2,0) = std::max(tiny*flux, mo0(2,0));
      if (mo0(1,1) > 0)
	mo0(1,1) = std::min(mo0(1,1), sqrt(mo0(0,2)*mo0(2,0)));
      else
	mo0(1,1) = std::max(mo0(1,1), -sqrt(mo0(0,2)*mo0(2,0)));

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
    mo = mo0;

  }

  void DEIMOSForward::computeMomentsFromGuess() {
    // 1) convolve guess with psf moments
    // 2) measure deweighted moments from each exposure
    //    with weight function shape based on guess of convolved moments
    // 3) get the convolved moment errors
    // 4) compute the contribution to chi^2
    DEIMOS::scale_factor = 1; // FIXME: WCS info needed here
    DEIMOS::mo.setOrder(DEIMOS::N);
    for (int i = 0; i < meo.size(); i++) {
      convolveExposure(i);
      Moments& mo0c = mem[i];
         
      // set ellipticities and sizes for weight functions in each exposure
      // FIXME: how to set the width (think varying PSF FWHM in exposures) 
      DEIMOS::scale = width;//sqrt((mo0c(0,2) + mo0c(2,0))/mo0c(0,0));
      DEIMOS::eps = shapelens::epsilon(mo0c);
      data_t abs_eps = abs(DEIMOS::eps);
      DEIMOS::scale *= DEIMOS::scale_factor/sqrt(1 + abs_eps*abs_eps - 2*abs_eps);
      //std::cout << scale << "\t" << meo[i].centroid << "\t" <<  eps << "\t" << shapelens::epsilon(mo0) << std::endl;
      DEIMOS::DEIMOSWeightFunction w(scale, meo[i].centroid, eps);
      Moments mo_w(meo[i], w, N+C);
      DEIMOS::deweight(mo_w);
      mem_[i] = DEIMOS::mo;
      // std::cout << "# " << mo << std::endl;

      // compute moment errors
      DEIMOS::setNoiseImage(meo[i]); // FIXME: large overhead
      DEIMOS::computeCovariances(mo_w);
      meS[i] = DEIMOS::S;
    }
  }

  void DEIMOSForward::convolveExposure(unsigned int k) {
    const Moments& p = mepsf[k];
    NumMatrix<data_t>& P = meP[k];

    // matrix representation of convolution eq. 9
    for (int n = 0; n <= mo.getOrder(); n++) {
      for (int i = 0; i <= n; i++) {
	int  j = n-i;
	int n = mo0.getIndex(i,j);  // convolved moment index
	for (int k = 0; k <= i; k++) {
	  for (int l = 0; l <= j; l++) {
	    int m = mo0.getIndex(k,l); // unconvolved moment index
	    P(n,m) = binomial(i,k) * binomial(j,l) * p(i-k,j-l);
	  }
	}
      }
    }
    P.gemv(mo0, mem[k]);
  }

}
