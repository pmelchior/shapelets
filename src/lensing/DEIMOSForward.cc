#include "../../include/lensing/DEIMOSForward.h"
#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/MathHelper.h"

namespace shapelens {

  DEIMOSForward::DEIMOSForward(const MultiExposureObject& meo_, const MultiExposureObject& mepsf_, int N, int C, data_t flux, data_t width) :
    meo(meo_), mepsf(mepsf_), K(meo_.size()) {
    DEIMOS::N = N;
    DEIMOS::C = C;
    DEIMOS::D = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+C+1)*(N+C+2))/2);
    DEIMOS::S = NumMatrix<data_t>(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
    DEIMOS::mo.setOrder(DEIMOS::N);
    DEIMOS::matching_scale = width;

    // initialize 0th and 2nd moments 
    // with circular Gaussian with given flux & scale
    mo(0,0) = flux;
    // s = sqrt(trQ/F) = sqrt(2*Q_ii/F)
    mo(0,2) = mo(2,0) = 0;//width*width*mo(0,0)/2;

    // set up containers
    {
      Moments tmp(N);
      NumMatrix<data_t> P(((N+1)*(N+2))/2, ((N+1)*(N+2))/2);
      DEIMOS::PSFMultiScale psfs;
      for (int k = 0; k < K; k++) {
	mem.push_back(tmp);
	meP.push_back(P);
	meD.push_back(*this);
	// set noise image from each exposure
	meD[k].setNoiseImage(meo[k]);

	// create initial PSFMultiScale with given width
	DEIMOS psf(mepsf[k], N, C, width);
	psfs[width] = psf.mo;
	mePSFMultiScale.push_back(psfs);
      }
    }
    
    // Minimize chi^2
    data_t last_chi2 = 0; 
    int t = 1;
    while (t < 10) {
      history << t << "\t" << mo << std::endl;
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
	data_t chi2_k = diff * (S_1 * (NumVector<data_t>) diff);
	chi2 += chi2_k;
	NumMatrix<data_t> X = meP[k].transpose() * S_1;
	//mo += X * (NumVector<data_t>) d.mo;
	NumVector<data_t> mo_k = X * (NumVector<data_t>) d.mo;
	mo += mo_k;
	NumMatrix<data_t> S_k = X*meP[k];
	S += X*meP[k];
	//history << t << "." << k << "\t" << (X*meP[k]).invert()*mo_k << "\t" << chi2_k << "\t" << mo_k(0)/sqrt(S_k(0,0)) << std::endl;
      }
      S = S.invert();
      mo = S * (NumVector<data_t>)mo;
      SN[width] = mo(0,0)/sqrt(S(0,0)/K);
      history << "->\t" << mo << "\t" << shapelens::epsilon(mo) << "\t" << chi2 << "\t" << SN[width] << std::endl;
     
      // non-sensical ellipticity check
      data_t tiny = 1e-4*flux;
      mo(0,0) = std::max(tiny, mo(0,0));
      mo(0,2) = std::max(tiny, mo(0,2));
      mo(2,0) = std::max(tiny, mo(2,0));
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

      // convergence: chi^2 does not change by more than relative 1e-4
      // note: chi^2 might increase slightly before convergence
      // so go for asbolute value of change
      if (last_chi2 > 0 && fabs(last_chi2 - chi2) < 1e-4*chi2)
	break;
      last_chi2 = chi2;
      t++;
    }

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
      // FIXME: need individual WCS scale_factors here!
      // i.e. d.scale_factor needs to be set from outside
      d.eps = shapelens::epsilon(mo_c);
      d.matching_scale = getWeightFunctionScale(k);
      d.scale = d.getEpsScale();
      history << "\t" << k << "\t" << mo_c << "\t" << d.matching_scale << "\t" << d.eps << std::endl;
      DEIMOS::DEIMOSWeightFunction w(d.scale, meo[k].centroid, d.eps);
      Moments mo_w(meo[k], w, N+C);
      d.deweight(mo_w);
      d.computeCovariances(mo_w);
    }
  }

  void DEIMOSForward::convolveExposure(unsigned int k) {
    DEIMOS& d = meD[k];
    PSFMultiScale& psfs = mePSFMultiScale[k];

    // pick closest PSF scale to d.matching_scale;
    data_t scale = psfs.getScaleClosestTo(d.matching_scale);
    // if this if too different, create a new PSF moment measurement
    if (fabs(scale-d.matching_scale) > 0.05*d.matching_scale) {
      scale = d.matching_scale;
      DEIMOS psf(mepsf[k], N, C, scale);
      psfs.insert(scale, psf.mo);
    }
    Moments& p = psfs[scale];

    // create matrix representation of convolution eq. 9
    NumMatrix<data_t>& P = meP[k];
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
    // convolve moments
    P.gemv(mo, mem[k]);
  }
  
  
  data_t DEIMOSForward::getWeightFunctionScale(unsigned int k) const {
    // simply use sqrt(trQ/F) as size: equivalen width of Gaussian 
    const NumMatrix<data_t>& P = meP[k];
    const Moments& m = mem[k];
    data_t s = sqrt((m(0,2) + m(2,0))/m(0,0));
    // get the same for the current PSF moments
    const PSFMultiScale& psfs = mePSFMultiScale[k];
    // use the largest, most reliable PSF moments
    // FIXME: with noise in the PSF images, the larges may not be the best!
    const Moments& p = psfs.rbegin()->second;
    data_t s_psf = sqrt((p(0,2) + p(2,0))/p(0,0));

    // Noise correction, expect noise variance to be set
    if (S(0,0) > 0) {
      NumMatrix<data_t> S_c = P*S*P.transpose();
      // error on sqrt(trQ), ignoring flux uncertainty
      data_t sigma_2 = 0;
      unsigned int i = m.getIndex(2,0);
      sigma_2 += S_c(i,i);
      i = m.getIndex(0,2);
      sigma_2 += S_c(i,i);
      data_t d_s = 0.822179*sqrt(sigma_2/m(0,0));
      std::complex<data_t> eps = shapelens::epsilon(m);
      s -= d_s;
    }
    return std::max(s, s_psf); // don't undercut PSF equivalent width
  }
}
