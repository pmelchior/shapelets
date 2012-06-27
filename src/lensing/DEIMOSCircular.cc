#include "../../include/lensing/DEIMOSCircular.h"
#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/IO.h"
#include "../../include/ShapeLensConfig.h"

namespace shapelens {

  bool DEIMOSCircular::FIX_CENTROID = false;

  DEIMOSCircular::DEIMOSCircular() : DEIMOS(), scale(0), R2(0), scale_factor(1) {}

  DEIMOSCircular::DEIMOSCircular(const Object& obj_, int N_, int C_, data_t matching_scale_) : 
    DEIMOS(N_), C(C_), R2(0), matching_scale(matching_scale_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2) {
     
    // dirty little trick to update obj.centroid for moment measurement
    Object& obj = const_cast<Object&>(obj_);
    Point<data_t> old_centroid = obj.centroid; // save for later
    scale_factor = obj.grid.getScaleFactor();

    // measure moments within optimized weighting function
    Moments mo_w(N+C);
    match(obj, mo_w);

    // compute noise statistics
    setNoiseImage(obj);
    computeCovariances(mo_w);
    SN[matching_scale] = computeSN(mo_w);

    // restore obj's original centroid
    obj.centroid = old_centroid;
 }

  DEIMOSCircular::DEIMOSCircular (const Object& obj_, const DEIMOSCircular::PSFMultiScale& psf, int N_, int C_, data_t matching_scale_) :
    DEIMOS(N_), C(C_), matching_scale(matching_scale_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2) {
    
    // dirty little trick to update obj.centroid for moment measurement
    Object& obj = const_cast<Object&>(obj_);
    Point<data_t> old_centroid = obj.centroid; // save for later
    scale_factor = obj.grid.getScaleFactor();

    // measure moments within optimized weighting function
    Moments mo_w(N+C);
    setNoiseImage(obj);
    match(obj, mo_w);

    // initial matching failed, try with smaller scale
    while (flags.any() && matching_scale > psf.getMinimumScale()) {
      matching_scale = psf.getScaleSmallerThan(matching_scale);
      history << "# Matching failed (" << flags << "), restarting with s = " << matching_scale << std::endl; 
      match(obj, mo_w);
    }

    // matching successfull, try with lower scale and see whether S/N improves
    if (flags.none()) {
      while (matching_scale > psf.getMinimumScale()) {
	data_t SN_current = SN[matching_scale];
	data_t matching_scale_current = matching_scale;
	Point<data_t> centroid_current = centroid;
	Moments mo_current = mo;
	matching_scale = psf.getScaleSmallerThan(matching_scale);
	history << "# Trying smaller scale s = " << matching_scale << " to find optimal S/N" << std::endl;
	match(obj, mo_w);
	// if S/N goes down, reset to best S/N case
	if (SN[matching_scale] < SN_current) {
	  matching_scale = matching_scale_current;
	  history << "# Reverting to scale s = " << matching_scale << std::endl;
	  mo = mo_current;
	  scale = matching_scale * scale_factor;
	  centroid = centroid_current;
	  break;
	}
      }
      history << "# Deweighted moments:\t" << mo << std::endl;

      // deconvolve and check moments
      Moments tmp = mo;
      R2 = shapelens::R2(mo, psf.getAtScale(matching_scale));
      deconvolve(psf.getAtScale(matching_scale));
      flags[3] = flagMoments(mo);
      while (flags.any() && matching_scale > psf.getMinimumScale()) {
	matching_scale = psf.getScaleSmallerThan(matching_scale);
	history << "# Deconvolution failed, repeat matching with scale s = " << matching_scale << std::endl;
	match(obj, mo_w);
	history << "# Deweighted moments:\t" << mo << std::endl;
	tmp = mo;
	R2 = shapelens::R2(mo, psf.getAtScale(matching_scale));
	deconvolve(psf.getAtScale(matching_scale));
	flags[3] = flagMoments(mo);
      }
      history << "# Deconvolved moments:\t" << mo << std::endl;
      if (flags.any()) {
	history << "# Deconvolution failed, minimum PSF scale reached. GAME OVER." << std::endl;
	mo = tmp;
      }
      else
	computeCovariances(mo_w);
    }
    else
      history << "# Deweighting failed (" << flags << "), minimum PSF scale reached. GAME OVER." << std::endl;
    
    // restore obj's original centroid
    obj.centroid = old_centroid;

  }

  DEIMOSCircular::DEIMOSCircular(std::string filename) {
    fitsfile* fptr = IO::openFITSFile(filename);
    NumMatrix<data_t> M;
    IO::readFITSImage(fptr,M);
    N = int(M.getRows()) - 1;
    mo.setOrder(N);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	mo(n1,n2) = M(n2,n1);

    // necessary parameters
    IO::readFITSKeyword(fptr,"SCALE_M",matching_scale);
    IO::readFITSKeyword(fptr,"C",C);

    try {
      IO::readFITSKeyword(fptr,"SCALEFAC",scale_factor);
    } catch (std::invalid_argument) {
      scale_factor = 1;
    }
    try {
      IO::readFITSKeyword(fptr,"SCALE",scale);
    } catch (std::invalid_argument) {
      scale = matching_scale * scale_factor;
    }
    try {
      IO::readFITSKeyword(fptr,"SN",SN[matching_scale]);
    } catch (std::invalid_argument) {}
    try {
      IO::readFITSKeyword(fptr,"R2",R2);
    } catch (std::invalid_argument) {
      R2 = 0;
    }
    try {
      unsigned long flagnum;
      IO::readFITSKeyword(fptr,"FLAGS",flagnum);
      flags = std::bitset<4>(flagnum);
    } catch (std::invalid_argument) {}

    // read covariance
    try {
      IO::moveToFITSExtension(fptr, 2);
      IO::readFITSImage(fptr,S);
    } catch (std::runtime_error) {}
	
    IO::closeFITSFile(fptr);
  }

  void DEIMOSCircular::save(std::string filename) const {
    fitsfile* fptr = IO::createFITSFile(filename);
    NumMatrix<data_t> M(N+1,N+1);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	M(n2,n1) = mo(n1,n2); // transpose to have correctly oriented image
    IO::writeFITSImage(fptr,M);
    IO::updateFITSKeyword(fptr,"C",C,"deweighting correction order");
    IO::updateFITSKeyword(fptr,"SCALE_M",matching_scale,"matching_scale [pixel]");
    IO::updateFITSKeyword(fptr,"SCALE",scale,"actual weighting function width");
    IO::updateFITSKeyword(fptr,"SCALEFAC",scale_factor,"avg. WCS units/pixel");
    std::map<data_t, data_t>::const_iterator iter = SN.find(matching_scale);
    if (iter != SN.end())
      IO::updateFITSKeyword(fptr,"SN",iter->second,"S/N");
    IO::updateFITSKeyword(fptr,"R2", R2, "resolution of the galaxy");
    IO::updateFITSKeyword(fptr,"FLAGS",flags.to_ulong(),"matching and deconvolution flags");
    IO::writeFITSImage(fptr,S,"VARIANCE");
    IO::closeFITSFile(fptr);
  }

  // determine optimal weighting parameters, centroid and size
  void DEIMOSCircular::match(Object& obj, Moments& mo_w) {
    // obj.centroid is in pixel coordinates by convention
    // while shape measurement (including scale, centroid, ellipticity)
    // needs to be done in WCS coordinates.
    centroid = obj.centroid;
    if (ShapeLensConfig::USE_WCS)
      obj.grid.getWCS().transform(centroid);
    scale = matching_scale * scale_factor;
    computeDeweightingMatrix(mo_w);

    flags.reset();
    int iter = 0, maxiter = 8, iter_initial = 0;
    data_t SN_initial = 0;
    Point<data_t> centroid_shift;
    history << "# Matching weight function: s = " << matching_scale << std::endl;
    history << "# iter\tscale\tcentroid\t\tS/N" << std::endl;
    history << "# " + std::string(70, '-') << std::endl;
    while (iter < maxiter) {
      // measure moments under weight function
      GaussianWeightFunction w(scale, centroid);
      mo_w = Moments (obj, w, N + C);
      
      deweight(mo_w);
      flags[1] = flagMoments(mo);

      if (FIX_CENTROID)
	break;

      SN[matching_scale] = computeSN(mo_w);
      history << "# " << iter+1 << "\t" << scale/scale_factor << "\t" << obj.centroid;
      if (noise.size() > 0) { // only then SN_ is meaningful
	history << "\t";
	history << SN[matching_scale];
      }
      history << std::endl;

      // centroid correction
      centroid_shift(0) = mo(1,0)/mo(0,0);
      centroid_shift(1) = mo(0,1)/mo(0,0);
      data_t shift = sqrt(centroid_shift(0)*centroid_shift(0) + centroid_shift(1)*centroid_shift(1)) / scale_factor;
      if (shift > 5) {
	flags[0] = 1;
	break;
      }

      // shift centroid
      centroid += centroid_shift;
      obj.centroid = centroid;
      // obj.centroid is in pixel units (by convention)
      if (ShapeLensConfig::USE_WCS)
	obj.grid.getWCS().inverse_transform(obj.centroid); 

      // centroiding has converged: stop it
      if (shift < 1e-2*scale_factor)
	break;
      iter++;
    }
  }

  void DEIMOSCircular::computeDeweightingMatrix(const Moments& mo_w) {
    data_t s2 = scale*scale;
    data_t s4 = s2*s2;

    unsigned int i,j;
    for (int n=0; n <= N; n++) {
      for (int m=0; m <= n; m++) {
	i = mo.getIndex(m,n-m);
	j = mo_w.getIndex(m,n-m);
	D(i,j) = 1;
	if (C >= 2) {
	  j = mo_w.getIndex(m+2,n-m);
	  D(i,j) = 1./(2*s2);
	  // since the moments are of same order n
	  // and ordered wrt to (first/last) index
	  // moment (m+1,n-m+1) is just directly following in mo
	  j+=2;
	  D(i,j) = 1./(2*s2);
	}
	if (C >= 4) {
	  j = mo_w.getIndex(m+4,n-m);
	  D(i,j) = 1./(8*s4);
	  j+=2;
	  D(i,j) = 1./(4*s4);
	  j+=2;
	  D(i,j) = 1./(8*s4);
	}
	if (C >= 6) {
	  data_t s6 = s2*s2*s2;
	  j = mo_w.getIndex(m+6,n-m);
	  D(i,j) = 1./(48*s6);
	  j+=2;
	  D(i,j) = 1./(16*s6);
	  j+=2;
	  D(i,j) = 1./(16*s6);
	  j+=2;
	  D(i,j) = 1./(48*s6);
	}
	if (C >= 8) {
	  data_t s8 = s4*s4;
	  j = mo_w.getIndex(m+8,n-m);
	  D(i,j) = 1./(384*s8);
	  j+=2;
	  D(i,j) = 1./(96*s8);
	  j+=2;
	  D(i,j) = 1./(64*s8);
	  j+=2;
	  D(i,j) = 1./(96*s8);
	  j+=2;
	  D(i,j) = 1./(384*s8);
	}
	if (C >= 10) {
	  data_t s10 = s4*s4*s2;
	  j = mo_w.getIndex(m+10,n-m);
	  D(i,j) = 1./(3840*s10);
	  j+=2;
	  D(i,j) = 1./(768*s10);
	  j+=2;
	  D(i,j) = 1./(384*s10);
	  j+=2;
	  D(i,j) = 1./(384*s10);
	  j+=2;
	  D(i,j) = 1./(768*s10);
	  j+=2;
	  D(i,j) = 1./(3840*s10);
	}
      }
    }
  }
  void DEIMOSCircular::deweight(const Moments& mo_w) {
    // mo = D*mo_w
    D.gemv(mo_w, mo);
  }

  data_t DEIMOSCircular::computeSN(const Moments& mo_w) {
    // measure S/N (of flux)
    data_t SN_ = 0;
    // only do it if we have a noise image
    if (noise.size() > 0) {
      Moments mo_w2;
      GaussianWeightFunction w2(scale/M_SQRT2, centroid);
      mo_w2 = Moments(noise, w2, 0);
      SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
    }
    return SN_;
  }

  void DEIMOSCircular::computeCovariances(const Moments& mo_w) {
    // only do something if we have a noise image
    if (noise.size() > 0) {
      // covariances of a weighted moment mo(i,j) is given by
      // <mo(i,j) mo_(k,l)> = sigma_n^2 \int d^2x W^2(x) x_1^{i+k} x_2^{j+l},
      // i.e. moments up to (i*2,j*2) measured with square of weighting function
      // (square of Gaussian: sigma -> sigma/sqrt(2));
      GaussianWeightFunction w2(scale/M_SQRT2, centroid);
      Moments mo_noise(noise, w2, 2*(N+C));
      NumMatrix<data_t> S_(mo_w.size(), mo_w.size());
      data_t det = 1;
      for (int n = 0; n < mo_w.size(); n++) {
	std::pair<int, int> p = mo_w.getPowers(n);
	for (int m = 0; m < mo_w.size(); m++) {
	  std::pair<int, int> p_ = mo_w.getPowers(m);
	  if (p.first + p.second + p_.first + p_.second <= 2*(N+C))
	    S_(n,m) = mo_noise(p.first + p_.first, p.second + p_.second);
	}
	det *= S_(n,n);
      }
      S = (D*S_*D.transpose());
    }
  }
} // end namespace

