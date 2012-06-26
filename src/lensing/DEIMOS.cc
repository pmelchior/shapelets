#include "../../include/lensing/DEIMOS.h"
#include "../../include/lensing/LensHelper.h"
#include "../../include/utils/IO.h"
#include "../../include/utils/MathHelper.h"
#include "../../include/ShapeLensConfig.h"

namespace shapelens {

  bool DEIMOS::FIX_CENTROID = false;
  bool DEIMOS::FIX_ELLIPTICITY = false;

  DEIMOS::DEIMOSWeightFunction::DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid_, const complex<data_t>& eps) :
    GaussianWeightFunction(scale,Point<data_t>(0,0)),
    centroid(centroid_),
    T(0,eps) { 
  }

  DEIMOS::DEIMOSWeightFunction::DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid_, const complex<data_t>& eps, const complex<data_t>& G) :
    GaussianWeightFunction(scale,Point<data_t>(0,0)),
    centroid(centroid_),
    T(0,eps,complex<data_t>(0,0), G) { 
  }

  data_t DEIMOS::DEIMOSWeightFunction::operator()(const Point<data_t>& P_) const {
    Point<data_t> P = P_ - centroid;
    T.inverse_transform(P);
    return GaussianWeightFunction::operator()(P);
  }
  
  void DEIMOS::PSFMultiScale::insert(data_t scale, const Moments& mo) {
    std::map<data_t, Moments>::insert(std::pair<data_t, Moments>(scale, mo));
  }
  const Moments& DEIMOS::PSFMultiScale::getAtScale(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::find(scale);
    if (iter == std::map<data_t, Moments>::end())
      throw std::invalid_argument("DEIMOS::PSFMultiScale: scale not available!");
    return iter->second;
  }

  data_t DEIMOS::PSFMultiScale::getScaleSmallerThan(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::find(scale);
    if (iter == std::map<data_t, Moments>::end())
      throw std::invalid_argument("DEIMOS::PSFMultiScale: scale not available!");
    if (iter == std::map<data_t, Moments>::begin())
      throw std::runtime_error("DEIMOS::PSFMultiScale: no smaller scale available!");
    iter--;
    return iter->first;

  }

  data_t DEIMOS::PSFMultiScale::getMinimumScale() const {
    return std::map<data_t, Moments>::begin()->first;
  }
  
  data_t DEIMOS::PSFMultiScale::getMaximumScale() const {
    return std::map<data_t, Moments>::rbegin()->first;
  }

  data_t DEIMOS::PSFMultiScale::getScaleClosestTo(data_t scale) const {
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::lower_bound(scale);
    // iter points to element >= scale;
    data_t diff = (iter->first)-scale;
    if (diff < 0) { // bigger than largest element
      return std::map<data_t, Moments>::rbegin()->first;
    }
    else {
      // check whethere there is a smaller scale
      if (iter != std::map<data_t, Moments>::begin()) {
	iter--;
	data_t diff_small = scale - (iter->first);
	if (diff < diff_small)
	  iter++;
      }
      return iter->first;
    }
  }

  DEIMOS::DEIMOS() : id(0), scale(0), eps(0,0), G(0,0), N(0),  C(0), R2(0), flexed(false), scale_factor(1) {}

  DEIMOS::DEIMOS (const Object& obj_, int N_, int C_, data_t matching_scale_, bool flexed_) :
    id(obj_.id), matching_scale(matching_scale_), eps(0,0), N(N_), C(C_), R2(0), flexed(flexed_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2),
    S(((N_+1)*(N_+2))/2, ((N_+1)*(N_+2))/2) {

    // dirty little trick to update obj.centroid for moment measurement
    Object& obj = const_cast<Object&>(obj_);
    Point<data_t> old_centroid = obj.centroid; // save for later
    setScaleFactor(obj);

    // measure moments within optimized weighting function
    mo.setOrder(N);
    Moments mo_w(N+C);
    match(obj, mo_w);

    setNoiseImage(obj);
    computeCovariances(mo_w);
    SN[matching_scale] = computeSN(mo_w);

    // restore obj's original centroid
    obj.centroid = old_centroid;
 }

  DEIMOS::DEIMOS (const Object& obj_, const DEIMOS::PSFMultiScale& psf, int N_, int C_, data_t matching_scale_, bool flexed_) :
    id(obj_.id), matching_scale(matching_scale_), eps(0,0), N(N_), C(C_), R2(0), flexed(flexed_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2),
    S(((N_+1)*(N_+2))/2, ((N_+1)*(N_+2))/2) {
    
    // dirty little trick to update obj.centroid for moment measurement
    Object& obj = const_cast<Object&>(obj_);
    Point<data_t> old_centroid = obj.centroid; // save for later
    setScaleFactor(obj);

    // measure moments within optimized weighting function
    mo.setOrder(N);
    Moments mo_w(N+C);
    match(obj, mo_w);
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
	complex<data_t> eps_current = eps, G_current = G;
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
	  eps = eps_current;
	  G = G_current;
	  scale = getEpsScale();
	  centroid = centroid_current;
	  break;
	}
      }
      history << "# Deweighted moments:\t" << mo << std::endl;
      Moments tmp = mo;
      deconvolve(psf.getAtScale(matching_scale));
      while (flags.any() && matching_scale > psf.getMinimumScale()) {
	matching_scale = psf.getScaleSmallerThan(matching_scale);
	history << "# Deconvolution failed, repeat matching with scale s = " << matching_scale << std::endl;
	// FIXME: if we store centroid, mo, eps, scale,
	// we may avoid matching in cases where smaller scales
	// have already been tried to improve S/N
	match(obj, mo_w);
	history << "# Deweighted moments:\t" << mo << std::endl;
	tmp = mo;
	deconvolve(psf.getAtScale(matching_scale));
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

  DEIMOS::DEIMOS(std::string filename) {
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
    // optional parameters
    try {
      IO::readFITSKeyword(fptr,"ID",id);
    } catch (std::invalid_argument) {
      id = 0;
    }
    try {
      IO::readFITSKeyword(fptr,"EPS",eps);
    } catch (std::invalid_argument) {
      eps = complex<data_t>(0,0);
    }
    try {
      IO::readFITSKeyword(fptr,"G",G);
    } catch (std::invalid_argument) {
      G = complex<data_t>(0,0);
    }
    try {
      IO::readFITSKeyword(fptr,"SCALEFAC",scale_factor);
    } catch (std::invalid_argument) {
      scale_factor = 1;
    }
    try {
      IO::readFITSKeyword(fptr,"SCALE",scale);
    } catch (std::invalid_argument) {
      scale = getEpsScale();
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

  void DEIMOS::save(std::string filename) const {
    fitsfile* fptr = IO::createFITSFile(filename);
    NumMatrix<data_t> M(N+1,N+1);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	M(n2,n1) = mo(n1,n2); // transpose to have correctly oriented image
    IO::writeFITSImage(fptr,M);
    IO::updateFITSKeyword(fptr,"ID",id,"object ID");
    IO::updateFITSKeyword(fptr,"C",C,"deweighting correction order");
    IO::updateFITSKeyword(fptr,"SCALE_M",matching_scale,"matching_scale [pixel]");
    IO::updateFITSKeyword(fptr,"SCALE",scale,"actual weighting function width");
    IO::updateFITSKeyword(fptr,"SCALEFAC",scale_factor,"avg. WCS units/pixel");
    IO::updateFITSKeyword(fptr,"EPS",eps,"weighting function ellipticity");
    IO::updateFITSKeyword(fptr,"G",G,"weighting function G-flexion");
    std::map<data_t, data_t>::const_iterator iter = SN.find(matching_scale);
    if (iter != SN.end())
      IO::updateFITSKeyword(fptr,"SN",iter->second,"S/N");
    IO::updateFITSKeyword(fptr,"R2", R2, "resolution of the galaxy");
    IO::updateFITSKeyword(fptr,"FLAGS",flags.to_ulong(),"matching and deconvolution flags");
    IO::writeFITSImage(fptr,S,"VARIANCE");
    IO::closeFITSFile(fptr);
  }

  // compute scale in WCS units:
  // average scale within objects bounding box
  void DEIMOS::setScaleFactor(const Object& obj) {
    scale_factor = 1;
    if (ShapeLensConfig::USE_WCS) {
      scale_factor = obj.grid.getSupport().getArea();
      scale_factor /= obj.grid.getBoundingBox().getArea();
      scale_factor = sqrt(scale_factor);
    }
  }

  data_t DEIMOS::getEpsScale() const {
    // this scale preserves the size of semi-minor axis of ellipsoid
    // as the smallest scale (= matching_scale)
    data_t abs_eps = abs(eps);
    return scale_factor*matching_scale/sqrt(1 + abs_eps*abs_eps - 2*abs_eps);
  }

  // determine optimal weighting parameters, centroid, and ellipticity
  void DEIMOS::match(Object& obj, Moments& mo_w) {
    // (re-)set centroid, ellipticity and scale
    centroid = obj.centroid;
    eps = complex<data_t>(0,0);
    flags.reset();

    // obj.centroid is in pixel coordinates by convention
    // while shape measurement (including scale, centroid, ellipticity)
    // needs to be done in WCS coordinates.
    if (ShapeLensConfig::USE_WCS)
      obj.grid.getWCS().transform(centroid);

    // use a smaller scale for centroiding
    // but maintain a sensible width (of 1 pixel) otherwise the
    // cenroid could become dominated by local noise fluctuations 
    scale = std::max(matching_scale*0.66, 1.0) * scale_factor;

    int iter = 0, maxiter = 18, maxiter_centroid = 4, iter_initial = 0, run = 0, maxrun = 1;
    bool centroiding = true;
    data_t SN_initial = 0;
    if (flexed) 
      maxrun++;
    
    Point<data_t> centroid_shift;
    history << "# Matching weight function: s = " << matching_scale << std::endl;
    history << "# iter\tscale\tcentroid\t\tepsilon\t\t\tS/N" << std::endl;
    history << "# " + std::string(70, '-') << std::endl;
    while (true) {
      // measure moments under weight function
      if (flexed) {
	 DEIMOSWeightFunction w(scale, centroid, eps, G);
	mo_w = Moments (obj, w, N + C);
      }	else {
	DEIMOSWeightFunction w(scale, centroid, eps);
	mo_w = Moments (obj, w, N + C);
      }

      data_t SN_ = computeSN(mo_w);
      history << "# " << iter+1 << "\t" << scale/scale_factor << "\t" << obj.centroid << "\t" << eps;
      if (noise.size() > 0) { // only then SN_ is meaningful
	history << "\t";
	if (centroiding)
	  history << "\t\t";
	history << SN_;
      }
      history << std::endl;

      // deweight and estimate new centroid and ellipticity
      deweight(mo_w);
      bool trouble = flagMoments(mo);

      if (flexed) {
	throw std::runtime_error("DEIMOS: flexion matching not implemented");
	/*
	complex<data_t> delta_ = delta();
	if (iter < maxiter/3 - 1)
	  centroid += centroid_shift;
	else if (iter < 2*maxiter/3 - 1) {
	  eps = eps_;
	  scale = getEpsScale(matching_scale, eps);
	} else if (iter < maxiter - 1) {
	  G = (4./3) * delta_;
	}
	*/
      }
      else {
	// centroiding
	if (centroiding && iter < maxiter_centroid - 1 && FIX_CENTROID == false)  {
	  centroid_shift(0) = mo(1,0)/mo(0,0);
	  centroid_shift(1) = mo(0,1)/mo(0,0);
	  data_t shift = sqrt(centroid_shift(0)*centroid_shift(0) + centroid_shift(1)*centroid_shift(1)) / scale_factor;
	  if (shift > 5) {
	    flags[0] = 1;
	    mo(0,0) = mo_w(0,0);
	    mo(0,2) = mo_w(0,2);
	    mo(1,1) = mo_w(1,1);
	    mo(2,0) = mo_w(2,0);
	    break;
	  }
	  if (iter == maxiter_centroid - 1 && trouble) {
	    flags[1] = trouble;
	    mo(0,0) = mo_w(0,0);
	    mo(0,2) = mo_w(0,2);
	    mo(1,1) = mo_w(1,1);
	    mo(2,0) = mo_w(2,0);
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
	    centroiding = false;
	}
	// match ellipticity
	else if (iter < maxiter - 1 && FIX_ELLIPTICITY == false) {
	  complex<data_t> eps_ = epsilon();

	  // abort for non-sensical ellipticities
	  if (trouble) {
	    flags[1] = trouble;
	    mo(0,0) = mo_w(0,0);
	    mo(0,2) = mo_w(0,2);
	    mo(1,1) = mo_w(1,1);
	    mo(2,0) = mo_w(2,0);
	    break;
	  }
	  // convergence test
	  if (iter_initial == 0) {
	    iter_initial = iter;
	    centroiding = false;
	  }
	  if (iter == iter_initial + 1)
	    SN_initial = SN_;
	  if (iter > iter_initial + 1) {
	    // strong non-convergence: don't cut too hard here, it saves time but objects
	    // with large apparent ellipticity see an decrease in S/N during matching
	    if (SN_ < (0.5*SN_initial)) { 
	      flags[2] = 1;
	      break;
	    }
	    if (fabs(SN_/SN[matching_scale] - 1) < 1e-5 && 
		abs(eps - eps_) < 1e-3) {
	      SN[matching_scale] = SN_;
	      break;
	    }
	  }
	  eps = eps_;
	  scale = getEpsScale();
	}
	SN[matching_scale] = SN_;
      }

      // repeat same series of iteration again
      // but only one other time
      if (iter == maxiter-1) {
	run++;
	if (run == maxrun)
	  break;
	else 
	  iter = 0;
      }
      iter++;
    }
  }

  void DEIMOS::deweight(const Moments& mo_w) {
    data_t e1 = real(eps);
    data_t e2 = imag(eps);
    data_t c1 = (1-e1)*(1-e1) + e2*e2;
    data_t c2 = (1+e1)*(1+e1) + e2*e2;
    data_t s2 = scale*scale;
    data_t s4 = s2*s2;
    data_t G1 = real(G);
    data_t G2 = imag(G);

    unsigned int i,j;
    for (int n=0; n <= N; n++) {
      for (int m=0; m <= n; m++) {
	i = mo.getIndex(m,n-m);
	j = mo_w.getIndex(m,n-m);
	D(i,j) = 1;
	if (C >= 2) {
	  j = mo_w.getIndex(m+2,n-m);
	  D(i,j) = c1/(2*s2);
	  // since the moments are of same order n
	  // and ordered wrt to (first/last) index
	  // moment (m+1,n-m+1) is just directly following in mo
	  j++;
	  D(i,j) = - 2*e2/s2;
	  j++;
	  D(i,j) = c2/(2*s2);
	}
	if (C >= 3 && flexed) {
	  j = mo_w.getIndex(m+3,n-m);
	  D(i,j) = (-G1 + e1*G1 + e2*G2)/(4*s4);
	  j++;
	  D(i,j) = (-e2*G1 - 3*G2 + e1*G2)/(4*s4);
	  j++;
	  D(i,j) = (3*G1 + e1*G1 + e2*G2)/(4*s4);
	  j++;
	  D(i,j) = (-e2*G1 + G2 + e1*G2)/(4*s4);
	}
	if (C >= 4) {
	  j = mo_w.getIndex(m+4,n-m);
	  D(i,j) = c1*c1/(8*s4);
	  if (flexed)
	    D(i,j) += (G1*G1 + G2*G2)/(8*s2);
	  j++;
	  D(i,j) = -c1*e2/s4;
	  j++;
	  D(i,j) = (c1*c2/4 + 2*e2*e2)/s4;
	  if (flexed)
	    D(i,j) += (G1*G1 + G2*G2)/(8*s2);
	  j++;
	  D(i,j) = -c2*e2/s4;
	  j++;
	  D(i,j) = c2*c2/(8*s4);
	  if (flexed)
	    D(i,j) += (G1*G1 + G2*G2)/(8*s2);
	}
	if (C >= 6) {
	  data_t s6 = s2*s2*s2;
	  j = mo_w.getIndex(m+6,n-m);
	  D(i,j) = c1*c1*c1/(48*s6);
	  j++;
	  D(i,j) = -c1*c1*e2/(4*s6);
	  j++;
	  D(i,j) = (c1*c1*c2/16 + c1*e2*e2)/s6;
	  j++;
	  D(i,j) = -(c1*c2*e2/2 + 4*e2*e2*e2/3)/s6;
	  j++;
	  D(i,j) = (c1*c2*c2/16 + c2*e2*e2)/s6;
	  j++;
	  D(i,j) = -c2*c2*e2/(4*s6);
	  j++;
	  D(i,j) = c2*c2*c2/(48*s6);
	}
	if (C >= 8) {
	  data_t s8 = s4*s4;
	  j = mo_w.getIndex(m+8,n-m);
	  D(i,j) = c1*c1*c1*c1/(384*s8);
	  j++;
	  D(i,j) = -c1*c1*c1*e2/(24*s8);
	  j++;
	  D(i,j) = (c1*c1*c1*c2/96 + c1*c1*e2*e2/4)/s8;
	  j++;
	  D(i,j) = -(c1*c1*c2*e2/8 + 2*c1*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = (c1*c1*c2*c2/64 + c1*c2*e2*e2/2 + 2*e2*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = -(c1*c2*c2*e2/8 + 2*c2*e2*e2*e2/3)/s8;
	  j++;
	  D(i,j) = (c1*c2*c2*c2/96 + c2*c2*e2*e2/4)/s8;
	  j++;
	  D(i,j) = -c2*c2*c2*e2/(24*s8);
	  j++;
	  D(i,j) = c2*c2*c2*c2/(384*s8);
	}
	if (C >= 10) {
	  data_t s10 = s4*s4*s2;
	  j = mo_w.getIndex(m+10,n-m);
	  D(i,j) = c1*c1*c1*c1*c1/(3840*s10);
	  j++;
	  D(i,j) = -c1*c1*c1*c1*e2/(192*s10);
	  j++;
	  D(i,j) = (c1*c1*c1*c1*c2/768 + c1*c1*c1*e2*e2/24)/s10;
	  j++;
	  D(i,j) = -(c1*c1*c1*c2*e2/48 + c1*c1*e2*e2*e2/6)/s10;
	  j++;
	  D(i,j) = (c1*c1*c1*c2*c2/384 + c1*c1*c2*e2*e2/8 + c1*e2*e2*e2*e2/3)/s10;
	  j++;
	  D(i,j) = -(c1*c1*c2*c2*e2/32 + c1*c2*e2*e2*e2/3 + 4*e2*e2*e2*e2*e2/15)/s10;
	  j++;
	  D(i,j) = (c1*c1*c2*c2*c2/384 + c1*c2*c2*e2*e2/8 + c2*e2*e2*e2*e2/3)/s10;
	  j++;
	  D(i,j) = -(c1*c2*c2*c2*e2/48 + c2*c2*e2*e2*e2/6)/s10;
	  j++;
	  D(i,j) = (c1*c2*c2*c2*c2/768 + c2*c2*c2*e2*e2/24)/s10;
	  j++;
	  D(i,j) = -c2*c2*c2*c2*e2/(192*s10);
	  j++;
	  D(i,j) = c2*c2*c2*c2*c2/(3840*s10);
	}
      }
    }

    // mo = D*mo_w
    D.gemv(mo_w, mo);
  }

  void DEIMOS::setNoiseImage(const Object& obj) {
    // only do something if noise properties are set
    if (obj.noise_rms > 0) {
      noise.resize(obj.size());
      noise.grid = obj.grid;
      noise.centroid = obj.centroid;
      noise.noise_rms = obj.noise_rms;
      for (unsigned int i=0; i < noise.size(); i++)
	noise(i) = obj.noise_rms*obj.noise_rms;
      if (obj.segmentation.size() == obj.size())
	noise.segmentation = obj.segmentation;
    }
  }

  data_t DEIMOS::computeSN(const Moments& mo_w) {
    // measure S/N (of flux)
    data_t SN_ = 1;
    // only do it if we have a noise image
    if (noise.size() > 0) {
      Moments mo_w2;
      if (flexed) {
	DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps, G);
	mo_w2 = Moments(noise, w2, 0);
      } else {
	DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps);
	mo_w2 = Moments(noise, w2, 0);
      }	
      SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
    }
    return SN_;
  }

  void DEIMOS::computeCovariances(const Moments& mo_w) {
    // only do something if we have a noise image
    if (noise.size() > 0) {

      // covariances of a weighted moment mo(i,j) is given by
      // <mo(i,j) mo_(k,l)> = sigma_n^2 \int d^2x W^2(x) x_1^{i+k} x_2^{j+l},
      // i.e. moments up to (i*2,j*2) measured with square of weighting function
      // (square of Gaussian: sigma -> sigma/sqrt(2));
      DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps);
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

  Moments DEIMOS::getMomentErrors() const { 
    Moments mo_noise(N);
    if (S(0,0) != 0) {
      NumMatrix<data_t> S_1 = S.invert();
      for (int i=0; i < mo_noise.size(); i++)
	mo_noise(i) = 1./S_1(i,i);
    }
    return mo_noise;
  }

  void DEIMOS::deconvolve(const Moments& p) {
    Moments& g = mo;
    int Nmin = std::min(p.getOrder(),g.getOrder());

    // compute resolution factor R_2 (Hirata04.1, eq. 8)
    R2 = 1 - ((p(2,0)+p(0,2))/p(0,0)) / ((g(2,0)+g(0,2))/g(0,0));

    //use explicit relations for up to 2nd moments
    g(0,0) /= p(0,0);
    if (Nmin >= 1) {
      g(0,1) -= g(0,0)*p(0,1);
      g(0,1) /= p(0,0);
      g(1,0) -= g(0,0)*p(1,0);
      g(1,0) /= p(0,0);
      if (Nmin >= 2) {
	g(0,2) -= g(0,0)*p(0,2) + 2*g(0,1)*p(0,1);
	g(0,2) /= p(0,0);
	g(1,1) -= g(0,0)*p(1,1) + g(0,1)*p(1,0) + g(1,0)*p(0,1);
	g(1,1) /= p(0,0);
	g(2,0) -= g(0,0)*p(2,0) + 2*g(1,0)*p(1,0);
	g(2,0) /= p(0,0);
	if (Nmin >= 3) {
	  //use general formula (9)
	    for (int n=3; n <= Nmin; n++) {
	      for (int i=0; i <= n; i++) {
		int j = n-i;
		for (int k=0; k <= i-1; k++)
		  for (int l=0; l <= j-1; l++)
		    g(i,j) -= binomial(i,k)*binomial(j,l)*g(k,l)*p(i-k,j-l);
		for (int k=0; k <= i-1; k++)
		  g(i,j) -= binomial(i,k)*g(k,j)*p(i-k,0);
		for (int l=0; l <= j-1; l++)
		  g(i,j) -= binomial(j,l)*g(i,l)*p(0,j-l);
		g(i,j) /= p(0,0);
	      }
	    }
	}
      }
    }
    // check deconvolved moments
    flags[3] = flagMoments(mo);
  }
  
  // check whether moments make sense
  bool DEIMOS::flagMoments(const Moments& M) const {
    if (M(0,0) < 0 || M(2,0) < 0 || M(0,2) < 0 || M(1,1)*M(1,1) > M(2,0)*M(0,2))
      return true;
    // FIXME: for higher order moments could also check M(4,0) ...
    else
      return false;
  }

  complex<data_t> DEIMOS::epsilon() const {
    return shapelens::epsilon(mo);
  }
  
  complex<data_t> DEIMOS::chi() const {
    return shapelens::chi(mo);
  }
 
  // HOLICS equations from OMU 2007
  complex<data_t> DEIMOS::zeta() const {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4); // eq. 25
      complex<data_t> zeta(mo(3,0) + mo(1,2),    // eq. 26
			   mo(2,1) + mo(0,3));
      zeta /= xi;
      return zeta;
    } else
      return complex<data_t>(0,0);
  }

  complex<data_t> DEIMOS::delta() const {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4);  // eq. 25
      complex<data_t> delta(mo(3,0) - 3*mo(1,2),  // eq. 27
			    3*mo(2,1) - mo(0,3));
      delta /= xi;
      return delta;
    } else
      return complex<data_t>(0,0);
  }

} // end namespace
