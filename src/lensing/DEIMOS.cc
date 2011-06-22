#include "../../include/lensing/DEIMOS.h"
#include "../../include/utils/IO.h"
#include "../../include/ShapeLensConfig.h"

namespace shapelens {

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
    std::map<data_t, Moments>::const_iterator iter = std::map<data_t, Moments>::end();
    iter--;
    return iter->first;
  }

  DEIMOS::DEIMOS() : id(0), scale(0), eps(0,0), G(0,0), N(0),  C(0), flexed(false) {}

  DEIMOS::DEIMOS (const Object& obj, int N_, int C_, data_t scale_, bool flexed_) :
    id(obj.id), scale(scale_), eps(0,0), N(N_), C(C_), flexed(flexed_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2),
    S(((N_+1)*(N_+2))/2, ((N_+1)*(N_+2))/2) {

    // compute scale in WCS units:
    // average scale within objects bounding box
    scale_factor = 1;
    if (ShapeLensConfig::USE_WCS) {
      scale_factor = obj.grid.getSupport().getArea();
      scale_factor /= obj.grid.getBoundingBox().getArea();
      scale_factor = sqrt(scale_factor);
      scale *= scale_factor;
    }

    // measure moments within optimized weighting function
    mo.setOrder(N);
    setNoiseImage(obj);
    data_t matching_scale = scale;
    match(obj, matching_scale);
    computeCovariances();
 }

  DEIMOS::DEIMOS (const Object& obj, const DEIMOS::PSFMultiScale& psf, int N_, int C_, data_t scale_, bool flexed_) :
    id(obj.id), scale(scale_), eps(0,0), N(N_), C(C_), flexed(flexed_),
    D(((N_+1)*(N_+2))/2, ((N_+C_+1)*(N_+C_+2))/2),
    S(((N_+1)*(N_+2))/2, ((N_+1)*(N_+2))/2) {

    // compute scale in WCS units:
    // average scale within objects bounding box
    scale_factor = 1;
    if (ShapeLensConfig::USE_WCS) {
      scale_factor = obj.grid.getSupport().getArea();
      scale_factor /= obj.grid.getBoundingBox().getArea();
      scale_factor = sqrt(scale_factor);
      scale *= scale_factor;
    }

    // measure moments within optimized weighting function
    mo.setOrder(N);
    data_t matching_scale = scale;
    setNoiseImage(obj);
    match(obj, matching_scale);

    // initial matching failed, try with smaller scale
    while (flags.any() && matching_scale > psf.getMinimumScale()) {
      matching_scale = psf.getScaleSmallerThan(matching_scale);
      history << "# Matching failed, restarting with scale s = " << matching_scale << std::endl; 
      match(obj, matching_scale);
    }

    // matching successfull, try with lower scale and see whether S/N improves
    if (flags.none()) {
      while (matching_scale/scale_factor > psf.getMinimumScale()) {
	data_t SN_current = SN[matching_scale];
	data_t matching_scale_current = matching_scale;
	complex<data_t> eps_current = eps, G_current = G;
	Point<data_t> centroid_current = centroid;
	Moments mo_current = mo;
	matching_scale = psf.getScaleSmallerThan(matching_scale/scale_factor)*scale_factor;
	history << "# Trying smaller scale s = " << matching_scale/scale_factor << " to find optimal S/N" << std::endl;
	match(obj, matching_scale);
	// if S/N goes down, reset to best S/N case
	if (SN[matching_scale] < SN_current) {
	  matching_scale = matching_scale_current;
	  history << "# Reverting to scale s = " << matching_scale/scale_factor << std::endl;
	  mo = mo_current;
	  eps = eps_current;
	  G = G_current;
	  scale = matching_scale/sqrt(1 + abs(eps)*abs(eps) - 2*abs(eps));
	  break;
	}
      }
      history << "# Deweighted moments:\t" << mo << std::endl;
      deconvolve(psf.getAtScale(matching_scale/scale_factor));
      while (flags.any() && 
	     matching_scale/scale_factor > psf.getMinimumScale()) {
	matching_scale = psf.getScaleSmallerThan(matching_scale/scale_factor)*scale_factor;
	history << "# Deconvolution failed, repeat matching with scale s = " << matching_scale/scale_factor << std::endl;
	// FIXME: if we store centroid, mo, eps, scale,
	// we may avoid matching in cases where smaller scales
	// have already been tried to improve S/N
	match(obj, matching_scale);
	history << "# Deweighted moments:\t" << mo << std::endl;
	deconvolve(psf.getAtScale(matching_scale/scale_factor));
      }
      history << "# Deconvolved moments:\t" << mo << std::endl;
      if (flags.any())
	history << "# Deconvolution failed, minimum PSF scale reached. GAME OVER." << std::endl;
      else
	computeCovariances();
    }
    else
      history << "# Deweighting failed, minimum PSF scale reached. GAME OVER." << std::endl;
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
    IO::readFITSKeyword(fptr,"WIDTH",scale);
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
      IO::readFITSKeyword(fptr,"SCALE",scale_factor);
    } catch (std::invalid_argument) {
      scale_factor = 1;
    }
    try {
      IO::readFITSKeyword(fptr,"SN",SN[getMatchingScale()]);
    } catch (std::invalid_argument) {}
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
    IO::updateFITSKeyword(fptr,"WIDTH",scale,"weighting function width");
    IO::updateFITSKeyword(fptr,"C",C,"deweighting correction order");
    IO::updateFITSKeyword(fptr,"SCALE",scale_factor,"avg. degrees/pixel");
    IO::updateFITSKeyword(fptr,"EPS",eps,"weighting function ellipticity");
    IO::updateFITSKeyword(fptr,"G",G,"weighting function G-flexion");
    std::map<data_t, data_t>::const_iterator iter = SN.find(getMatchingScale());
    if (iter != SN.end())
      IO::updateFITSKeyword(fptr,"SN",*iter,"S/N");
    IO::updateFITSKeyword(fptr,"FLAGS",flags.to_ulong(),"matching and deconvolution flags");
    IO::writeFITSImage(fptr,S,"VARIANCE");
    IO::closeFITSFile(fptr);
  }

  data_t DEIMOS::getEpsScale(data_t matching_scale, const std::complex<data_t>& eps_) const {
    // this scale preserves the size of semi-minor axis of ellipsoid
    // as the smallest scale (= matching_scale)
    return matching_scale/sqrt(1 + abs(eps_)*abs(eps_) - 2*abs(eps_));
  }

  data_t DEIMOS::getMatchingScale() const {
    return scale / scale_factor * sqrt(1 + abs(eps)*abs(eps) - 2*abs(eps));
  }
  
  // determine optimal weighting parameters, centroid, and ellipticity
  void DEIMOS::match(const Object& obj, data_t matching_scale) {
    // (re-)set centroid, ellipticity and scale (in WCS units if required)
    eps = complex<data_t>(0,0);
    centroid = obj.centroid;
    if (ShapeLensConfig::USE_WCS)
      obj.grid.getWCS().transform(centroid);
    flags.reset();

    // use a smaller scale for centroiding
    scale = std::max(matching_scale*0.66, 1.0*scale_factor);

    int iter = 0, maxiter = 18, maxiter_centroid = 4, iter_initial = 0, run = 0, maxrun = 1;
    bool centroiding = true;
    data_t SN_initial = 0;
    if (flexed) 
      maxrun++;
    
    Point<data_t> centroid_shift, centroid_pix;
    history << "# Matching weight function:" << std::endl;
    history << "# iter\tscale\tcentroid\t\tepsilon\t\t\tS/N" << std::endl;
    history << "# " + std::string(70, '-') << std::endl;
    Moments mo_w(N+C);
    while (true) {
      // measure moments under weight function
      if (flexed) {
	 DEIMOSWeightFunction w(scale, centroid, eps, G);
	mo_w = Moments (obj, w, N + C);
      }	else {
	DEIMOSWeightFunction w(scale, centroid, eps);
	mo_w = Moments (obj, w, N + C);
      }

      // measure S/N (of flux)
      data_t SN_;
      if (flexed) {
	DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps, G);
	Moments mo_w2(noise, w2, 0);
	SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
      } else {
	DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps);
	Moments mo_w2(noise, w2, 0);
	SN_ = mo_w(0,0) / sqrt(mo_w2(0,0));
      }
 
      // output current weighting parameters (in pixel units)
      centroid_pix = centroid;
      obj.grid.getWCS().inverse_transform(centroid_pix); 
      history << "# " << iter+1 << "\t" << scale/scale_factor << "\t" << centroid_pix << "\t" << eps << "\t";
      if (centroiding)
	history << "\t\t";
      history << SN_ << "\t" << SN_initial << std::endl;

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
	if (centroiding && iter < maxiter_centroid - 1)  {
	  centroid_shift(0) = mo(1,0)/mo(0,0);
	  centroid_shift(1) = mo(0,1)/mo(0,0);
	  data_t shift = sqrt(centroid_shift(0)*centroid_shift(0) + centroid_shift(1)*centroid_shift(1)) / scale_factor;
	  if (shift > 10) {
	    flags[0] = 1;
	    break;
	  }
	  if (iter == maxiter_centroid - 1 && trouble) {
	    flags[1] = trouble;
	    break;
	  }
	  centroid -= centroid_shift;
	  // centroiding has converged: stop it
	  if (shift < 1e-2*scale_factor)
	    centroiding = false;
	}
	// match ellipticity
	else if (iter < maxiter - 1) {
	  complex<data_t> eps_ = epsilon();
	  // abort for non-sensical ellipticities
	  if (trouble) {
	    flags[1] = trouble;
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
	    if (SN_ < (0.95*SN_initial)) { // strong non-convergence
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
	  scale = getEpsScale(matching_scale, eps);
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
    data_t s6 = s2*s2*s2;
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
      }
    }

    // mo = D*mo_w
    D.gemv(mo_w, mo);
  }

  void DEIMOS::setNoiseImage(const Object& obj) {
    noise.resize(obj.size());
    noise.grid = obj.grid;
    noise.centroid = obj.centroid;
    if (obj.weight.size() == 0)
      for (unsigned int i=0; i < noise.size(); i++)
	noise(i) = obj.noise_rms*obj.noise_rms;
    else
      for (unsigned int i=0; i < noise.size(); i++)
	noise(i) = obj.weight(i);
  }

  void DEIMOS::computeCovariances() {
    // variance of weighted moment (i,j) is propto
    // moment (i*2,j*2) measured with square of weighting function
    // (square of Gaussian: sigma -> sigma/sqrt(2));
    DEIMOSWeightFunction w2(scale/M_SQRT2, centroid, eps);
    Moments mo_noise(noise,w2,2*(N+C));

    // copy terms from (2*i, 2*j) to (i,j)
    for (int n=1; n <= N+C; n++)
      for (int m=0; m <= n; m++)
      mo_noise(m,n-m) = mo_noise(2*m,2*(n-m));
    mo_noise.setOrder(N+C);

    S.clear();
    for (int i=0; i < mo.size(); i++)
      for (int j=0; j < mo.size(); j++)
	for (unsigned int k=0; k < mo_noise.size(); k++)
	  S(i,j) += D(i,k)*D(j,k)*mo_noise(k);
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


  // CAUTION: this overflows for n > 10!!!
  unsigned long factorial(int n) {
    unsigned long f = 1;
    for (int m=2; m <= n; m++)
      f *= m;
    return f;
  }
  unsigned long binomial(int n, int m) {
    return factorial(n)/(factorial(m)*factorial(n-m));
  }

  void DEIMOS::deconvolve(const Moments& p) {
    Moments& g = mo;
    int Nmin = std::min(p.getOrder(),g.getOrder());
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
    if (mo.getOrder() >= 2) {
      complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
      e/= (complex<data_t>(mo(2,0) + mo(0,2)) + 2.*sqrt(complex<data_t>(mo(0,2)*mo(2,0) - mo(1,1)*mo(1,1))));
      return e;
    } else
      return complex<data_t>(0,0);
  }
  
  complex<data_t> DEIMOS::chi() const {
    if (mo.getOrder() >= 2) {
      complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
      e/= mo(2,0) + mo(0,2);
      return e;
    } else
      return complex<data_t>(0,0);
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

  
  // #### DEIMOSList ####
  DEIMOSList::DEIMOSList() : std::vector<boost::shared_ptr<DEIMOS> >() {}

#ifdef HAS_SQLiteDB
  DEIMOSList::DEIMOSList(SQLiteDB& sql, std::string table, std::string where) :
    std::vector<boost::shared_ptr<DEIMOS> > () {
    
  }
  void DEIMOSList::save(SQLiteDB& sql, std::string table) const {
    // drop table: faster than deleting its entries
    std::string query = "DROP TABLE IF EXISTS " + table + ";";
    sql.query(query);

    // create it
    query = "CREATE TABLE " + table + "(";
    query += "`id` int NOT NULL PRIMARY KEY,";
    query += "`width` float NOT NULL,";
    query += "`eps1` float NOT NULL,";
    query += "`eps2` float NOT NULL,";
    query += "`G1` float NOT NULL,";
    query += "`G2` float NOT NULL,";
    query += "`c` int NOT NULL,";
    query += "`mo` blob NOT NULL);";
    sql.query(query); 
    
    // create prepared statement
    sqlite3_stmt *stmt;
    query = "INSERT INTO `" + table + "` VALUES (?,?,?,?,?,?,?,?);";
    sql.checkRC(sqlite3_prepare_v2(sql.db, query.c_str(), query.size(), &stmt, NULL));
    for(DEIMOSList::const_iterator iter = DEIMOSList::begin(); iter != DEIMOSList::end(); iter++) {
      const DEIMOS& d = *(*iter);
      sql.checkRC(sqlite3_bind_int(stmt,1,d.id));
      sql.checkRC(sqlite3_bind_double(stmt,2,d.scale));
      sql.checkRC(sqlite3_bind_double(stmt,3,real(d.eps)));
      sql.checkRC(sqlite3_bind_double(stmt,4,imag(d.eps)));
      sql.checkRC(sqlite3_bind_double(stmt,5,real(d.G)));
      sql.checkRC(sqlite3_bind_double(stmt,6,imag(d.G)));
      sql.checkRC(sqlite3_bind_int(stmt,7,d.C));
      sql.checkRC(sqlite3_bind_blob(stmt,8,d.mo.c_array(),d.mo.size()*sizeof(data_t),SQLITE_STATIC));
      if(sqlite3_step(stmt)!=SQLITE_DONE)
	throw std::runtime_error("DEIMOSList::save() insertion failed: " + std::string(sqlite3_errmsg(sql.db)));
      sql.checkRC(sqlite3_reset(stmt));
    }
    sql.checkRC(sqlite3_finalize(stmt));
  }
#endif  

} // end namespace
