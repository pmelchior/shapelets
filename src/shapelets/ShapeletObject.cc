#include <shapelets/ShapeletObject.h>
#include <shapelets/SIFFile.h>
#include <stdio.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace shapelens;
using namespace std;
typedef complex<data_t> Complex;

ShapeletObject::ShapeletObject() : 
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov) {
  chisquare = id = 0;
  fits = false;
  updatePolar = true;
  flags.reset();
  basefilename = "";
}

ShapeletObject::ShapeletObject(const ShapeletObject& source) : 
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov) {
  ShapeletObject::operator=(source);
}

ShapeletObject& ShapeletObject::operator=(const ShapeletObject& source) {
  Composite2D::operator=(source);
  // copy every variable of this class
  // base class copy constructor already called
  polarCoeffs = source.polarCoeffs;
  c2p = source.c2p;
  trafo = source.trafo;
  chisquare = source.chisquare;
  id = source.id;
  fits = source.fits;
  updatePolar = source.updatePolar;
  flags = source.flags;
  basefilename = source.basefilename;
  history = source.history;
  prop = source.prop;
  return *this;
}

ShapeletObject::ShapeletObject(string sifFile, bool preserve_config) : 
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov) {
  // get infos from file
  load(sifFile,preserve_config);
  updatePolar = true;
  fits = false;
}

ShapeletObject::ShapeletObject(const CoefficientVector<data_t>& incoeffs, data_t beta, const Point2D<data_t>& xcentroid, const Grid& grid) :
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov) {
  chisquare = id = 0;
  fits = false;
  updatePolar = true;
  flags.reset();
  basefilename = "";
  coeffs = incoeffs;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  Composite2D::setGrid(grid);
}

ShapeletObject::ShapeletObject(const CoefficientVector<Complex>& incoeffs, data_t beta, const Point2D<data_t>& xcentroid, const Grid& grid):
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov)  {
  chisquare = id = 0 ;
  fits = false;
  updatePolar = false;
  flags.reset();
  basefilename = "";
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  Composite2D::setGrid(grid);
  polarCoeffs = incoeffs;
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
}

ShapeletObject::ShapeletObject(const Object& obj) : 
Composite2D(), coeffs(Composite2D::coeffs), cov(Composite2D::cov) {
  fits = true;
  updatePolar = true;
  id = obj.id;
  basefilename = obj.basefilename;
  // in legacy mode: uncomment the following line
  // prop["classifier"] = obj.classifier;

  // optimized decomposition
  // transforms obj into Composite2D *this
  OptimalDecomposite2D optimalDecomp(obj,*this);
  chisquare = optimalDecomp.getOptimalChiSquare();
  // needs more testing before automatically used ...
  // correctCovarianceMatrix();

  history.clear();
  history.setSilent();
  history = obj.history;
  history << optimalDecomp.getHistory();
  history.unsetSilent();
  // joint detection and decomposition flags to form a 16 bit set
  flags.reset();
  const bitset<8>& fitsFlags = obj.flags;
  const bitset<8>& decompFlags = optimalDecomp.getDecompositionFlags();
  for (int i = 0; i < 8; i++) {
    flags[i] = fitsFlags[i];
    flags[8+i] = decompFlags[i];
  }
  // since Composite2D::model ist correctly populated,
  // we currently do not need Composite2D::M anymore (also Composite2D::Mint if present);
  // so we throw it away for saving space
  Composite2D::M.resize(0,0);
  if (ShapeLensConfig::PIXEL_INTEGRATION)
    Composite2D::MInt.resize(0,0);
  Composite2D::changeM = true;
}  


void ShapeletObject::setCoeffs(const CoefficientVector<data_t>& incoeffs) {
  Composite2D::setCoeffs(incoeffs);
  updatePolar = true;
}

void ShapeletObject::setPolarCoeffs(const CoefficientVector<Complex>& newpolarCoeffs) {
  polarCoeffs = newpolarCoeffs;
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::changeM = Composite2D::changeModel = true;
  updatePolar = false;
}
 
void ShapeletObject::setErrors(const CoefficientVector<data_t>& errors) {
  if (errors.size() != coeffs.getNCoeffs()) {
    std::cerr << "ShapeletObject: coefficient errors do not have correct dimension!" << std::endl;
    std::terminate();
  } else {
    if (cov.getRows() != coeffs.getNCoeffs() || cov.getColumns() != coeffs.getNCoeffs())
      cov.resize(coeffs.getNCoeffs(),coeffs.getNCoeffs());
    cov.clear();
    for (unsigned int i=0; i < errors.getNCoeffs(); i++)
      cov(i,i) = errors(i)*errors(i);
  }
}

CoefficientVector<Complex> ShapeletObject::getPolarCoeffs() {
  if (updatePolar) {
    c2p.getPolarCoeffs(coeffs,polarCoeffs);
    updatePolar = false;
  }
  return polarCoeffs;
}

CoefficientVector<Complex> ShapeletObject::getPolarCoeffs() const {
  if (updatePolar) {
    CoefficientVector<Complex> pC;
    PolarTransformation trafo;
    trafo.getPolarCoeffs(coeffs,pC);
    return pC;
  }
  else 
    return polarCoeffs;
}

void ShapeletObject::correctCovarianceMatrix() {
  if (Composite2D::getNMax() > 0) {
    std::vector<NumVector<data_t> > lc;
    lc.push_back(coeffs.getNumVector());
    // 100 slightly modified versions of this
    // changes of the centroid: gaussian, rms of deviation 0.5 pixels
    // changes in beta: gaussian, std = ShapeLensConfig::DELTA_BETA*beta
    data_t beta  = Composite2D::getBeta();
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    for (unsigned int i=0; i< 100; i++) {
      trafo.translate(coeffs,beta,gsl_ran_gaussian(r,M_SQRT1_2*0.5), gsl_ran_gaussian(r,M_SQRT1_2*0.5));
      trafo.rescale(coeffs,beta,beta + gsl_ran_gaussian(r,ShapeLensConfig::DELTA_BETA*beta));
      lc.push_back(coeffs.getNumVector());
      coeffs = lc[0];
    }
    gsl_rng_free (r);

    // build cov matrix from beta/centroid uncertainty:
    // <(lc[i] - av)*(lc[i]-av)^T>
    // use optimal coeffs (=lc[0]) as fiducial value (av)
    NumMatrix<data_t> covU (coeffs.getNCoeffs(),coeffs.getNCoeffs());
    for (unsigned int c=1; c < lc.size(); c++) {
      for (unsigned int i=0; i < coeffs.getNCoeffs(); i++) { 
	for (unsigned int j=0; j < coeffs.getNCoeffs(); j++) {
	  covU(i,j) +=  ((lc[c])(i)-lc[0](i)) * ((lc[c])(j)-lc[0](j));
	}
      }
    }
    covU /= lc.size()-1;
  
    // quadratically add pixel noise contribution from OptimalDecomposite2D
    for (unsigned int i=0; i < cov.getRows(); i++)
      for (unsigned int j=0; j < cov.getColumns(); j++)
	cov(i,j) = sqrt(covU(i,j)*covU(i,j) + cov(i,j)*cov(i,j));
  }
}

CoefficientVector<data_t> ShapeletObject::getErrors() const {
  CoefficientVector<data_t> deltaC(coeffs.getNMax());
  for (unsigned int i=0; i < deltaC.getNCoeffs(); i++) {
    for (unsigned int j=0; j < cov.getColumns(); j++)
      deltaC(i) += cov(i,j);
    deltaC(i) = sqrt(deltaC(i));
  }
  return deltaC;
}

data_t ShapeletObject::getChiSquare() const {
  return chisquare;
}

void ShapeletObject::converge(data_t kappa) {
  trafo.converge(coeffs,kappa,&cov,&history);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::shear(Complex gamma) {
  trafo.shear(coeffs,gamma,&cov,&history);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::flex(Complex F, Complex G) {
  trafo.flex(coeffs,F,G,&cov,&history);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::lens(data_t kappa, Complex gamma, Complex F, Complex G) {
  // lensing is applying convergence, shear and flexion simultaneously
  // since operations work as f' = (1 + x*T) f for any of these transformations
  // we have to do it this way
  // f' = (1+kappa*K) f + (1+gamma*S) f + (1+dGamma*Sij) f - 2*f;
  CoefficientVector<data_t> oldCoeffs = coeffs;
  CoefficientVector<data_t> converged(coeffs.getNMax()), 
    sheared(coeffs.getNMax()), flexed(coeffs.getNMax());
  if (kappa != 0) {
    converge(kappa);
    converged = coeffs;
    Composite2D::setCoeffs(oldCoeffs);
  }

  // apply shear
  shear(gamma);
  sheared = coeffs;
  Composite2D::setCoeffs(oldCoeffs);

  // apply flexion
  flex(F,G);
  flexed = coeffs;

  // sum all up
  if (kappa !=0) {
    flexed += converged;
    oldCoeffs *=2;
  }
  flexed += sheared;
  flexed -= oldCoeffs;
  Composite2D::setCoeffs(flexed);
  Composite2D::changeModel = true;
  updatePolar = true;
} 
  
void ShapeletObject::translate(data_t dx0, data_t dx1) {
  trafo.translate(coeffs,Composite2D::getBeta(),dx0,dx1,&cov,&history);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::rotate(data_t rho) {
  if (updatePolar)
    c2p.getPolarCoeffs(coeffs,polarCoeffs,&cov,&polarCov);
  trafo.rotate(polarCoeffs,rho,&polarCov,&history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs,&polarCov,&cov);
  Composite2D::changeModel = true;
  updatePolar = false;
}

void ShapeletObject::circularize() {
  if (updatePolar)
    c2p.getPolarCoeffs(coeffs,polarCoeffs,&cov,&polarCov);
  trafo.circularize(polarCoeffs,&polarCov,&history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs,&polarCov,&cov);
  Composite2D::changeModel = true;
  updatePolar = false;
}

void ShapeletObject::flipX() {
  if (updatePolar)
    c2p.getPolarCoeffs(coeffs,polarCoeffs,&cov,&polarCov);
  trafo.flipX(polarCoeffs,&polarCov,&history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs,&polarCov,&cov);
  Composite2D::changeModel = true;
  updatePolar = false;
}

void ShapeletObject::brighten(data_t factor) {
  trafo.brighten(coeffs,factor,&cov,&history);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::convolve(const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel) {
  data_t beta = Composite2D::getBeta();
  trafo.convolve(coeffs,beta,kernelCoeffs,beta_kernel,&cov,&history);
  Composite2D::setBeta(beta);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::deconvolve(const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel) {
  data_t beta = Composite2D::getBeta();
  trafo.deconvolve(coeffs,beta,kernelCoeffs,beta_kernel,&cov,&history);
  Composite2D::setBeta(beta);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::rescale(data_t newBeta) {
  data_t beta = Composite2D::getBeta();
  trafo.rescale(coeffs,beta,newBeta,&cov,&history);
  Composite2D::setBeta(beta);
  Composite2D::changeModel = true;
  updatePolar = true;
}

void ShapeletObject::save(string sifFile) const {
  SIFFile sfile(sifFile);
  sfile.save(*this);
}

void ShapeletObject::load(string sifFile, bool preserve_config) {
  SIFFile sfile(sifFile);
  sfile.load(*this, preserve_config);
}

string ShapeletObject::getHistory() const {
  return history.str();
}

void ShapeletObject::setHistory(std::string comment) {
  history.clear();
  history << comment;
}

std::string ShapeletObject::getBaseFilename() const {
  return basefilename;
}

unsigned long ShapeletObject::getObjectID() const {
  return id;
}

const std::bitset<16>& ShapeletObject::getFlags() const {
  return flags;
}

// begin legacy functions //
data_t ShapeletObject::getObjectClassifier() const {
  Property::const_iterator iter = prop.find("classifier");
  data_t d;
  if (iter != prop.end() && iter->second.type() == typeid(d))
    return boost::get<data_t>(iter->second);
  else
    return 0;
}

void ShapeletObject::setName(std::string name) {
  prop["name"] = name;
}

std::string ShapeletObject::getName() const {
  Property::const_iterator iter = prop.find("name");
  std::string s;
  if (iter != prop.end() && iter->second.type() == typeid(s))
    return boost::get<std::string>(iter->second);
  else
    return "";
}

void ShapeletObject::setTag(data_t tag) {
  prop["tag"] = tag;
}

data_t ShapeletObject::getTag() const {
  Property::const_iterator iter = prop.find("tag");
  data_t d;
  if (iter != prop.end() && iter->second.type() == typeid(d))
    return boost::get<data_t>(iter->second);
  else
    return 0;
}
// end legacy function //
