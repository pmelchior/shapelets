#include <shapelets/ShapeletObject.h>
#include <shapelets/SIFFile.h>
#include <stdio.h>

using namespace std;
typedef complex<data_t> Complex;

ShapeletObject::ShapeletObject() : 
Composite2D(), coeffs(Composite2D::coeffs) {
   tag = classifier = chisquare = R = noise_mean = noise_rms = nr = id = fits = 0;
  flags.reset();
  name = basefilename = "";
  unreg = NULL;
}

ShapeletObject::ShapeletObject(const ShapeletObject& source) : 
Composite2D(), coeffs(Composite2D::coeffs) {
  ShapeletObject::operator=(source);
}

ShapeletObject& ShapeletObject::operator=(const ShapeletObject& source) {
  Composite2D::operator=(source);
  // copy every variable of this class
  // base class copy constructor already called
  errors = source.errors;
  polarCoeffs = source.polarCoeffs;
  c2p = source.c2p;
  trafo = source.trafo;
  tag = source.tag;
  classifier = source.classifier;
  chisquare = source.chisquare;
  R = source.R;
  noise_mean = source.noise_mean;
  noise_rms = source.noise_rms;
  nr = source.nr;
  id = source.id;
  fits = source.fits;
  flags = source.flags;
  name = source.name;
  basefilename = source.basefilename;
  history = source.history;
  // important: to avoid double deletion, unreg must not be copied!
  unreg = NULL;
  return *this;
}

ShapeletObject::~ShapeletObject() {
  if (unreg != NULL) 
    delete unreg;
}

ShapeletObject::ShapeletObject(string sifFile, bool preserve_config) : 
Composite2D(), coeffs(Composite2D::coeffs) {
  unreg = NULL;
  // get infos from file
  load(sifFile,preserve_config);
  c2p = PolarTransformation(coeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (coeffs.getRows(),coeffs.getColumns());
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  trafo = ImageTransformation();
  fits = 0;
}

ShapeletObject::ShapeletObject(const NumMatrix<data_t>& incoeffs, data_t beta, const Point2D& xcentroid, const Grid& grid) :
Composite2D(), coeffs(Composite2D::coeffs) {
  tag = classifier = chisquare = R = noise_mean = noise_rms = nr = id = fits = 0;
  flags.reset();
  name = basefilename = "";
  unreg = NULL;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  Composite2D::setGrid(grid);
  // put the incoming cartesian coeffs into a triangular matrix of appropriate size
  triangularizeCoeffs(incoeffs,coeffs,0);
  c2p = PolarTransformation(coeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (coeffs.getRows(),coeffs.getColumns());
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  trafo = ImageTransformation();
}

ShapeletObject::ShapeletObject(const NumMatrix<Complex>& inpolarCoeffs, data_t beta, const Point2D& xcentroid, const Grid& grid) :
Composite2D(), coeffs(Composite2D::coeffs) {
  tag = classifier = chisquare = R = noise_mean = noise_rms = nr = id = fits = 0;
  flags.reset();
  name = basefilename = "";
  unreg = NULL;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  Composite2D::setGrid(grid);
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  coeffs = NumMatrix<data_t> (polarCoeffs.getRows(), polarCoeffs.getColumns());
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  trafo = ImageTransformation();
  Composite2D::setCoeffs(coeffs);
}

ShapeletObject::ShapeletObject(const Object& obj) : 
Composite2D(), coeffs(Composite2D::coeffs) {
  fits = 1;
  tag = R = 0;
  name = "";
  unreg = NULL;
  const Grid& FitsGrid = obj.getGrid();
  const Point2D& xcentroid = obj.getCentroid();
  noise_mean = obj.getNoiseMean();
  noise_rms = obj.getNoiseRMS();
  id = obj.getID();
  nr = obj.getNumber();
  classifier = obj.getClassifier();
  basefilename = obj.getBaseFilename();
  // decomposing with given constraits on shapelet decomposition parameters
  OptimalDecomposite2D optimalDecomp(obj, ShapeLensConfig::NMAX_LOW,ShapeLensConfig::NMAX_HIGH,ShapeLensConfig::BETA_LOW,ShapeLensConfig::BETA_HIGH);

  // if set, save the unregularized model to sif file with given name
  if (ShapeLensConfig::REGULARIZE && ShapeLensConfig::SAVE_UNREG) {
    // first get all necessary data for model
    optimalDecomp.getShapeletCoeffs(coeffs);
    optimalDecomp.getShapeletErrors(errors);
    chisquare = optimalDecomp.getOptimalChiSquare();
    history.setSilent();
    history << obj.history.str();
    history << optimalDecomp.getHistory().str();
    history.unsetSilent();
    // joing detection and decomposition flags to form a 16 bit set
    const bitset<8>& fitsFlags = obj.getDetectionFlags();
    const bitset<8>& decompFlags = optimalDecomp.getDecompositionFlags();
    for (int i = 0; i < 8; i++) {
      flags[i] = fitsFlags[i];
      flags[8+i] = decompFlags[i];
    }
    c2p = PolarTransformation(coeffs.getRows()-1);
    polarCoeffs = NumMatrix<Complex> (coeffs.getRows(),coeffs.getRows());
    c2p.getPolarCoeffs(coeffs,polarCoeffs);
    trafo = ImageTransformation();
    Composite2D::setBeta(optimalDecomp.getOptimalBeta());
    Composite2D::setCentroid(xcentroid);
    Composite2D::setGrid(FitsGrid);
    // create copy of *this as new ShapeletObject entity
    name = "UNREG";
    unreg = new ShapeletObject(*this);
    name = "";
  }

  // use regularization if specified
  if (ShapeLensConfig::REGULARIZE)
    R = optimalDecomp.regularize(ShapeLensConfig::REG_LIMIT);

  optimalDecomp.getShapeletCoeffs(coeffs);
  optimalDecomp.getShapeletErrors(errors);
  chisquare = optimalDecomp.getOptimalChiSquare();
  history.clear();
  history.setSilent();
  history << obj.history.str();
  history << optimalDecomp.getHistory().str();
  history.unsetSilent();
  // joing detection and decomposition flags to form a 16 bit set
  flags.reset();
  const bitset<8>& fitsFlags = obj.getDetectionFlags();
  const bitset<8>& decompFlags = optimalDecomp.getDecompositionFlags();
  for (int i = 0; i < 8; i++) {
    flags[i] = fitsFlags[i];
    flags[8+i] = decompFlags[i];
  }

  c2p = PolarTransformation(coeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (coeffs.getRows(),coeffs.getRows());
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  trafo = ImageTransformation();
  Composite2D::setBeta(optimalDecomp.getOptimalBeta());
  Composite2D::setCentroid(xcentroid);
  Composite2D::setGrid(FitsGrid);
  Composite2D::accessModel() = optimalDecomp.getModel();
}  


void ShapeletObject::setCoeffs(const NumMatrix<data_t>& incoeffs) {
  // put the incoming cartesian coeffs into a triangular matrix of appropriate size
  triangularizeCoeffs(incoeffs,coeffs,0);
  Composite2D::change = 1;
  c2p = PolarTransformation(coeffs.getRows()-1);
  polarCoeffs.resize(coeffs.getRows(),coeffs.getColumns());
  polarCoeffs.clear();
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
}

void ShapeletObject::setCoeffErrors(const NumMatrix<data_t>& newerrors) {
  if (coeffs.getRows() == newerrors.getRows() && coeffs.getColumns() == newerrors.getColumns())
    errors = newerrors;
  else {
    std::cerr << "ShapeletObject: errors given do not have correct dimensions!" << std::endl;
    std::terminate();
  }
}
void ShapeletObject::setPolarCoeffs(const NumMatrix<Complex>& inpolarCoeffs) {
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  coeffs.resize(polarCoeffs.getRows(), polarCoeffs.getColumns());
  coeffs.clear();
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::change = 1;
}

const NumMatrix<Complex>& ShapeletObject::getPolarCoeffs() const {
  return polarCoeffs;
}
 
data_t ShapeletObject::getDecompositionChiSquare() const {
  return chisquare;
}

const NumMatrix<data_t>& ShapeletObject::getDecompositionErrors() const {
  return errors;
}

void ShapeletObject::rotate(data_t rho) {
  trafo.rotate(polarCoeffs,rho,history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::change = 1;
 }

void ShapeletObject::converge(data_t kappa) {
  data_t beta;
  trafo.converge(polarCoeffs,beta,kappa,history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::setBeta(beta);
  Composite2D::change = 1;
  // another method would be to change beta without changing the coeffs
  //Composite2D::setBeta((1+kappa)*beta);
}

void ShapeletObject::shear(Complex gamma) {
  // shearing mixes terms with n+-2, so we have extend matrix
  //history << "# Increasing order nmax by 2 for shearing." << endl;
  polarCoeffs.resize_clear(polarCoeffs.getRows()+2,polarCoeffs.getColumns()+2);
  coeffs.resize_clear(coeffs.getRows()+2,coeffs.getColumns()+2);
  c2p.setOrder(polarCoeffs.getRows()-1);
  
  trafo.shear(polarCoeffs,real(gamma),imag(gamma),history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::change = 1;
}

void ShapeletObject::flex(const NumMatrix<data_t>& Dgamma) {
  // since flexing mixes term with n+-3 we have to extend matrix dimension by 3
  //history << "# Increasing order nmax by 3 for flexing." << endl;
  polarCoeffs.resize_clear(polarCoeffs.getRows()+3,polarCoeffs.getColumns()+3);
  coeffs.resize_clear(coeffs.getRows()+3,coeffs.getColumns()+3);
  c2p.setOrder(coeffs.getRows()-1);
  
  trafo.flex(coeffs,Dgamma,history);
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  Composite2D::change = 1;
}

void ShapeletObject::lens(data_t kappa, Complex gamma, Complex F, Complex G) {
  // lensing is applying convergence, shear and flexion simultaneously
  // since operations work as f' = (1 + x*T) f for any of these transformations
  // we have to do it this way
  // f' = (1+kappa*K) f + (1+gamma*S) f + (1+dGamma*Sij) f - 2*f;
  NumMatrix<data_t> oldCoeffs = coeffs;
  NumMatrix<data_t> converged, sheared, flexed;
  if (kappa != 0) {
    converge(kappa);
    converged = coeffs;
    setCoeffs(oldCoeffs);
  }
  // apply shear
  shear(gamma);
  sheared = coeffs;
  setCoeffs(oldCoeffs);
  // apply flexion
  NumMatrix<data_t> dGamma(2,2);
  dGamma(0,0) = 0.5*(real(F) + real(G));
  dGamma(0,1) = 0.5*(imag(G) - imag(F));
  dGamma(1,0) = 0.5*(imag(F) + imag(G));
  dGamma(1,1) = 0.5*(real(F) - real(G));
  flex(dGamma);
  flexed = coeffs;

  // we have to extend the coeffs to fit biggest matrix
  // which is the flexed one
  oldCoeffs.resize_clear(oldCoeffs.getRows()+3,oldCoeffs.getColumns()+3);
  converged.resize_clear(converged.getRows()+3,converged.getColumns()+3);
  sheared.resize_clear(sheared.getRows()+1,sheared.getColumns()+1);

  if (kappa !=0) {
    flexed += converged;
    oldCoeffs *=2;
  }
  flexed += sheared;
  flexed -= oldCoeffs;
  setCoeffs(flexed);
} 
  
void ShapeletObject::translate(data_t dx0, data_t dx1) {
  trafo.translate(coeffs,Composite2D::getBeta(),dx0,dx1,history);
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  Composite2D::change = 1;
}

void ShapeletObject::circularize() {
  trafo.circularize(polarCoeffs,history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::change = 1;
}

void ShapeletObject::flipX() {
  trafo.flipX(polarCoeffs,history);
  c2p.getCartesianCoeffs(polarCoeffs,coeffs);
  Composite2D::change = 1;
}

void ShapeletObject::brighten(data_t factor) {
  trafo.brighten(coeffs,polarCoeffs,factor,history);
  Composite2D::change = 1;
}

void ShapeletObject::convolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel) {
  // first triangularize coeffs
  NumMatrix<data_t> triKernelCoeffs;
  data_t beta = Composite2D::getBeta();
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.convolve(coeffs,beta,triKernelCoeffs,beta_kernel,history);
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  Composite2D::change = 1;
  Composite2D::setBeta(beta);
}

void ShapeletObject::deconvolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel) {
  // first triangularize coeffs
  NumMatrix<data_t> triKernelCoeffs;
  data_t beta = Composite2D::getBeta();
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.deconvolve(coeffs,beta,triKernelCoeffs,beta_kernel,history);
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  Composite2D::change = 1;
  Composite2D::setBeta(beta);
}

void ShapeletObject::rescale(data_t newBeta) {
  // since the rescaling transformation is only valid for infinite expansions
  // raising the order here helps to obtain proper results
  // But: The selected nmax is arbitrary; therefore we keep nmax fixed.
  trafo.rescale(coeffs,Composite2D::getBeta(),newBeta,history);
  c2p.getPolarCoeffs(coeffs,polarCoeffs);
  Composite2D::change = 1;
}

void ShapeletObject::save(string sifFile) const {
  SIFFile sfile(sifFile);
  // only store actual ShapeletObject.
  if (unreg == NULL)
    sfile.save(*this);
  // since there is also the unregularized object and we have to save it, too
  else
    sfile.save(*this, *unreg);
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

data_t ShapeletObject::getRegularizationR() const {
  return R;
}
std::string ShapeletObject::getBaseFilename() const {
  return basefilename;
}

unsigned long ShapeletObject::getObjectID() const {
  return id;
}

unsigned long ShapeletObject::getObjectNumber() const {
  return nr;
}

data_t ShapeletObject::getObjectClassifier() const {
  return classifier;
}

const std::bitset<16>& ShapeletObject::getFlags() const {
  return flags;
}

data_t ShapeletObject::getNoiseMean() const {
  return noise_mean;
}

data_t ShapeletObject::getNoiseRMS() const {
  return noise_rms;
}

void ShapeletObject::setName(std::string inname) {
  name = inname;
}

std::string ShapeletObject::getName() const {
  return name;
}

void ShapeletObject::setTag(data_t intag) {
  tag = intag;
}

data_t ShapeletObject::getTag() const {
  return tag;
}
  
