#include <shapelets/ShapeletObject.h>
#include <shapelets/SIFFile.h>
#include <stdio.h>

using namespace std;
typedef complex<data_t> Complex;

ShapeletObject::ShapeletObject() : Composite2D() {
}

ShapeletObject::ShapeletObject(string sifFile, bool preserve_config) : Composite2D() {
  // get infos from file
  load(sifFile,preserve_config);
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  trafo = ImageTransformation();
  fits = 0;
  Composite2D::setCoeffs(cartesianCoeffs);
}

ShapeletObject::ShapeletObject(const NumMatrix<data_t>& incartesianCoeffs, data_t beta, const Point2D& xcentroid) :
Composite2D() {
  chisquare = R = noise_mean = noise_rms = nr = id = flags = fits = 0;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  // put the incoming cartesian coeffs into a triangular matrix of appropriate size
  triangularizeCoeffs(incartesianCoeffs,cartesianCoeffs,0);
  // now set the triangular coeffs
  Composite2D::setCoeffs(cartesianCoeffs);
 
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  trafo = ImageTransformation();
}

ShapeletObject::ShapeletObject(const NumMatrix<Complex>& inpolarCoeffs, data_t beta, const Point2D& xcentroid) :
Composite2D() {
  chisquare = R = noise_mean = noise_rms = nr = id = flags = fits = 0;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  cartesianCoeffs = NumMatrix<data_t> (polarCoeffs.getRows(), polarCoeffs.getColumns());
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  trafo = ImageTransformation();
  Composite2D::setCoeffs(cartesianCoeffs);
}

ShapeletObject::ShapeletObject(const Object& obj) : Composite2D() {
  fits = 1;
  R = 0;
  const Grid& FitsGrid = obj.getGrid();
  const Point2D& xcentroid = obj.getCentroid();
  unsigned int fitsFlag = obj.getDetectionFlag(), decompFlag = 0;
  data_t beta;
  noise_mean = obj.getNoiseMean();
  noise_rms = obj.getNoiseRMS();
  id = obj.getID();
  nr = obj.getNumber();
  basefilename = obj.getBaseFilename();
  // decomposing with given constraits on shapelet decomposition parameters
  OptimalDecomposite2D optimalDecomp(obj, ShapeLensConfig::NMAX_LOW,ShapeLensConfig::NMAX_HIGH,ShapeLensConfig::BETA_LOW,ShapeLensConfig::BETA_HIGH);

  // if set, save the unregularized model to sif file with given name
  if (ShapeLensConfig::REGULARIZE && ShapeLensConfig::SAVE_UNREG) {
    // first get all necessary data for model
    optimalDecomp.getShapeletCoeffs(cartesianCoeffs);
    optimalDecomp.getShapeletErrors(errors);
    beta = optimalDecomp.getOptimalBeta();
    chisquare = optimalDecomp.getOptimalChiSquare();
    history.setSilent();
    history << obj.history.str();
    history << optimalDecomp.getHistory().str();
    history.unsetSilent();
    decompFlag = optimalDecomp.getDecompositionFlag();
    if (decompFlag > 0)
      flags = decompFlag+256+fitsFlag;
    else
      flags = fitsFlag;
    // save temporary file here
    // TODO: use filesystem calls from boost to move is when save() is called
    SIFFile sfile("obj_unreg.sif");
    ShapeLensConfig::REGULARIZE = 0;
    sfile.save(*this);
    ShapeLensConfig::REGULARIZE = 1;
  }

  // use regularization if specified
  if (ShapeLensConfig::REGULARIZE)
    R = optimalDecomp.regularize(ShapeLensConfig::REG_LIMIT);

  optimalDecomp.getShapeletCoeffs(cartesianCoeffs);
  optimalDecomp.getShapeletErrors(errors);
  beta = optimalDecomp.getOptimalBeta();
  chisquare = optimalDecomp.getOptimalChiSquare();
  history.clear();
  history.setSilent();
  history << obj.history.str();
  history << optimalDecomp.getHistory().str();
  history.unsetSilent();
  decompFlag = optimalDecomp.getDecompositionFlag();
  if (decompFlag > 0)
    flags = decompFlag+256+fitsFlag;
  else
    flags = fitsFlag;

  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getRows());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  trafo = ImageTransformation();
  Composite2D::setBeta(beta);
  Composite2D::setCentroid(xcentroid);
  Composite2D::setCoeffs(cartesianCoeffs);
  Composite2D::setGrid(FitsGrid);
  Composite2D::accessModel() = optimalDecomp.getModel();
}  


void ShapeletObject::setCartesianCoeffs(const NumMatrix<data_t>& incartesianCoeffs) {
  // put the incoming cartesian coeffs into a triangular matrix of appropriate size
  triangularizeCoeffs(incartesianCoeffs,cartesianCoeffs,0);
  Composite2D::setCoeffs(cartesianCoeffs);
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
}

void ShapeletObject::setCartesianCoeffErrors(const NumMatrix<data_t>& newerrors) {
  if (cartesianCoeffs.getRows() == newerrors.getRows() && cartesianCoeffs.getColumns() == newerrors.getColumns())
    errors = newerrors;
  else {
    std::cout << "ShapeletObject: errors given do not have correct dimensions!" << std::endl;
    std::terminate();
  }
}
void ShapeletObject::setPolarCoeffs(const NumMatrix<Complex>& inpolarCoeffs) {
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  cartesianCoeffs = NumMatrix<data_t> (polarCoeffs.getRows(), polarCoeffs.getColumns());
  Composite2D::setCoeffs(cartesianCoeffs);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
}

const NumMatrix<data_t>& ShapeletObject::getCartesianCoeffs() const {
  return cartesianCoeffs;
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
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
 }

void ShapeletObject::converge(data_t kappa) {
  trafo.converge(polarCoeffs,Composite2D::accessBeta(),kappa,history);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  // another method would be to change beta without changing the coeffs
  //Composite2D::setBeta((1+kappa)*beta);
}

void ShapeletObject::shear(Complex gamma) {
  // shearing mixes terms with n+-2, so we have extend matrix
  //history << "# Increasing order nmax by 2 for shearing." << endl;
  polarCoeffs.resize_clear(polarCoeffs.getRows()+2,polarCoeffs.getColumns()+2);
  cartesianCoeffs.resize_clear(cartesianCoeffs.getRows()+2,cartesianCoeffs.getColumns()+2);
  c2p.setOrder(polarCoeffs.getRows()-1);
  
  trafo.shear(polarCoeffs,real(gamma),imag(gamma),history);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::flex(const NumMatrix<data_t>& Dgamma) {
  // since flexing mixes term with n+-3 we have to extend matrix dimension by 3
  //history << "# Increasing order nmax by 3 for flexing." << endl;
  polarCoeffs.resize_clear(polarCoeffs.getRows()+3,polarCoeffs.getColumns()+3);
  cartesianCoeffs.resize_clear(cartesianCoeffs.getRows()+3,cartesianCoeffs.getColumns()+3);
  c2p.setOrder(cartesianCoeffs.getRows()-1);
  
  trafo.flex(cartesianCoeffs,Dgamma,history);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::lens(data_t kappa, Complex gamma, Complex F, Complex G) {
  // lensing is applying convergence, shear and flexion simultaneously
  // since operations work as f' = (1 + x*T) f for any of these transformations
  // we have to do it this way
  // f' = (1+kappa*K) f + (1+gamma*S) f + (1+dGamma*Sij) f - 2*f;
  NumMatrix<data_t> oldCoeffs = cartesianCoeffs;
  NumMatrix<data_t> converged, sheared, flexed;
  if (kappa != 0) {
    converge(kappa);
    converged = cartesianCoeffs;
    setCartesianCoeffs(oldCoeffs);
  }
  // apply shear
  shear(gamma);
  sheared = cartesianCoeffs;
  setCartesianCoeffs(oldCoeffs);
  // apply flexion
  NumMatrix<data_t> dGamma(2,2);
  dGamma(0,0) = 0.5*(real(F) + real(G));
  dGamma(0,1) = 0.5*(imag(G) - imag(F));
  dGamma(1,0) = 0.5*(imag(F) + imag(G));
  dGamma(1,1) = 0.5*(real(F) - real(G));
  flex(dGamma);
  flexed = cartesianCoeffs;

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
  setCartesianCoeffs(flexed);
} 
  
void ShapeletObject::translate(data_t dx0, data_t dx1) {
  trafo.translate(cartesianCoeffs,Composite2D::getBeta(),dx0,dx1,history);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::circularize() {
  trafo.circularize(polarCoeffs,history);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::flipX() {
  trafo.flipX(polarCoeffs,history);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::brighten(data_t factor) {
  trafo.brighten(cartesianCoeffs,polarCoeffs,factor,history);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::convolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel) {
  // first triangularize coeffs
  NumMatrix<data_t> triKernelCoeffs;
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.convolve(cartesianCoeffs,Composite2D::accessBeta(),triKernelCoeffs,beta_kernel,history);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  //  Composite2D::setBeta(beta);
}

void ShapeletObject::deconvolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel) {
  // first triangularize coeffs
  NumMatrix<data_t> triKernelCoeffs;
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.deconvolve(cartesianCoeffs,Composite2D::accessBeta(),triKernelCoeffs,beta_kernel,history);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  //Composite2D::setBeta(beta);
}

void ShapeletObject::rescale(data_t newBeta) {
  // since the rescaling transformation is only valid for infinite expansions
  // raising the order here helps to obtain proper results
  // But: The selected nmax is arbitrary; therefore we keep nmax fixed.
  trafo.rescale(cartesianCoeffs,Composite2D::getBeta(),newBeta,history);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::save(string sifFile) const {
  SIFFile sfile(sifFile);
  sfile.save(*this);
  // if the sif file for the unregularized object has been saved:
  // append it to the current sif file
  if (ShapeLensConfig::REGULARIZE && ShapeLensConfig::SAVE_UNREG) {
    int status = 0;
    fitsfile* infptr = openFITSFile("obj_unreg.sif",0);
    fitsfile* outfptr = openFITSFile(sifFile,1);
    fits_copy_file(infptr,outfptr,1,1,1,&status);

    // set EXTNAME ="UNREG" and "UNREG_ERRORS" for new HDUs
    int numHDUs;
    fits_get_num_hdus(outfptr,&numHDUs,&status);
    // move the the last two HDUs
    fits_movabs_hdu(outfptr,numHDUs-1,IMAGE_HDU,&status);
    status = updateFITSKeywordString(outfptr,"EXTNAME","UNREG");
    fits_movabs_hdu(outfptr,numHDUs,IMAGE_HDU,&status);
    fits_set_hdustruc(outfptr, &status);
    status = updateFITSKeywordString(outfptr,"EXTNAME","UNREG_ERRORS");
    fits_set_hdustruc(outfptr, &status);
    fits_close_file(outfptr, &status);
    fits_close_file(infptr, &status);
    
    // remove temporary sif file
    remove("obj_unreg.sif");
  }
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

unsigned int ShapeletObject::getObjectID() const {
  return id;
}

data_t ShapeletObject::getObjectNumber() const {
  return nr;
}

std::bitset<16> ShapeletObject::getFlags() const {
  return std::bitset<16>(flags);
}

data_t ShapeletObject::getNoiseMean() const {
  return noise_mean;
}

data_t ShapeletObject::getNoiseRMS() const {
  return noise_rms;
}
