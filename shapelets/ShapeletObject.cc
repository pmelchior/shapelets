#include <ShapeletObject.h>

using namespace std;
typedef complex<double> Complex;

ShapeletObject::ShapeletObject() : Composite2D() {
}

ShapeletObject::ShapeletObject(string sifFile) : Composite2D() {
  // get infos from file
  load(sifFile);
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  trafo = ImageTransformation();
  fits = 0;
  Composite2D::setCoeffs(cartesianCoeffs);
}

ShapeletObject::ShapeletObject(const NumMatrix<double>& incartesianCoeffs, double beta, const Point2D& xcentroid) :
Composite2D() {
  fits = regularized = chisquare = 0;
  fitsFlag = decompFlag = 0;
  R = 0;
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
  text.str("");
}

ShapeletObject::ShapeletObject(const NumMatrix<Complex>& inpolarCoeffs, double beta, const Point2D& xcentroid) :
Composite2D() {
  fits = regularized = chisquare = 0;
  fitsFlag = decompFlag = 0;
  R = 0;
  Composite2D::setCentroid(xcentroid);
  Composite2D::setBeta(beta);
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  cartesianCoeffs = NumMatrix<double> (polarCoeffs.getRows(), polarCoeffs.getColumns());
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  trafo = ImageTransformation();
  Composite2D::setCoeffs(cartesianCoeffs);
  text.str("");
}

// Default values for the shapelet decomposition
unsigned int ShapeletObject::DEFAULTS::NMAX_LOW = 0;
unsigned int ShapeletObject::DEFAULTS::NMAX_HIGH = 100;
double ShapeletObject::DEFAULTS::BETA_LOW = 0;
double ShapeletObject::DEFAULTS::BETA_HIGH = INFINITY;
bool ShapeletObject::DEFAULTS::REGULARIZE = 0;
double ShapeletObject::DEFAULTS::REG_LIMIT = 1e-5;

ShapeletObject::ShapeletObject(const Object& obj) : Composite2D() {//,int nmaxLow, int nmaxHigh, double betaLow, double betaHigh, bool regularize) : Composite2D(){
  fits = 1;
  const Grid& FitsGrid = obj.getGrid();
  fitsFlag = obj.getDetectionFlag();
  // decomposing with given constraits on shapelet decomposition parameters
  OptimalDecomposite2D *optimalDecomp =  new OptimalDecomposite2D(obj, DEFAULTS::NMAX_LOW,DEFAULTS::NMAX_HIGH,DEFAULTS::BETA_LOW,DEFAULTS::BETA_HIGH);
  // use regularization if specified
  if (DEFAULTS::REGULARIZE) {
    regularized = 1;
    R = optimalDecomp->regularize(DEFAULTS::REG_LIMIT);
  } else
    R = regularized = 0;

  optimalDecomp->getShapeletCoeffs(cartesianCoeffs);
  optimalDecomp->getShapeletErrors(errors);
  double beta = optimalDecomp->getOptimalBeta();
  const Point2D& xcentroid = obj.getCentroid();
  chisquare = optimalDecomp->getOptimalChiSquare();
  text << obj.history.getContent();
  text << optimalDecomp->getHistory().getContent();
  history = History(text.str());
  text.str("");
  decompFlag = optimalDecomp->getDecompositionFlag();
  
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getRows());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  trafo = ImageTransformation();
  Composite2D::setBeta(beta);
  Composite2D::setCentroid(xcentroid);
  Composite2D::setCoeffs(cartesianCoeffs);
  Composite2D::setGrid(FitsGrid);
  Composite2D::accessModel() = optimalDecomp->getModel();
  delete optimalDecomp;
}  


void ShapeletObject::setCartesianCoeffs(const NumMatrix<double>& incartesianCoeffs) {
  // put the incoming cartesian coeffs into a triangular matrix of appropriate size
  triangularizeCoeffs(incartesianCoeffs,cartesianCoeffs,0);
  Composite2D::setCoeffs(cartesianCoeffs);
  c2p = PolarTransformation(cartesianCoeffs.getRows()-1);
  polarCoeffs = NumMatrix<Complex> (cartesianCoeffs.getRows(),cartesianCoeffs.getColumns());
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
}

void ShapeletObject::setPolarCoeffs(const NumMatrix<Complex>& inpolarCoeffs) {
  triangularizeCoeffs(inpolarCoeffs,polarCoeffs,0);
  c2p = PolarTransformation(polarCoeffs.getRows()-1);
  cartesianCoeffs = NumMatrix<double> (polarCoeffs.getRows(), polarCoeffs.getColumns());
  Composite2D::setCoeffs(cartesianCoeffs);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
}

const NumMatrix<double>& ShapeletObject::getCartesianCoeffs() {
  return cartesianCoeffs;
}

const NumMatrix<Complex>& ShapeletObject::getPolarCoeffs() {
  return polarCoeffs;
}
 
double ShapeletObject::getDecompositionChiSquare() {
  return chisquare;
}

const NumMatrix<double>& ShapeletObject::getDecompositionErrors() {
  return errors;
}

void ShapeletObject::rotate(double rho) {
  trafo.rotate(polarCoeffs,rho,text);
  history.append(text);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
 }

void ShapeletObject::converge(double kappa) {
  trafo.converge(polarCoeffs,Composite2D::accessBeta(),kappa,text);
  history.append(text);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  // another method would be to change beta without changing the coeffs
  //Composite2D::setBeta((1+kappa)*beta);
}

void ShapeletObject::shear(Complex gamma) {
  // shearing mixes terms with n+-2, so we have extend matrix
  //history << "# Increasing order nmax by 2 for shearing." << endl;
  polarCoeffs.resize(polarCoeffs.getRows()+2,polarCoeffs.getColumns()+2);
  cartesianCoeffs.resize(cartesianCoeffs.getRows()+2,cartesianCoeffs.getColumns()+2);
  c2p.setOrder(polarCoeffs.getRows()-1);
  
  trafo.shear(polarCoeffs,real(gamma),imag(gamma),text);
  history.append(text);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::flex(const NumMatrix<double>& Dgamma) {
  // since flexing mixes term with n+-3 we have to extend matrix dimension by 3
  //history << "# Increasing order nmax by 3 for flexing." << endl;
  polarCoeffs.resize(polarCoeffs.getRows()+3,polarCoeffs.getColumns()+3);
  cartesianCoeffs.resize(cartesianCoeffs.getRows()+3,cartesianCoeffs.getColumns()+3);
  c2p.setOrder(cartesianCoeffs.getRows()-1);
  
  trafo.flex(cartesianCoeffs,Dgamma,text);
  history.append(text);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::lens(double kappa, Complex gamma, Complex F, Complex G) {
  // lensing is applying convergence, shear and flexion simultaneously
  // since operations work as f' = (1 + x*T) f for any of these transformations
  // we have to do it this way
  // f' = (1+kappa*K) f + (1+gamma*S) f + (1+dGamma*Sij) f - 2*f;
  NumMatrix<double> oldCoeffs = cartesianCoeffs;
  NumMatrix<double> converged, sheared, flexed;
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
  NumMatrix<double> dGamma(2,2);
  dGamma(0,0) = 0.5*(real(F) + real(G));
  dGamma(0,1) = 0.5*(imag(G) - imag(F));
  dGamma(1,0) = 0.5*(imag(F) + imag(G));
  dGamma(1,1) = 0.5*(real(F) - real(G));
  flex(dGamma);
  flexed = cartesianCoeffs;

  // we have to extend the coeffs to fit biggest matrix
  // which is the flexed one
  oldCoeffs.resize(oldCoeffs.getRows()+3,oldCoeffs.getColumns()+3);
  converged.resize(converged.getRows()+3,converged.getColumns()+3);
  sheared.resize(sheared.getRows()+1,sheared.getColumns()+1);

  if (kappa !=0) {
    flexed += converged;
    oldCoeffs *=2;
  }
  flexed += sheared;
  flexed -= oldCoeffs;
  setCartesianCoeffs(flexed);
} 
  
void ShapeletObject::translate(double dx0, double dx1) {
  trafo.translate(cartesianCoeffs,Composite2D::getBeta(),dx0,dx1,text);
  history.append(text);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::circularize() {
  trafo.circularize(polarCoeffs,text);
  history.append(text);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::flipX() {
  trafo.flipX(polarCoeffs,text);
  history.append(text);
  c2p.getCartesianCoeffs(polarCoeffs,cartesianCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::brighten(double factor) {
  trafo.brighten(cartesianCoeffs,polarCoeffs,factor,text);
  history.append(text);
  Composite2D::setCoeffs(cartesianCoeffs);
}

void ShapeletObject::convolve(const NumMatrix<double>& KernelCoeffs, double beta_kernel) {
  // first triangularize coeffs
  NumMatrix<double> triKernelCoeffs;
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.convolve(cartesianCoeffs,Composite2D::accessBeta(),triKernelCoeffs,beta_kernel,text);
  history.append(text);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  //  Composite2D::setBeta(beta);
}

void ShapeletObject::deconvolve(const NumMatrix<double>& KernelCoeffs, double beta_kernel) {
  // first triangularize coeffs
  NumMatrix<double> triKernelCoeffs;
  triangularizeCoeffs(KernelCoeffs,triKernelCoeffs,0);
  trafo.deconvolve(cartesianCoeffs,Composite2D::accessBeta(),triKernelCoeffs,beta_kernel,text);
  history.append(text);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  //Composite2D::setBeta(beta);
}

void ShapeletObject::rescale(double newBeta) {
  trafo.rescale(cartesianCoeffs,Composite2D::getBeta(),newBeta,text);
  history.append(text);
  c2p.getPolarCoeffs(cartesianCoeffs,polarCoeffs);
  Composite2D::setCoeffs(cartesianCoeffs);
  Composite2D::setBeta(newBeta);
}

void ShapeletObject::getProfile(const Point2D& start, NumVector<double>& values, int axsize) {
  const Point2D& xcentroid = Composite2D::getCentroid();
  const Point2D end = Point2D((xcentroid(0) - start(0))*2 + start(0), (xcentroid(1) - start(1))*2 + start(1));
  Profile profile(start,end);
  profile.calculate(values,axsize);
}

void ShapeletObject::getProfile(const Point2D& start, const Point2D& end, NumVector<double>& values, int axsize) {
   Profile profile (start,end);
   profile.calculate(values,axsize);
}

void ShapeletObject::save(string sifFile) {
  SIFFile sfile(sifFile);
  sfile.save(history.getContent(),cartesianCoeffs,errors,Composite2D::getGrid(),Composite2D::getBeta(),Composite2D::getCentroid(),chisquare,fitsFlag,decompFlag,regularized,R);
}

void ShapeletObject::load(string sifFile) {
  std::string historyString;
  Grid grid;
  double beta;
  Point2D xcentroid;
  SIFFile sfile(sifFile);
  sfile.load(historyString,cartesianCoeffs,errors,grid,beta,xcentroid,chisquare,fitsFlag,decompFlag,regularized,R);
  history = History(historyString);
  Composite2D::setGrid(grid);
  Composite2D::setBeta(beta);
  Composite2D::setCentroid(xcentroid);
  Composite2D::setCoeffs(cartesianCoeffs);
}

string ShapeletObject::getHistory() {
  return history.getContent();
}

void ShapeletObject::setHistory(std::string comment) {
  history = History(comment);
}
