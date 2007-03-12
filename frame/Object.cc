#include <Object.h>
#include <sstream>
#include <gsl/gsl_math.h>

Object::Object(unsigned int inid) {
  id = inid;
  flag = 0;
  blend = s_g = -1;
  flux = centroid(0) = centroid(1) = -1;
}

unsigned int Object::getID() const {
  return id;
}

const Grid& Object::getGrid() const {
  return grid;
}

Grid& Object::accessGrid() {
  return grid;
}

const NumVector<double>& Object::getData() const {
  return data;
}

NumVector<double>& Object::accessData() {
  return data;
}

const NumVector<double>& Object::getBackground() const {
  return bg_mean;
}

NumVector<double>& Object::accessBackground() {
  return bg_mean;
}

const NumVector<double>& Object::getBackgroundRMS() const {
  return bg_rms;
}

NumVector<double>& Object::accessBackgroundRMS() {
  return bg_rms;
}

int Object::getSize(bool direction) const {
  return grid.getSize(direction);
}

const Point2D& Object::getCentroid() const {
  return centroid;
}

void Object::setCentroid(const Point2D& xc) {
  centroid = xc;
  std::ostringstream text;
  text << "# Centroid set to ("<< centroid(0) << "/" << centroid(1) <<  ") by user." << std::endl;
  history.append(text);
}

double Object::getFlux() const {
  return flux;
}

void Object::setFlux(double F) {
  flux = F;
  std::ostringstream text;
  text << "# Flux set to " << flux << " by user." <<std::endl;
  history.append(text);
}

NumMatrix<double> Object::get2ndBrightnessMoments() {
  NumMatrix<double> Q(2,2);
  for (int i=0; i< grid.size(); i++) {
    Q(0,0) += gsl_pow_2(grid(i,0)-centroid(0)) * data(i);
    Q(0,1) += (grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i);
    Q(1,1) += gsl_pow_2(grid(i,1)-centroid(1)) * data(i);
  }
  Q(0,0) /= flux;
  Q(0,1) /= flux;
  Q(1,0) = Q(0,1);
  Q(1,1) /= flux;
  return Q;
}

// ellipticity as defined in Bartelmann & Schneider (2001)
complex<double> Object::getEllipticity() {
  NumMatrix<double> Q = get2ndBrightnessMoments();
  complex<double> I(0,1);
  complex<double> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<double> epsilon = (Q11 - Q22 + 2.*I*Q12)/(Q11+Q22 + 2.*sqrt(Q11*Q22-Q12*Q12));
  return epsilon;
  // axis ratio
  //double e = abs(epsilon);
  //double r = (1-e)/(1+e);
}

void Object::computeFluxCentroid() {
  flux = 0;
  centroid(0) = centroid(1) = 0;
  for (int i=0; i< grid.size(); i++) {
    flux += data(i);
    centroid(0) += grid(i,0) * data(i);
    centroid(1) += grid(i,1) * data(i);
  }
  centroid(0) /= flux;
  centroid(1) /= flux;
  
  std::ostringstream text;
  text << "# Flux = " << flux << ", Centroid = ("<< centroid(0) << "/" << centroid(1);
  text <<  ")." <<std::endl;
  history.append(text);
}

unsigned short Object::getDetectionFlag() const {
  return flag;
}

void Object::setDetectionFlag(unsigned short inflag) {
  flag = inflag;
}

double Object::getNoiseMean() const {
  return noise_mean;
}

double Object::getNoiseRMS() const {
  return noise_rms;
}

void Object::setNoiseMeanRMS(double mean, double rms) {
  noise_mean = mean;
  noise_rms = rms;
}

void Object::setNoiseModel(std::string innoisemodel) {
  noisemodel = innoisemodel;
}

std::string Object::getNoiseModel() const {
  return noisemodel;
}

double Object::getBlendingProbability() const {
  return blend;
}

void Object::setBlendingProbability(double b) {
  blend = b;
}

double Object::getStarGalaxyProbability() const {
  return s_g;
}

void Object::setStarGalaxyProbability(double sg) {
  s_g = sg;
}

void Object::setBaseFilename(std::string filename) {
  basefilename = filename;
}

void Object::save(std::string fitsfile) {
  std::map<std::string, std::string> keywords;
  keywords["BASEFILE"] = basefilename;
  std::ostringstream value;
  value << grid.getStartPosition(0);
  keywords["XMIN"] = value.str();
  value.str("");
  value << grid.getStopPosition(0);
  keywords["XMAX"] = value.str();
  value.str("");
  value << grid.getStartPosition(1);
  keywords["YMIN"] = value.str();
  value.str("");
  value << grid.getStopPosition(1);
  keywords["YMAX"] = value.str();
  value.str("");
  value << id;
  keywords["ID"] = value.str();
  value.str("");
  value << noise_mean;
  keywords["BG_MEAN"] = value.str();
  value.str("");
  value << noise_rms;
  keywords["BG_RMS"] = value.str();
  keywords["NOISE"] = noisemodel;
  value.str("");
  value << flag;
  keywords["FLAG"] = value.str();
  value.str("");
  value << blend;
  keywords["BLEND"] = value.str();
  value.str("");
  value << s_g;
  keywords["S_G"] = value.str();
  writeFITSFile(fitsfile,grid,data,keywords);

  //if background maps are provided, save them too
  keywords.clear();
  if (bg_mean.size() != 0)
    addFITSExtension(fitsfile,"BG_MEAN",grid,bg_mean,keywords);
  if (bg_rms.size() != 0)
    addFITSExtension(fitsfile,"BG_RMS",grid,bg_rms,keywords);
}
