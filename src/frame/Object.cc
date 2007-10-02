#include <frame/Object.h>
#include <sstream>
#include <gsl/gsl_math.h>
#include <fitsio.h>

Object::Object(unsigned int inid) : Image<data_t>(), segMap(*this) {
  id = inid;
  flag = 0;
  blend = s_g = -1;
  flux = centroid(0) = centroid(1) = -1;
  number = 0;
}

Object::Object(std::string objfile) : Image<data_t>(), segMap(*this) {
  fitsfile *fptr;
  int status, nkeys, keypos, hdutype;
  char card[FLEN_CARD];
  char comment[FLEN_CARD];
  status = 0; 

  history << "# Loading object from Fits file " << objfile << ":" << std::endl;
  fits_open_file(&fptr, objfile.c_str(), READONLY, &status);

  // open pHDU and read header
  fits_movabs_hdu(fptr, 1, &hdutype, &status);
  // recover object information from header keywords
  status = readFITSKeywordString(fptr,"BASEFILE",basefilename);
  status = readFITSKeyword(fptr,"ID",id);
  status = readFITSKeyword(fptr,"NUMBER",number);
  data_t xmin,xmax,ymin,ymax;
  status = readFITSKeyword(fptr,"XMIN",xmin);
  status = readFITSKeyword(fptr,"XMAX",xmax);
  status = readFITSKeyword(fptr,"YMIN",ymin);
  status = readFITSKeyword(fptr,"YMAX",ymax);
  Image<data_t>::accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);
  complex<data_t> xc;
  status = readFITSKeyword(fptr,"CENTROID",xc);
  centroid(0) = real(xc);
  centroid(1) = imag(xc);
  status = readFITSKeyword(fptr,"BG_MEAN",noise_mean);
  status = readFITSKeyword(fptr,"BG_RMS",noise_rms);
  status = readFITSKeyword(fptr,"FLAG",flag);
  status = readFITSKeyword(fptr,"BLEND",blend);
  status = readFITSKeyword(fptr,"S_G",s_g);
  
  // read history
  std::string hstr;
  readFITSKeyCards(fptr,"HISTORY",hstr);
  history.clear();
  bool verb = history.getVerbosity();
  history.setVerbosity(0);
  history << hstr;
  history.setVerbosity(verb);

  // if status is not 0 we have faced an error
  if (status != 0) {
    std::cerr << "Object: header keywords erroneous!" << std::endl;
    std::terminate();
  }

  int naxis;
  fits_get_img_dim(fptr, &naxis, &status);
  if (naxis!=2) {
    std::cerr << "Object: naxis != 2. This is not a FITS image!" << std::endl;
    std::terminate();
  } 
  // first obtain the size of the image array
  long naxes[2] = {1,1};
  fits_get_img_size(fptr, naxis, naxes, &status);
  unsigned int axsize0, axsize1;
  axsize0 = naxes[0];
  axsize1 = naxes[1];
  long npixels = axsize0*axsize1;
  // grid is defined by XMIN..YMAX; this should have the same
  // number of pixels as the actual image
  if (npixels != Image<data_t>::getGrid().size()) {
    std::cerr << "Object: Grid size from header keywords wrong" << std::endl;
    std::terminate();
  }
  
  history << "# Reading object's pixel data" << std::endl;
  Image<data_t>::resize(npixels);
  long firstpix[2] = {1,1};
  data_t vald;
  int imageformat = getFITSImageFormat(vald);
  int datatype = getFITSDataType(vald);
  fits_read_pix(fptr, datatype, firstpix, npixels, NULL,Image<data_t>::c_array(), NULL, &status);
  
  history << ", segmentation map";
  // move to 1st extHDU for the segmentation map
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  segMap.resize(npixels);
  int vali;
  imageformat = getFITSImageFormat(vali);
  datatype = getFITSDataType(vali);
  fits_read_pix(fptr, datatype, firstpix, npixels, NULL,segMap.c_array(), NULL, &status);
  segMap.accessGrid() = Image<data_t>::getGrid();

  // check if there is 2nd extHDU: the weight map
  if (!fits_movabs_hdu(fptr, 3, &hdutype, &status)) {
    history << " and weight map";
    weight.resize(npixels);
    imageformat = getFITSImageFormat(vald);
    datatype = getFITSDataType(vald);
    fits_read_pix(fptr, datatype, firstpix, npixels, NULL,weight.c_array(), NULL, &status);
  }
  history << std::endl;
  fits_close_file(fptr, &status);

  // since the grid is change, the centroid has to be recomputed
  computeFluxCentroid();
}

unsigned int Object::getID() const {
  return id;
}

void Object::setNumber(unsigned int num) {
  number = num;
}

unsigned int Object::getNumber() const {
  return number;
}

const NumVector<data_t>& Object::getWeightMap() const {
  return weight;
}

NumVector<data_t>& Object::accessWeightMap() {
  return weight;
}

const Point2D& Object::getCentroid() const {
  return centroid;
}

void Object::setCentroid(const Point2D& xc) {
  centroid = xc;
  history << "# Centroid set to ("<< centroid(0) << "/" << centroid(1) <<  ") by user." << std::endl;
}

data_t Object::getFlux() const {
  return flux;
}

void Object::setFlux(data_t F) {
  flux = F;
  history << "# Flux set to " << flux << " by user." <<std::endl;
}

NumMatrix<data_t> Object::get2ndBrightnessMoments() {
  const NumVector<data_t>& data = Image<data_t>::getData();
  const Grid& grid = Image<data_t>::getGrid();
  NumMatrix<data_t> Q(2,2);
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
complex<data_t> Object::getEllipticity() {
  NumMatrix<data_t> Q = get2ndBrightnessMoments();
  complex<data_t> I(0,1);
  complex<data_t> Q11(Q(0,0),0),Q22(Q(1,1),0),Q12(Q(0,1),0);
  complex<data_t> epsilon = (Q11 - Q22 + data_t(2)*I*Q12)/(Q11+Q22 + data_t(2)*sqrt(Q11*Q22-Q12*Q12));
  return epsilon;
  // axis ratio
  //data_t e = abs(epsilon);
  //data_t r = (1-e)/(1+e);
}

void Object::computeFluxCentroid() {
  const NumVector<data_t>& data = Image<data_t>::getData();
  const Grid& grid = Image<data_t>::getGrid();
  flux = 0;
  centroid(0) = centroid(1) = 0;
  for (int i=0; i< grid.size(); i++) {
    flux += data(i);
    centroid(0) += grid(i,0) * data(i);
    centroid(1) += grid(i,1) * data(i);
  }
  centroid(0) /= flux;
  centroid(1) /= flux;
  
  history << "# Flux = " << flux << ", Centroid = ("<< centroid(0) << "/" << centroid(1);
  history <<  ")." <<std::endl;
}

unsigned short Object::getDetectionFlag() const {
  return flag;
}

void Object::setDetectionFlag(unsigned short inflag) {
  flag = inflag;
}

data_t Object::getNoiseMean() const {
  return noise_mean;
}

data_t Object::getNoiseRMS() const {
  return noise_rms;
}

void Object::setNoiseMeanRMS(data_t mean, data_t rms) {
  noise_mean = mean;
  noise_rms = rms;
  history << "# Setting noise to (" << mean << ") +- (" << rms << ")" << std::endl;
}

data_t Object::getBlendingProbability() const {
  return blend;
}

void Object::setBlendingProbability(data_t b) {
  blend = b;
}

data_t Object::getStarGalaxyProbability() const {
  return s_g;
}

void Object::setStarGalaxyProbability(data_t sg) {
  s_g = sg;
}

void Object::setBaseFilename(std::string filename) {
  basefilename = filename;
}
std::string Object::getBaseFilename() const {
  return basefilename;
}
const SegmentationMap& Object::getSegmentationMap() const {
  return segMap;
}

SegmentationMap& Object::accessSegmentationMap() {
  return segMap;
}

const PixelCovarianceMatrix& Object::getPixelCovarianceMatrix() const {
  return cov;
}

PixelCovarianceMatrix& Object::accessPixelCovarianceMatrix() {
  return cov;
}

const CorrelationFunction& Object::getCorrelationFunction() const {
  return xi;
}

CorrelationFunction& Object::accessCorrelationFunction() {
  return xi;
}

void Object::save(std::string filename) {
  // write pixel data
  const NumVector<data_t>& data = Image<data_t>::getData();
  const Grid& grid = Image<data_t>::getGrid();
  int status = 0;
  fitsfile *outfptr = createFITSFile(filename);
  status = writeFITSImage(outfptr,grid,data);

  // add object information to header
  status = updateFITSKeywordString(outfptr,"BASEFILE",basefilename,"name of source file");
  status = updateFITSKeyword(outfptr,"ID",id,"object id");
  status = updateFITSKeyword(outfptr,"NUMBER",number,"object number");
  status = updateFITSKeyword(outfptr,"XMIN",grid.getStartPosition(0),"min(X) in image pixels");
  status = updateFITSKeyword(outfptr,"XMAX",grid.getStopPosition(0),"max(X) in image pixels");
  status = updateFITSKeyword(outfptr,"YMIN",grid.getStartPosition(1),"min(Y) in image pixels");
  status = updateFITSKeyword(outfptr,"YMAX",grid.getStopPosition(1),"min(Y) in image pixels");
  complex<data_t> xc(centroid(0),centroid(1));
  status = updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");
  status = updateFITSKeyword(outfptr,"BG_MEAN",noise_mean,"mean of background noise");
  status = updateFITSKeyword(outfptr,"BG_RMS",noise_rms,"rms of background noise");
  status = updateFITSKeyword(outfptr,"FLAG",flag,"extraction flag");
  status = updateFITSKeyword(outfptr,"BLEND",blend,"blending probability");
  status = updateFITSKeyword(outfptr,"S_G",s_g,"stellarity");
  status = appendFITSHistory(outfptr,history.str());

  // save segMap 
  status = addFITSExtension(outfptr,"SEGMAP",grid,segMap.getData());
  status = appendFITSHistory(outfptr,(segMap.getHistory()).str());

  //if weight map provided, save it too
  if (weight.size() != 0)
    status = addFITSExtension(outfptr,"WEIGHT",grid,weight);

  fits_close_file(outfptr, &status);
}
