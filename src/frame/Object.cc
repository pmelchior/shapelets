#include <frame/Object.h>
#include <sstream>
#include <gsl/gsl_math.h>
#include <fitsio.h>

Object::Object(unsigned long inid) : Image<data_t>(), segMap() {
  id = inid;
  flag = 0;
  classifier = 0;
  flux = centroid(0) = centroid(1) = 0;
}

Object::Object(std::string objfile) : Image<data_t>(), segMap() {
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
  status = readFITSKeyword(fptr,"FLUX",flux);
  status = readFITSKeyword(fptr,"BG_MEAN",noise_mean);
  status = readFITSKeyword(fptr,"BG_RMS",noise_rms);
  unsigned long flags;
  status = readFITSKeyword(fptr,"FLAG",flags);
  flag = std::bitset<8>(flags);
  status = readFITSKeyword(fptr,"CLASSIFIER",classifier);
  
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
  
  history << "# Reading object's pixel data";
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
  segMap.resize(npixels);
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
  status = updateFITSKeyword(outfptr,"XMIN",grid.getStartPosition(0),"min(X) in image pixels");
  status = updateFITSKeyword(outfptr,"XMAX",grid.getStopPosition(0),"max(X) in image pixels");
  status = updateFITSKeyword(outfptr,"YMIN",grid.getStartPosition(1),"min(Y) in image pixels");
  status = updateFITSKeyword(outfptr,"YMAX",grid.getStopPosition(1),"min(Y) in image pixels");
  status = updateFITSKeyword(outfptr,"FLUX",flux,"flux in ADUs");
  complex<data_t> xc(centroid(0),centroid(1));
  status = updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");
  status = updateFITSKeyword(outfptr,"BG_MEAN",noise_mean,"mean of background noise");
  status = updateFITSKeyword(outfptr,"BG_RMS",noise_rms,"rms of background noise");
  status = updateFITSKeyword(outfptr,"FLAG",flag.to_ulong(),"extraction flag");
  status = updateFITSKeyword(outfptr,"CLASSIFIER",classifier,"object classifier");
  status = appendFITSHistory(outfptr,history.str());

  // save segMap
  if (segMap.size() != 0) {
    status = writeFITSImage(outfptr,grid,segMap.getData(),"SEGMAP");
    status = appendFITSHistory(outfptr,(segMap.getHistory()).str());
  }

  //if weight map provided, save it too
  if (weight.size() != 0)
    status = writeFITSImage(outfptr,grid,weight,"WEIGHT");

  status = closeFITSFile(outfptr);
}

void Object::setID(unsigned long inid) {
  id = inid;
}

unsigned long Object::getID() const {
  return id;
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
  
  // check if weights are available: if yes, use them
  bool weights = false;
  if (weight.size() != 0)
    weights = true;
  
  data_t datapoint, unnormed_flux = 0;
  for (int i=0; i< grid.size(); i++) {
    if (weights) {
      if (weight(i) != 0)
	datapoint = data(i) * sqrt(weight(i));
      else
	datapoint = 0;
    }
    else
      datapoint = data(i);
    
    unnormed_flux += datapoint;
    Q(0,0) += gsl_pow_2(grid(i,0)-centroid(0)) * datapoint;
    Q(0,1) += (grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * datapoint;
    Q(1,1) += gsl_pow_2(grid(i,1)-centroid(1)) * datapoint;
  }
  
  // since unnormed_flux and Qs have same normalization, so it drops out
  Q(0,0) /= unnormed_flux;
  Q(0,1) /= unnormed_flux;
  Q(1,0) = Q(0,1);
  Q(1,1) /= unnormed_flux;
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

  // check if weights are available: if yes, use them
  bool weights = false;
  if (weight.size() != 0)
    weights = true;
  
  flux = 0;
  centroid(0) = centroid(1) = 0;
  data_t datapoint, sum_weights = 0;
  unsigned long sum_pixels = 0;
  for (int i=0; i< grid.size(); i++) {
    if (weights) {
      if (weight(i) != 0) {
	datapoint = data(i) * sqrt(weight(i));
	sum_weights += sqrt(weight(i));
	sum_pixels++;
      } else
	datapoint = 0;
    }
    else
      datapoint = data(i);
    flux += datapoint;
    centroid(0) += grid(i,0) * datapoint;
    centroid(1) += grid(i,1) * datapoint;
  }
  
  // even for weighted centroid: sum_weights drops out
  centroid(0) /= flux;
  centroid(1) /= flux;
  if (weights)
    flux /= sum_weights/sum_pixels;
  
  history << "# Flux = " << flux << ", Centroid = ("<< centroid(0) << "/" << centroid(1) << ")" << std::endl;
}

std::bitset<8> Object::getDetectionFlags() const {
  return flag;
}

void Object::setDetectionFlags(const std::bitset<8>& inflag) {
  flag = std::bitset<8>(inflag.to_ulong());
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

data_t Object::getClassifier() const {
  return classifier;
}

void Object::setClassifier(data_t c) {
  classifier = c;
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

std::string Object::getHistory() const {
  return history.str();
}

