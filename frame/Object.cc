#include <Object.h>
#include <sstream>
#include <gsl/gsl_math.h>
#include <fitsio.h>

Object::Object(unsigned int inid) : Image<double>(), segMap(*this) {
  id = inid;
  flag = 0;
  blend = s_g = -1;
  flux = centroid(0) = centroid(1) = -1;
}

Object::Object(std::string objfile) : Image<double>(), segMap(*this) {
  fitsfile *fptr;
  int status, nkeys, keypos, hdutype;
  char card[FLEN_CARD];
  char comment[FLEN_CARD];
  status = 0; 

  std::ostringstream text;
  history.append("# Loading object from Fits file "+objfile+":\n");
  fits_open_file(&fptr, objfile.c_str(), READONLY, &status);

  // open pHDU and read header
  fits_movabs_hdu(fptr, 1, &hdutype, &status);
  // recover object information from header keywords
  fits_read_key (fptr,TSTRING,"BASEFILE",card,comment, &status);
  basefilename = std::string(card);
  fits_read_key (fptr,TINT,"ID",&id,comment, &status);
  double grid_min_0, grid_max_0, grid_min_1, grid_max_1;
  fits_read_key (fptr,TDOUBLE,"XMIN",&grid_min_0,comment, &status);
  fits_read_key (fptr,TDOUBLE,"XMAX",&grid_max_0,comment, &status);
  fits_read_key (fptr,TDOUBLE,"YMIN",&grid_min_1,comment, &status);
  fits_read_key (fptr,TDOUBLE,"YMAX",&grid_max_1,comment, &status);
  Image<double>::accessGrid() = Grid(grid_min_0,grid_max_0,1,grid_min_1,grid_max_1,1);
  fits_read_key (fptr,TDOUBLE,"BG_MEAN",&noise_mean,comment, &status);
  fits_read_key (fptr,TDOUBLE,"BG_RMS",&noise_rms,comment, &status);
  fits_read_key (fptr,TSTRING,"NOISE",card,comment, &status);
  noisemodel = std::string(card);
  fits_read_key (fptr,TUSHORT,"FLAG",&flag,comment, &status);
  fits_read_key (fptr,TDOUBLE,"BLEND",&blend,comment, &status);
  fits_read_key (fptr,TDOUBLE,"S_G",&s_g,comment, &status);

  // if status is not 0 we have faced an error
  if (status != 0) {
    std::cout << "Object: header keywords erroneous!" << std::endl;
    std::terminate();
  }

  // write header information to history
  text << "# Object " << id << " from " + basefilename + "\n";
  text << "# Segment area =  (" << grid_min_0 << "/" << grid_min_1 << ") .. (";
  text << grid_max_0 << "/" << grid_max_1 << ")" << std::endl;
  text << "# Noise model = " + noisemodel << std::endl;
  text << "# Background noise = (" << noise_mean << ") +- (" << noise_rms << ")" << std::endl;
  text << "# Segmentation flag = "<< flag << ", Blending probability = ";
  text << blend << ", Stellarity = " << s_g << std::endl;
  history.append(text);

  int naxis;
  fits_get_img_dim(fptr, &naxis, &status);
  if (naxis!=2) {
    std::cout << "Object: naxis != 2. This is not a FITS image!" << std::endl;
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
  if (npixels != Image<double>::getGrid().size()) {
    std::cout << "Object: Grid size from header keywords wrong" << std::endl;
    std::terminate();
  }
  
  history.append("# Reading object's pixel data");
  Image<double>::resize(npixels);
  long firstpix[2] = {1,1};
  double vald;
  int imageformat, datatype;
  setFITSTypes(vald,imageformat,datatype);
  fits_read_pix(fptr, datatype, firstpix, npixels, NULL,Image<double>::c_array(), NULL, &status);
  
  history.append(", segmentation map");
  // move to 1st extHDU for the segmentation map
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  segMap.resize(npixels);
  int vali;
  setFITSTypes(vali,imageformat,datatype);
  fits_read_pix(fptr, datatype, firstpix, npixels, NULL,segMap.c_array(), NULL, &status);
  segMap.accessGrid() = Image<double>::getGrid();

  // check if there is 2nd extHDU: the weight map
  if (!fits_movabs_hdu(fptr, 3, &hdutype, &status)) {
    history.append(" and weight map");
    weight.resize(npixels);
    setFITSTypes(vald,imageformat,datatype);
    fits_read_pix(fptr, datatype, firstpix, npixels, NULL,weight.c_array(), NULL, &status);
  }
  history.append("\n");
  fits_close_file(fptr, &status);

  // since the grid is change, the centroid has to be recomputed
  computeFluxCentroid();
}

unsigned int Object::getID() const {
  return id;
}

const NumVector<double>& Object::getWeightMap() const {
  return weight;
}

NumVector<double>& Object::accessWeightMap() {
  return weight;
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
  const NumVector<double>& data = Image<double>::getData();
  const Grid& grid = Image<double>::getGrid();
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
  const NumVector<double>& data = Image<double>::getData();
  const Grid& grid = Image<double>::getGrid();
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
  std::ostringstream text;
  text << "# Setting noise to (" << mean << ") +- (" << rms << ")" << std::endl;
  history.append(text);
}

void Object::setNoiseModel(std::string innoisemodel) {
  noisemodel = innoisemodel;
  history.append("# Setting noise model to "+noisemodel+"\n");
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
  const NumVector<double>& data = Image<double>::getData();
  const Grid& grid = Image<double>::getGrid();
  writeFITSFile(filename,grid,data);

  // add object information to header
  int status = 0;
  fitsfile *outfptr;
  // open pHUD of fits file
  fits_open_data(&outfptr,filename.c_str(),READWRITE,&status);
  fits_write_key (outfptr,TSTRING,"BASEFILE",const_cast<char*>(basefilename.c_str()),"name of source file", &status);
  fits_write_key (outfptr,TINT,"ID",&id,"object id", &status);
  double buffer;
  buffer = grid.getStartPosition(0);
  fits_write_key (outfptr,TDOUBLE,"XMIN",&buffer,"min(X) in image pixels", &status);
  buffer = grid.getStopPosition(0);
  fits_write_key (outfptr,TDOUBLE,"XMAX",&buffer,"max(X) in image pixels", &status);
  buffer = grid.getStartPosition(1);
  fits_write_key (outfptr,TDOUBLE,"YMIN",&buffer,"min(Y) in image pixels", &status);
  buffer = grid.getStopPosition(1);
  fits_write_key (outfptr,TDOUBLE,"YMAX",&buffer,"max(Y) in image pixels", &status);
  fits_write_key (outfptr,TDOUBLE,"BG_MEAN",&noise_mean,"mean of background noise", &status);
  fits_write_key (outfptr,TDOUBLE,"BG_RMS",&noise_rms,"rms of background noise", &status);
  fits_write_key (outfptr,TSTRING,"NOISE",const_cast<char*>(noisemodel.c_str()),"noise model", &status);
  fits_write_key (outfptr,TUSHORT,"FLAG",&flag,"extraction flag", &status);
  fits_write_key (outfptr,TDOUBLE,"BLEND",&blend,"blending probability", &status);
  fits_write_key (outfptr,TDOUBLE,"S_G",&s_g,"stellarity", &status);
  fits_close_file(outfptr, &status);

  // save segMap 
  //keywords.clear();
  addFITSExtension(filename,"SEGMAP",grid,segMap.getData());
  //if weight map provided, save it too
  if (weight.size() != 0)
    addFITSExtension(filename,"WEIGHT",grid,weight);
}
