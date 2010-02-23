#include "../../include/frame/Object.h"
#include <sstream>
#include <gsl/gsl_math.h>

using namespace shapelens;

Object::Object() : Image<data_t>(), segMap() {
  id = 0;
  flags = 0;
  classifier = 0;
  flux = centroid(0) = centroid(1) = 0;
}

Object::Object (const Image<data_t>& base) : Image<data_t>(base), segMap()  {
  id = 0;
  flags = 0;
  classifier = 0;
  flux = centroid(0) = centroid(1) = 0;
}

void Object::operator=(const Image<data_t>& base) {
  Image<data_t>::operator=(base);
}

Object::Object(std::string objfile) : Image<data_t>(), segMap() {
  int status, nkeys, keypos, hdutype;
  char card[FLEN_CARD];
  char comment[FLEN_CARD];

  history << "# Loading object from Fits file " << objfile << ":" << std::endl;
  fitsfile* fptr = IO::openFITSFile(objfile);

  // reading objects pixel data
  history << "# Reading object's pixel data";
  IO::readFITSImage(fptr,*this);
  
  // recover object information from header keywords
  IO::readFITSKeywordString(fptr,"BASEFILE",Image<data_t>::basefilename);
  IO::readFITSKeyword(fptr,"ID",id);
  int xmin,ymin;
  IO::readFITSKeyword(fptr,"XMIN",xmin);
  IO::readFITSKeyword(fptr,"YMIN",ymin);
  Image<data_t>::grid = Grid(xmin,ymin,grid.getSize(0),grid.getSize(1));
  complex<data_t> xc;
  IO::readFITSKeyword(fptr,"CENTROID",xc);
  centroid(0) = real(xc);
  centroid(1) = imag(xc);
  IO::readFITSKeyword(fptr,"FLUX",flux);
  IO::readFITSKeyword(fptr,"BG_MEAN",noise_mean);
  IO::readFITSKeyword(fptr,"BG_RMS",noise_rms);
  unsigned long f;
  IO::readFITSKeyword(fptr,"FLAG",f);
  flags = std::bitset<8>(f);
  IO::readFITSKeyword(fptr,"CLASSIFIER",classifier);
  
  // read history
  std::string hstr;
  IO::readFITSKeyCards(fptr,"HISTORY",hstr);

  // check whether grid has same size as object
  if (Object::size() != Image<data_t>::grid.size())
    throw std::invalid_argument("Object: Grid size from header keywords wrong!");
  
  history << ", segmentation map";
  // move to 1st extHDU for the segmentation map
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  IO::readFITSImage(fptr, segMap);

  // check if there is 2nd extHDU: the weight map or correlation
  if (!fits_movabs_hdu(fptr, 3, &hdutype, &status)) {
    std::string extname;
    IO::readFITSKeywordString(fptr, "EXTNAME", extname);
    if (extname == "WEIGHT") {
      history << " and weight map";
      IO::readFITSImage(fptr, weight);
    } else if (extname == "CORRELATION") {
      history << " and correlation function";
      NumMatrix<data_t> corr;
      IO::readFITSImage(fptr, corr);
      xi = CorrelationFunction(corr);
    }
  }
  history << std::endl;

  // append pHDUs history
  history.setSilent();
  history << hstr;
  history.unsetSilent();

  IO::closeFITSFile(fptr);
}

void Object::save(std::string filename) {
  // write pixel data
  fitsfile *outfptr = IO::createFITSFile(filename);
  IO::writeFITSImage(outfptr,*this);

  // add object information to header
  IO::updateFITSKeywordString(outfptr,"BASEFILE",Image<data_t>::basefilename,"name of source file");
  IO::updateFITSKeyword(outfptr,"ID",id,"object id");
  IO::updateFITSKeyword(outfptr,"XMIN",Image<data_t>::grid.getStartPosition(0),"min(X) in image pixels");
  IO::updateFITSKeyword(outfptr,"YMIN",Image<data_t>::grid.getStartPosition(1),"min(Y) in image pixels");
  IO::updateFITSKeyword(outfptr,"FLUX",flux,"flux in ADUs");
  complex<data_t> xc(centroid(0),centroid(1));
  IO::updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");
  IO::updateFITSKeyword(outfptr,"BG_MEAN",noise_mean,"mean of background noise");
  IO::updateFITSKeyword(outfptr,"BG_RMS",noise_rms,"rms of background noise");
  IO::updateFITSKeyword(outfptr,"FLAG",flags.to_ulong(),"extraction flags");
  IO::updateFITSKeyword(outfptr,"CLASSIFIER",classifier,"object classifier");
  IO::appendFITSHistory(outfptr,Image<data_t>::history.str());

  // save segMap
  if (segMap.size() != 0) {
    IO::writeFITSImage(outfptr,segMap,"SEGMAP");
    IO::appendFITSHistory(outfptr,segMap.history.str());
  }

  //if weight map provided, save it too
  if (weight.size() != 0)
    IO::writeFITSImage(outfptr,weight,"WEIGHT");
  //if correlationFunction is provided, save it
  if (xi.getCorrelationFunction().size() > 0)
    IO::writeFITSImage(outfptr,xi.getCorrelationMatrix(),"CORRELATION");

  IO::closeFITSFile(outfptr);
}

void Object::computeFlux() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  flux = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	flux += data(i) * weight(i);
	sum_weights += weight(i);
      }
    }
    flux /= sum_weights;
  }
  else // unweigthed
    for (long i=0; i< grid.size(); i++)
      flux += data(i);
}

void Object::computeCentroid() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  centroid(0) = centroid(1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	centroid(0) += data(i) * grid(i,0) * weight(i);
	centroid(1) += data(i) * grid(i,1) * weight(i);
	sum_weights += weight(i);
      }
    }
    centroid(0) /= flux * sum_weights;
    centroid(1) /= flux * sum_weights;
  }
  else { // unweighted
    for (long i=0; i< grid.size(); i++) {
      centroid(0) += grid(i,0) * data(i);
      centroid(1) += grid(i,1) * data(i);
    }
    centroid(0) /= flux;
    centroid(1) /= flux;
  }
}

void Object::computeFluxCentroid() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  flux = 0;
  centroid(0) = centroid(1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (long i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	flux += data(i) * weight(i);
	centroid(0) += data(i) * grid(i,0) * weight(i);
	centroid(1) += data(i) * grid(i,1) * weight(i);
	sum_weights += weight(i);
      }
    }
    flux /= sum_weights;
    centroid(0) /= flux * sum_weights;
    centroid(1) /= flux * sum_weights;
  }
  else { // unweigthed
    for (long i=0; i< grid.size(); i++) {
      flux += data(i);
      centroid(0) += grid(i,0) * data(i);
      centroid(1) += grid(i,1) * data(i);
    }
    centroid(0) /= flux;
    centroid(1) /= flux;
  }
}

void Object::computeCorrelationFunction(data_t threshold) {
  if (segMap.size()) // if a segMap is provided, mask object pixels
    xi = CorrelationFunction(*this,segMap,threshold);
  else
    xi = CorrelationFunction(*this,threshold);
}

void Object::computeFFT() {
  FFT::transform(*this,fourier);
}

void Object::convolve(const Object& kernel) {

  // compute fourier if not done yet
  if(fourier.getRealSize(0) == 0)
    computeFFT();

  // kernel.fourier does not exist but size is ok
  if (kernel.fourier.getRealSize(0) == 0 && 
      kernel.getSize(0) == getSize(0) && kernel.getSize(1) == getSize(1)) {
    int M = getSize(0);
    int N = getSize(1);
    FourierTransform2D kernel_transformed(M,N);
    FFT::transform(kernel,kernel_transformed);
    FFT::conv_multiply(fourier,kernel_transformed,fourier);
    FFT::transform(fourier,*this);
  }

  // size does not fit
  else if (kernel.getSize(0) != getSize(0) || kernel.getSize(1) != getSize(1))
    FFT::convolve(*this,fourier,kernel);

  // everything is already set up
  else {
    FFT::conv_multiply(fourier,kernel.fourier,fourier);
    FFT::transform(fourier,*this);
  }
}
