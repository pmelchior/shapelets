#include <frame/Object.h>
#include <sstream>
#include <gsl/gsl_math.h>
#include <IO.h>

Object::Object(unsigned long inid) : Image<data_t>(), segMap() {
  id = inid;
  flags = 0;
  classifier = 0;
  flux = centroid(0) = centroid(1) = 0;
}

Object::Object (const Image<data_t>& base) : Image<data_t>(base), segMap() {
}

Object::Object(std::string objfile) : Image<data_t>(), segMap() {
  int status, nkeys, keypos, hdutype;
  char card[FLEN_CARD];
  char comment[FLEN_CARD];
  status = 0;

  history << "# Loading object from Fits file " << objfile << ":" << std::endl;
  fitsfile* fptr = IO::openFITSFile(objfile);

  // reading objects pixel data
  history << "# Reading object's pixel data";
  IO::readFITSImage(fptr,grid,Object::accessNumVector());
  
  // recover object information from header keywords
  status = IO::readFITSKeywordString(fptr,"BASEFILE",basefilename);
  status = IO::readFITSKeyword(fptr,"ID",id);
  grid_t xmin,ymin;
  status = IO::readFITSKeyword(fptr,"XMIN",xmin);
  status = IO::readFITSKeyword(fptr,"YMIN",ymin);
  Image<data_t>::grid = Grid(xmin,ymin,grid.getSize(0),grid.getSize(1));
  complex<data_t> xc;
  status = IO::readFITSKeyword(fptr,"CENTROID",xc);
  centroid(0) = real(xc);
  centroid(1) = imag(xc);
  status = IO::readFITSKeyword(fptr,"FLUX",flux);
  status = IO::readFITSKeyword(fptr,"BG_MEAN",noise_mean);
  status = IO::readFITSKeyword(fptr,"BG_RMS",noise_rms);
  unsigned long f;
  status = IO::readFITSKeyword(fptr,"FLAG",f);
  flags = std::bitset<8>(f);
  status = IO::readFITSKeyword(fptr,"CLASSIFIER",classifier);
  
  // read history
  std::string hstr;
  IO::readFITSKeyCards(fptr,"HISTORY",hstr);

  // check whether grid has same size as object
  if (Object::size() != Image<data_t>::grid.size()) {
    std::cerr << "Object: Grid size from header keywords wrong" << std::endl;
    std::terminate();
  }
  
  history << ", segmentation map";
  // move to 1st extHDU for the segmentation map
  fits_movabs_hdu(fptr, 2, &hdutype, &status);
  IO::readFITSImage(fptr, segMap.grid, segMap);

  // check if there is 2nd extHDU: the weight map or correlation
  if (!fits_movabs_hdu(fptr, 3, &hdutype, &status)) {
    std::string extname;
    IO::readFITSKeywordString(fptr, "EXTNAME", extname);
    if (extname == "WEIGTH") {
      history << " and weight map";
      Grid weightgrid;
      IO::readFITSImage(fptr, weightgrid, weight);
    } else if (extname == "CORRELATION") {
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
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  int status = 0;
  fitsfile *outfptr = IO::createFITSFile(filename);
  status = IO::writeFITSImage(outfptr,grid,data);

  // add object information to header
  status = IO::updateFITSKeywordString(outfptr,"BASEFILE",basefilename,"name of source file");
  status = IO::updateFITSKeyword(outfptr,"ID",id,"object id");
  status = IO::updateFITSKeyword(outfptr,"XMIN",grid.getStartPosition(0),"min(X) in image pixels");
  status = IO::updateFITSKeyword(outfptr,"YMIN",grid.getStartPosition(1),"min(Y) in image pixels");
  status = IO::updateFITSKeyword(outfptr,"FLUX",flux,"flux in ADUs");
  complex<data_t> xc(centroid(0),centroid(1));
  status = IO::updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");
  status = IO::updateFITSKeyword(outfptr,"BG_MEAN",noise_mean,"mean of background noise");
  status = IO::updateFITSKeyword(outfptr,"BG_RMS",noise_rms,"rms of background noise");
  status = IO::updateFITSKeyword(outfptr,"FLAG",flags.to_ulong(),"extraction flags");
  status = IO::updateFITSKeyword(outfptr,"CLASSIFIER",classifier,"object classifier");
  status = IO::appendFITSHistory(outfptr,history.str());

  // save segMap
  if (segMap.size() != 0) {
    status = IO::writeFITSImage(outfptr,grid,segMap,"SEGMAP");
    status = IO::appendFITSHistory(outfptr,segMap.history.str());
  }

  //if weight map provided, save it too
  if (weight.size() != 0)
    status = IO::writeFITSImage(outfptr,grid,weight,"WEIGHT");
  //if correlationFunction is provided, save it
  if (xi.getCorrelationFunction().size() > 0)
    status = IO::writeFITSImage(outfptr,xi.getCorrelationMatrix(),"CORRELATION");

  status = IO::closeFITSFile(outfptr);
}

NumMatrix<data_t> Object::get2ndBrightnessMoments() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
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
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;

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
