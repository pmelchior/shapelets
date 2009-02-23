#include <frame/Object.h>
#include <sstream>
#include <gsl/gsl_math.h>

Object::Object() : Image<data_t>(), segMap() {
  id = 0;
  flags = 0;
  classifier = 0;
  flux = centroid(0) = centroid(1) = 0;
}

Object::Object (const Image<data_t>& base) : Image<data_t>(base), segMap() {
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
  status = 0;

  history << "# Loading object from Fits file " << objfile << ":" << std::endl;
  fitsfile* fptr = IO::openFITSFile(objfile);

  // reading objects pixel data
  history << "# Reading object's pixel data";
  IO::readFITSImage(fptr,grid,Object::accessNumVector());
  
  // recover object information from header keywords
  status = IO::readFITSKeywordString(fptr,"BASEFILE",Image<data_t>::basefilename);
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
  status = IO::updateFITSKeywordString(outfptr,"BASEFILE",Image<data_t>::basefilename,"name of source file");
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

void Object::computeFlux() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  flux = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (int i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	flux += data(i) * weight(i);
	sum_weights += weight(i);
      }
    }
    flux /= sum_weights;
  }
  else // unweigthed
    for (int i=0; i< grid.size(); i++)
      flux += data(i);
}

void Object::computeCentroid() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  centroid(0) = centroid(1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (int i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	centroid(0) += data(i) * grid(i,0) * weight(i);
	centroid(0) += data(i) * grid(i,1) * weight(i);
	sum_weights += weight(i);
      }
    }
    centroid(0) /= flux * sum_weights;
    centroid(1) /= flux * sum_weights;
  }
  else { // unweighted
    for (int i=0; i< grid.size(); i++) {
      centroid(0) += grid(i,0) * data(i);
      centroid(1) += grid(i,1) * data(i);
    }
    centroid(0) /= flux;
    centroid(1) /= flux;
  }
}

void Object::computeQuadrupole() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  Q(0,0) = Q(0,1) = Q(1,1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (int i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	Q(0,0) += gsl_pow_2(grid(i,0)-centroid(0)) * data(i) * weight(i);
	Q(0,1) += (grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i) * weight(i);
	Q(1,1) += gsl_pow_2(grid(i,1)-centroid(1)) * data(i) * weight(i);
	sum_weights += weight(i);
      }
    }
    Q(0,0) /= flux * sum_weights;
    Q(0,1) /= flux * sum_weights;
    Q(1,1) /= flux * sum_weights;
  }
  else { // unweighted
    for (int i=0; i< grid.size(); i++) {
      Q(0,0) += gsl_pow_2(grid(i,0)-centroid(0)) * data(i);
      Q(0,1) += (grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i);
      Q(1,1) += gsl_pow_2(grid(i,1)-centroid(1)) * data(i);
    }
    Q(0,0) /= flux;
    Q(0,1) /= flux;
    Q(1,1) /= flux;
  }
}

void Object::computeOctupole() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  O(0,0,0) = O(0,0,1) = O(0,1,1) = O(1,1,1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (int i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	O(0,0,0) += gsl_pow_3(grid(i,0)-centroid(0)) * data(i) * weight(i);
	O(0,0,1) += gsl_pow_2(grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i) * weight(i);
	O(0,1,1) += (grid(i,0)-centroid(0))*gsl_pow_2(grid(i,1)-centroid(1)) * data(i) * weight(i);
	O(1,1,1) += gsl_pow_3(grid(i,1)-centroid(1)) * data(i) * weight(i);
	sum_weights += weight(i);
      }
    }
    O(0,0,0) /= flux * sum_weights;
    O(0,0,1) /= flux * sum_weights;
    O(0,1,1) /= flux * sum_weights;
    O(1,1,1) /= flux * sum_weights;
  }
  else { // unweighted
    for (int i=0; i< grid.size(); i++) {
      O(0,0,0) += gsl_pow_3(grid(i,0)-centroid(0)) * data(i);
      O(0,0,1) += gsl_pow_2(grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i);
      O(0,1,1) += (grid(i,0)-centroid(0))*gsl_pow_2(grid(i,1)-centroid(1)) * data(i);
      O(1,1,1) += gsl_pow_3(grid(i,1)-centroid(1)) * data(i);
    }
    O(0,0,0) /= flux;
    O(0,0,1) /= flux;
    O(0,1,1) /= flux;
    O(1,1,1) /= flux;
  }
}

void Object::computeHexadecupole() {
  const NumVector<data_t>& data = *this;
  const Grid& grid = Image<data_t>::grid;
  H(0,0,0,0) = H(0,0,0,1) = H(0,0,1,1) = H(0,1,1,1) = H(1,1,1,1) = 0;
  // check if weights are available: if yes, use them
  if (weight.size() != 0) {
    data_t sum_weights = 0;
    for (int i=0; i< grid.size(); i++) {
      if (weight(i) > 0) {
	H(0,0,0,0) += gsl_pow_4(grid(i,0)-centroid(0)) * data(i) * weight(i);
	H(0,0,0,1) += gsl_pow_3(grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i) * weight(i);
	H(0,0,1,1) += gsl_pow_2(grid(i,0)-centroid(0))*gsl_pow_2(grid(i,1)-centroid(1)) * data(i) * weight(i);
	H(0,1,1,1) += (grid(i,0)-centroid(0))*gsl_pow_3(grid(i,1)-centroid(1)) * data(i) * weight(i);
	H(1,1,1,1) = gsl_pow_4(grid(i,1)-centroid(1)) * data(i) * weight(i);
	sum_weights += weight(i);
      }
    }
    H(0,0,0,0) /= flux * sum_weights;
    H(0,0,0,1) /= flux * sum_weights;
    H(0,0,1,1) /= flux * sum_weights;
    H(0,1,1,1) /= flux * sum_weights;
    H(1,1,1,1) /= flux * sum_weights;
  }
  else { // unweighted
    for (int i=0; i< grid.size(); i++) {
      H(0,0,0,0) += gsl_pow_4(grid(i,0)-centroid(0)) * data(i);
      H(0,0,0,1) += gsl_pow_3(grid(i,0)-centroid(0))*(grid(i,1)-centroid(1)) * data(i);
      H(0,0,1,1) += gsl_pow_2(grid(i,0)-centroid(0))*gsl_pow_2(grid(i,1)-centroid(1)) * data(i);
      H(0,1,1,1) += (grid(i,0)-centroid(0))*gsl_pow_3(grid(i,1)-centroid(1)) * data(i);
      H(1,1,1,1) = gsl_pow_4(grid(i,1)-centroid(1)) * data(i);
    }
    H(0,0,0,0) /= flux;
    H(0,0,0,1) /= flux;
    H(0,0,1,1) /= flux;
    H(0,1,1,1) /= flux;
    H(1,1,1,1) /= flux;
  }
}

void Object::computeMoments() {
  computeQuadrupole();
  computeOctupole();
  computeHexadecupole();
}

void Object::computeCorrelationFunction(data_t threshold) {
  if (segMap.size()) // if a segMap is provided, mask object pixels
    xi = CorrelationFunction(*this,segMap,threshold);
  else
    xi = CorrelationFunction(*this,2,threshold);
}
