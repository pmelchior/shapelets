#include <shapelets/SIFFile.h>
#include <utils/IO.h>
#include <ShapeLensConfig.h>
#include <bitset>
#include <fstream>

using namespace std;

SIFFile::SIFFile(std::string infilename) {
  filename = infilename;
}

void SIFFile::save(const ShapeletObject& sobj) {
  fitsfile *outfptr = IO::createFITSFile(filename);
  saveSObj(outfptr, sobj);
  IO::closeFITSFile(outfptr);
}

void SIFFile::save(const ShapeletObject& sobj, const ShapeletObject& unreg) {
  fitsfile *outfptr = IO::createFITSFile(filename);
  saveSObj(outfptr, sobj);
  saveSObj(outfptr, unreg);
  IO::closeFITSFile(outfptr);
}

void SIFFile::save(ShapeletObjectList& sl) {
  fitsfile *outfptr = IO::createFITSFile(filename);
  for(ShapeletObjectList::iterator iter = sl.begin(); iter != sl.end(); iter++) 
    saveSObj(outfptr, *(*iter));
  IO::closeFITSFile(outfptr);
}

void SIFFile::saveSObj(fitsfile* outfptr, const ShapeletObject& sobj) {
  int status = IO::writeFITSImage(outfptr,sobj.getCoeffs().getCoeffMatrix());

  // when present, save errors also
  bool saveErrors = 0;
  const NumMatrix<data_t>& cov = sobj.getCovarianceMatrix();
  if (cov.getRows() != 0 && cov.getColumns() != 0)
    saveErrors = 1;
  

  // add shapelet parameters and other necessary information in pHDU header
  IO::updateFITSKeyword(outfptr,"VERSION",(unsigned int) 1,"SIF version");

  // ** Shapelet parameters **
  fits_write_record(outfptr,"",&status);
  fits_write_record(outfptr,"        / Shapelet parameters  /",&status);
  IO::updateFITSKeyword(outfptr,"BETA",sobj.getBeta(),"scale size in pixel units");
  IO::updateFITSKeyword(outfptr,"DIM",sobj.coeffs.getNMax()+1,"dimensions in shapelet space (nmax+1)");
  IO::updateFITSKeyword(outfptr,"CHI2",sobj.getChiSquare(),"decomposition quality");
  IO::updateFITSKeyword(outfptr,"ERRORS",saveErrors,"whether coefficient errors are available");
  IO::updateFITSKeyword(outfptr,"FLAGS",sobj.getFlags().to_ulong(),"extraction and decomposition flags");
  IO::updateFITSKeywordString(outfptr, "EXTNAME", sobj.getName(), "shapelet object name");
  IO::updateFITSKeyword(outfptr, "TAG", sobj.getTag(), "shapelet object tag");

  // ** Frame parameters **
  fits_write_record(outfptr,"        / Frame parameters     /",&status);
  IO::updateFITSKeywordString(outfptr,"BASEFILE",sobj.getBaseFilename(),"originating data file");
  IO::updateFITSKeyword(outfptr,"ID",sobj.getObjectID(),"object id in BASEFILE");
  IO::updateFITSKeyword(outfptr,"CLASSIFIER",sobj.getObjectClassifier(),"object classifier");
  const Grid& grid = sobj.getGrid();
  IO::updateFITSKeyword(outfptr,"XMIN",grid.getStartPosition(0),"min(X) in image pixels");
  IO::updateFITSKeyword(outfptr,"XMAX",grid.getStopPosition(0),"max(X) in image pixels");
  IO::updateFITSKeyword(outfptr,"YMIN",grid.getStartPosition(1),"min(Y) in image pixels");
  IO::updateFITSKeyword(outfptr,"YMAX",grid.getStopPosition(1),"min(Y) in image pixels");
  const Point2D<data_t>& centroid = sobj.getCentroid();
  complex<data_t> xc(centroid(0),centroid(1));
  IO::updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");

  // ** ShapeLensConfig parameters **
  fits_write_record(outfptr,"        / ShapeLensConfig parameters /",&status);
  IO::updateFITSKeywordString(outfptr,"NOISEMODEL",ShapeLensConfig::NOISEMODEL,"noise model");
  IO::updateFITSKeyword(outfptr," NMAX_LOW",ShapeLensConfig::NMAX_LOW,"lower bound for n_max");
  IO::updateFITSKeyword(outfptr,"NMAX_HIGH",ShapeLensConfig::NMAX_HIGH,"upper bound for n_max");
  IO::updateFITSKeyword(outfptr," BETA_LOW",ShapeLensConfig::BETA_LOW,"lower bound for beta");
  IO::updateFITSKeyword(outfptr,"BETA_HIGH",ShapeLensConfig::BETA_HIGH,"upper bound for beta");
  IO::updateFITSKeyword(outfptr,"ALLOW_FLATTENING",ShapeLensConfig::ALLOW_FLATTENING,"whether flattening of chi^2 is allowed");
  IO::updateFITSKeyword(outfptr,"FILTER_SPURIOUS",ShapeLensConfig::FILTER_SPURIOUS,"whether spurious detection are removed");
  IO::updateFITSKeyword(outfptr,"ADD_BORDER",ShapeLensConfig::ADD_BORDER,"amount of border around detected objects");
  IO::updateFITSKeyword(outfptr,"MIN_PIXELS",ShapeLensConfig::MIN_PIXELS,"minimum number of pixels for an object");
  IO::updateFITSKeyword(outfptr,"MIN_THRESHOLD",ShapeLensConfig::MIN_THRESHOLD,"threshold for MIN_PIXELS");
  IO::updateFITSKeyword(outfptr,"DETECT_THRESHOLD",ShapeLensConfig::DETECT_THRESHOLD,"detection threshold for objects");

  fits_write_record(outfptr,"",&status);
  IO::appendFITSHistory(outfptr,sobj.getHistory());
  
  if (saveErrors) {
    // if noisemodel is GAUSSIAN and shapelet model is untransformed
    // cov is diagonal and constant.
    // so save only 1x1 matrix
    bool diag = true;
    if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
      // check whether cov is diagonal by running along first rows
      for (unsigned int i=0; i < 4; i++) {
	for (unsigned int j=i+1; j < sobj.coeffs.size(); j++) {
	  if (sobj.cov(i,j) != 0) {
	    diag = false;
	    break;
	  }
	}
      }
      if (diag) {
	NumMatrix<data_t> one(1,1);
	one(0,0) = cov(0,0);
	IO::writeFITSImage(outfptr,one,sobj.getName()+"ERRORS");
      }
    }
    if (ShapeLensConfig::NOISEMODEL != "GAUSSIAN" || diag == false)
      IO::writeFITSImage(outfptr,cov,sobj.getName()+"ERRORS");
  }
}

void SIFFile::load(ShapeletObject& sobj, bool preserve_config) {
  int status = 0;
  // check whether filename exists
  ifstream siffile;
  siffile.open(filename.c_str());
  if( !siffile ) {
    cerr << "SIFFile: error opening " + filename << endl;
    siffile.close();
    terminate();
  }
  siffile.close();
  
  // read shapelet coeff from pHDU
  fitsfile *fptr = IO::openFITSFile(filename);
  NumMatrix<data_t> coeffs;
  status = IO::readFITSImage(fptr,coeffs);
  sobj.coeffs.setCoeffs(coeffs);

  // read shapelet parameters
  // make use of friendship of Composite2D and ShapeletObject
  sobj.changeModel = sobj.changeM = true;
  data_t tmp;
  status = IO::readFITSKeyword(fptr,"BETA",tmp);
  sobj.setBeta(tmp);
  status = IO::readFITSKeyword(fptr,"CHI2",sobj.chisquare);
  bool errors;
  status = IO::readFITSKeyword(fptr,"ERRORS",errors);
  unsigned long flags;
  status = IO::readFITSKeyword(fptr,"FLAGS",flags);
  sobj.flags = std::bitset<16>(flags);
  status = IO::readFITSKeywordString(fptr,"EXTNAME",sobj.name);
  if (status != 0)
    sobj.name = "";
  status = IO::readFITSKeyword(fptr,"TAG",sobj.tag);
  if (status != 0)
    sobj.tag = 0;

  // read frame parameters
  status = IO::readFITSKeywordString(fptr,"BASEFILE",sobj.basefilename);
  status = IO::readFITSKeyword(fptr,"ID",sobj.id);
  status = IO::readFITSKeyword(fptr,"CLASSIFIER",sobj.classifier);
  if (status != 0)
    sobj.classifier = 0;

  grid_t xmin,xmax,ymin,ymax;
  status = IO::readFITSKeyword(fptr,"XMIN",xmin);
  status = IO::readFITSKeyword(fptr,"XMAX",xmax);
  status = IO::readFITSKeyword(fptr,"YMIN",ymin);
  status = IO::readFITSKeyword(fptr,"YMAX",ymax);
  sobj.model.grid = Grid(xmin,ymin,xmax-xmin,ymax-ymin);
  complex<data_t> xc;
  status = IO::readFITSKeyword(fptr,"CENTROID",xc);
  sobj.xcentroid(0) = real(xc);
  sobj.xcentroid(1) = imag(xc);
 
  // shapelensconfig parameters
  // set them only if preserve_config == 0
  if (!preserve_config) {
    status = IO::readFITSKeywordString(fptr,"NOISEMODEL",ShapeLensConfig::NOISEMODEL);
    status = IO::readFITSKeyword(fptr," NMAX_LOW",ShapeLensConfig::NMAX_LOW);
    status = IO::readFITSKeyword(fptr,"NMAX_HIGH",ShapeLensConfig::NMAX_HIGH);
    status = IO::readFITSKeyword(fptr," BETA_LOW",ShapeLensConfig::BETA_LOW);
    status = IO::readFITSKeyword(fptr,"BETA_HIGH",ShapeLensConfig::BETA_HIGH);
    status = IO::readFITSKeyword(fptr,"ALLOW_FLATTENING",ShapeLensConfig::ALLOW_FLATTENING);
    status = IO::readFITSKeyword(fptr,"FILTER_SPURIOUS",ShapeLensConfig::FILTER_SPURIOUS);
    status = IO::readFITSKeyword(fptr,"ADD_BORDER",ShapeLensConfig::ADD_BORDER);
    status = IO::readFITSKeyword(fptr,"MIN_PIXELS",ShapeLensConfig::MIN_PIXELS);
    status = IO::readFITSKeyword(fptr,"MIN_THRESHOLD",ShapeLensConfig::MIN_THRESHOLD);
    status = IO::readFITSKeyword(fptr,"DETECT_THRESHOLD",ShapeLensConfig::DETECT_THRESHOLD);
  }

  // read history
  std::string hstr;
  IO::readFITSKeyCards(fptr,"HISTORY",hstr);
  sobj.history.clear();
  sobj.history.setSilent();
  sobj.history << hstr;
  sobj.history.unsetSilent();

  // if errors have been saved, load it
  if (errors) {
    fits_movrel_hdu(fptr,1,IMAGE_HDU,&status);
    status = IO::readFITSImage(fptr,sobj.cov);
    // legacy mode: if errors are coefficient errors instead of full cov. matrix:
    // set them on diagonal of cov
    if (sobj.cov.getRows() == sobj.coeffs.getNMax() + 1) {
      status = IO::readFITSImage(fptr,coeffs);
      CoefficientVector<data_t> errors(coeffs);
      sobj.setErrors(errors);
    }
    // if noise is Gaussian and shapelet model is untransformed
    // expand 1x1 matrix to full size
    else if (sobj.cov.getRows() == 1 && sobj.cov.getColumns() == 1 && sobj.coeffs.size() > 1) {
      data_t sigma = sobj.cov(0,0);
      sobj.cov.resize(sobj.coeffs.size(),sobj.coeffs.size());
      sobj.cov.clear();
      for (unsigned int i=0; i < sobj.coeffs.size(); i++)
	sobj.cov(i,i) = sigma;
    }
  }
  IO::closeFITSFile(fptr);
}
