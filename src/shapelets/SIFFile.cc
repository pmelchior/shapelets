#include <shapelets/SIFFile.h>
#include <IO.h>
#include <ShapeLensConfig.h>
#include <bitset>

using namespace std;

SIFFile::SIFFile(std::string infilename) {
  filename = infilename;
}

void SIFFile::save(const ShapeletObject& sobj) {
  fitsfile *outfptr = createFITSFile(filename);
  saveSObj(outfptr, sobj);
  closeFITSFile(outfptr);
}

void SIFFile::save(const ShapeletObject& sobj, const ShapeletObject& unreg) {
  fitsfile *outfptr = createFITSFile(filename);
  saveSObj(outfptr, sobj);
  saveSObj(outfptr, unreg);
  closeFITSFile(outfptr);
}

void SIFFile::save(ShapeletObjectList& sl) {
  fitsfile *outfptr = createFITSFile(filename);
  for(ShapeletObjectList::iterator iter = sl.begin(); iter != sl.end(); iter++) 
    saveSObj(outfptr, *(*iter));
  closeFITSFile(outfptr);
}

void SIFFile::saveSObj(fitsfile* outfptr, const ShapeletObject& sobj) {
  int status = writeFITSImage(outfptr,sobj.getCoeffs());

  // when present, save errors also
  bool saveErrors = 0;
  const NumMatrix<data_t>& errors = sobj.getDecompositionErrors();
  if (errors.getRows() != 0 && errors.getColumns() != 0)
    saveErrors = 1;
  

  // add shapelet parameters and other necessary information in pHDU header
  updateFITSKeyword(outfptr,"VERSION",(unsigned int) 1,"SIF version");

  // ** Shapelet parameters **
  fits_write_record(outfptr,"",&status);
  fits_write_record(outfptr,"        / Shapelet parameters  /",&status);
  updateFITSKeyword(outfptr,"BETA",sobj.getBeta(),"scale size in pixel units");
  updateFITSKeyword(outfptr,"DIM",sobj.getCoeffs().getRows(),"dimensions in shapelet space (nmax+1)");
  updateFITSKeyword(outfptr,"CHI2",sobj.getDecompositionChiSquare(),"decomposition quality");
  updateFITSKeyword(outfptr,"ERRORS",saveErrors,"whether coefficient errors are available");
  updateFITSKeyword(outfptr,"R",sobj.getRegularizationR(),"negative flux / positive flux");
  updateFITSKeyword(outfptr,"FLAGS",sobj.getFlags().to_ulong(),"extraction and decomposition flags");
  updateFITSKeywordString(outfptr, "EXTNAME", sobj.getName(), "shapelet object name");
  updateFITSKeyword(outfptr, "TAG", sobj.getTag(), "shapelet object tag");

  // ** Frame parameters **
  fits_write_record(outfptr,"        / Frame parameters     /",&status);
  updateFITSKeywordString(outfptr,"BASEFILE",sobj.getBaseFilename(),"originating data file");
  updateFITSKeyword(outfptr,"ID",sobj.getObjectID(),"object id in BASEFILE");
  updateFITSKeyword(outfptr,"NR",sobj.getObjectNumber(),"object nr in BASEFILE");
  updateFITSKeyword(outfptr,"CLASSIFIER",sobj.getObjectClassifier(),"object classifier");
  updateFITSKeyword(outfptr,"NOISE_MEAN",sobj.getNoiseMean(),"mean of pixel noise");
  updateFITSKeyword(outfptr,"NOISE_RMS",sobj.getNoiseRMS(),"rms of pixel noise");
  const Grid& grid = sobj.getGrid();
  updateFITSKeyword(outfptr,"XMIN",grid.getStartPosition(0),"min(X) in image pixels");
  updateFITSKeyword(outfptr,"XMAX",grid.getStopPosition(0),"max(X) in image pixels");
  updateFITSKeyword(outfptr,"YMIN",grid.getStartPosition(1),"min(Y) in image pixels");
  updateFITSKeyword(outfptr,"YMAX",grid.getStopPosition(1),"min(Y) in image pixels");
  const Point2D& centroid = sobj.getCentroid();
  complex<data_t> xc(centroid(0),centroid(1));
  updateFITSKeyword(outfptr,"CENTROID",xc,"centroid position in image pixels");

  // ** ShapeLensConfig parameters **
  fits_write_record(outfptr,"        / ShapeLensConfig parameters /",&status);
  updateFITSKeywordString(outfptr,"NOISEMODEL",ShapeLensConfig::NOISEMODEL,"noise model");
  updateFITSKeyword(outfptr," NMAX_LOW",ShapeLensConfig::NMAX_LOW,"lower bound for n_max");
  updateFITSKeyword(outfptr,"NMAX_HIGH",ShapeLensConfig::NMAX_HIGH,"upper bound for n_max");
  updateFITSKeyword(outfptr," BETA_LOW",ShapeLensConfig::BETA_LOW,"lower bound for beta");
  updateFITSKeyword(outfptr,"BETA_HIGH",ShapeLensConfig::BETA_HIGH,"upper bound for beta");
  updateFITSKeyword(outfptr,"REGULARIZE",ShapeLensConfig::REGULARIZE,"whether regularization has been employed");
  updateFITSKeyword(outfptr,"REG_LIMIT",ShapeLensConfig::REG_LIMIT,"demanded upper limit for R");
  updateFITSKeyword(outfptr,"SAVE_UNREG",ShapeLensConfig::SAVE_UNREG,"whether unregularized model is saved");
  updateFITSKeyword(outfptr,"ALLOW_FLATTENING",ShapeLensConfig::ALLOW_FLATTENING,"whether flattening of chi^2 is allowed");
  updateFITSKeyword(outfptr,"FILTER_SPURIOUS",ShapeLensConfig::FILTER_SPURIOUS,"whether spurious detection are removed");
  updateFITSKeyword(outfptr,"ADD_BORDER",ShapeLensConfig::ADD_BORDER,"amount of border around detected objects");
  updateFITSKeyword(outfptr,"MIN_PIXELS",ShapeLensConfig::MIN_PIXELS,"minimum number of pixels for an object");
  updateFITSKeyword(outfptr,"MIN_THRESHOLD",ShapeLensConfig::MIN_THRESHOLD,"threshold for MIN_PIXELS");
  updateFITSKeyword(outfptr,"DETECT_THRESHOLD",ShapeLensConfig::DETECT_THRESHOLD,"detection threshold for objects");

  fits_write_record(outfptr,"",&status);
  appendFITSHistory(outfptr,sobj.getHistory());
  
  if (saveErrors)
    writeFITSImage(outfptr,errors,sobj.getName()+"ERRORS");

}

void SIFFile::load(ShapeletObject& sobj, bool preserve_config) {
  int status = 0;
  fitsfile *fptr;
  // open pHUD of fits file
  fits_open_file(&fptr,filename.c_str(),READONLY,&status);
  // read shapelet coeff from pHDU
  status = readFITSImage(fptr,sobj.coeffs);

  // read shapelet parameters
  // make use of friendship of Composite2D and ShapeletObject
  sobj.change = 1;
  data_t tmp;
  status = readFITSKeyword(fptr,"BETA",tmp);
  sobj.setBeta(tmp);
  status = readFITSKeyword(fptr,"CHI2",sobj.chisquare);
  bool errors;
  status = readFITSKeyword(fptr,"ERRORS",errors);
  status = readFITSKeyword(fptr,"R",sobj.R);
  unsigned long flags;
  status = readFITSKeyword(fptr,"FLAGS",flags);
  sobj.flags = std::bitset<16>(flags);
  status = readFITSKeywordString(fptr,"EXTNAME",sobj.name);
  if (status != 0)
    sobj.name = "";
  status = readFITSKeyword(fptr,"TAG",sobj.tag);
  if (status != 0)
    sobj.tag = 0;

  // read frame parameters
  status = readFITSKeywordString(fptr,"BASEFILE",sobj.basefilename);
  status = readFITSKeyword(fptr,"ID",sobj.id);
  status = readFITSKeyword(fptr,"NR",sobj.nr);
  status = readFITSKeyword(fptr,"CLASSIFIER",sobj.classifier);
  if (status != 0)
    sobj.classifier = 0;

  status = readFITSKeyword(fptr,"NOISE_MEAN",sobj.noise_mean);
  status = readFITSKeyword(fptr,"NOISE_RMS",sobj.noise_rms);
  data_t xmin,xmax,ymin,ymax;
  status = readFITSKeyword(fptr,"XMIN",xmin);
  status = readFITSKeyword(fptr,"XMAX",xmax);
  status = readFITSKeyword(fptr,"YMIN",ymin);
  status = readFITSKeyword(fptr,"YMAX",ymax);
  sobj.grid = Grid(xmin,xmax,1,ymin,ymax,1);
  complex<data_t> xc;
  status = readFITSKeyword(fptr,"CENTROID",xc);
  sobj.xcentroid(0) = real(xc);
  sobj.xcentroid(1) = imag(xc);
  
  // shapelensconfig parameters
  // set them only if preserve_config == 0
  if (!preserve_config) {
    status = readFITSKeywordString(fptr,"NOISEMODEL",ShapeLensConfig::NOISEMODEL);
    status = readFITSKeyword(fptr," NMAX_LOW",ShapeLensConfig::NMAX_LOW);
    status = readFITSKeyword(fptr,"NMAX_HIGH",ShapeLensConfig::NMAX_HIGH);
    status = readFITSKeyword(fptr," BETA_LOW",ShapeLensConfig::BETA_LOW);
    status = readFITSKeyword(fptr,"BETA_HIGH",ShapeLensConfig::BETA_HIGH);
    status = readFITSKeyword(fptr,"REGULARIZE",ShapeLensConfig::REGULARIZE);
    status = readFITSKeyword(fptr,"REG_LIMIT",ShapeLensConfig::REG_LIMIT);
    status = readFITSKeyword(fptr,"SAVE_UNREG",ShapeLensConfig::SAVE_UNREG);
    status = readFITSKeyword(fptr,"ALLOW_FLATTENING",ShapeLensConfig::ALLOW_FLATTENING);
    status = readFITSKeyword(fptr,"FILTER_SPURIOUS",ShapeLensConfig::FILTER_SPURIOUS);
    status = readFITSKeyword(fptr,"ADD_BORDER",ShapeLensConfig::ADD_BORDER);
    status = readFITSKeyword(fptr,"MIN_PIXELS",ShapeLensConfig::MIN_PIXELS);
    status = readFITSKeyword(fptr,"MIN_THRESHOLD",ShapeLensConfig::MIN_THRESHOLD);
    status = readFITSKeyword(fptr,"DETECT_THRESHOLD",ShapeLensConfig::DETECT_THRESHOLD);
  }

  // read history
  std::string hstr;
  readFITSKeyCards(fptr,"HISTORY",hstr);
  sobj.history.clear();
  sobj.history.setSilent();
  sobj.history << hstr;
  sobj.history.unsetSilent();

  // if errors have been saved, load it
  if (errors) {
    fits_movnam_hdu(fptr,IMAGE_HDU, "ERRORS",0, &status);
    status = readFITSImage(fptr,sobj.errors);
  }

  fits_close_file(fptr, &status);

}
