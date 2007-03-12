#include <SExFrame.h>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <set>
#include <vector>
using namespace std;
using namespace boost;

typedef unsigned int uint;

SExFrame::SExFrame (std::string fitsfile) : FitsImage<double>(fitsfile) {
  SExCatFormat empty = {0,0,0,0,0,0,0,0,0,0};
  sf = empty;
  catChecked = catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of unserlying FitsImage copied since often used
  axsize0 = FitsImage<double>::getSize(0);
  axsize1 = FitsImage<double>::getSize(1);

  text << "# Reading FITS file " << fitsfile << endl;
  text << "# Image properties: size = "<< axsize0 << "/" << axsize1 << std::endl; 
  history.append(text);

  // for artificial noise
  const gsl_rng_type * T;
  //gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
}

SExFrame::~SExFrame() {
  gsl_rng_free (r);
}

void SExFrame::readCatalog(std::string catfile) {
  // first inser empty object 0, since SExtractor starts with NUMBER 1
  objectList.clear();
  SExCatObject s0 = {0,0,0,0,0,0,0,0,0,0};
  objectList.push_back(s0);
  // open cat file
  ifstream catalog (catfile.c_str());
  if (catalog.fail()) {
    history.append("SExFrame: catalog file does not exist!\n");
    terminate();
  }
  catalog.clear();
  // read in cat file
  // 1) parse the format definition lines
  // 2) fill object information into SExCatObjects -> objectList;
  string line;
  while(getline(catalog, line)) {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    // split entries at empty chars
    boost::char_separator<char> sep(" ");
    Tok tok(line, sep);
    // first of all we copy the token into string vector
    // though this is not too sophisticated it does not hurt and is more 
    // convenient
    std::vector<std::string> column;
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
      column.push_back(*tok_iter);
    // comment line: contains format definition at columns 2,3
    if (column[0].compare("#") == 0)
      insertFormatField(column[2],column[1]);
    // at this point we should have a complete format definition
    // check it!
    // from here on we expect the list of object to come along.
    else {
      if (!catChecked) checkFormat();
      // then set up a true SExCatObject
      SExCatObject so;
      so.NUMBER = (unsigned int) atoi(column[sf.NUMBER-1].c_str());
      // the sextractor corrdinates start with (1,1), ours with (0,0)
      so.XMIN_IMAGE = atoi(column[sf.XMIN_IMAGE-1].c_str())-1;
      so.XMAX_IMAGE = atoi(column[sf.XMAX_IMAGE-1].c_str())-1;
      so.YMIN_IMAGE = atoi(column[sf.YMIN_IMAGE-1].c_str())-1;
      so.YMAX_IMAGE = atoi(column[sf.YMAX_IMAGE-1].c_str())-1;
      so.XWIN_IMAGE = atof(column[sf.XWIN_IMAGE-1].c_str())-1;
      so.YWIN_IMAGE = atof(column[sf.YWIN_IMAGE-1].c_str())-1;
      so.FLUX_AUTO = atof(column[sf.FLUX_AUTO-1].c_str())-1;
      so.FLAGS = (unsigned char) atoi(column[sf.FLAGS-1].c_str());
      so.CLASS_STAR = (double) atof(column[sf.CLASS_STAR-1].c_str());
     // then push it on objectList
      objectList.push_back(so);
    }
  }
  catalog.close();
  catRead = 1;
}

void SExFrame::readSegmentationMap(std::string segmentfile) {
  FitsImage<int>* seg = new FitsImage<int>(segmentfile);
  segMap = seg->getData();
  delete seg;
  segmapRead = 1;
}

unsigned int SExFrame::getNumberOfObjects() {
  if (catRead)
    return objectList.size() - 1;
  else
    return 0;
}


void SExFrame::fillObject(Object& O) {
  if (!estimatedBG) estimateNoise();
  unsigned int nr = O.getID();
  int xmin, xmax, ymin, ymax;
  xmin = objectList[nr].XMIN_IMAGE;
  xmax = objectList[nr].XMAX_IMAGE;
  ymin = objectList[nr].YMIN_IMAGE;
  ymax = objectList[nr].YMAX_IMAGE;

  O.history = history;
  text << "# Extracting Object " << O.getID() << " (NUMBER = " <<  objectList[nr].NUMBER << "), ";
  text << "found in the area (" << xmin << "/" << ymin << ") to (";
  text << xmax << "/" << ymax << ")" << std::endl;
  O.history.append(text);
  
  // check if outer sizes of the object are identical to the image
  // boundary, since then the objects is cutted 
  bool cutflag = 0;
  if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
    O.history.append("# Object cut off at the image boundary!\n");
    cutflag = 1;
  }
  
  addFrameBorder(0.5, xmin,xmax,ymin,ymax);
  text << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
  text << xmax << "/" << ymax << ")" << std::endl;
  O.history.append(text);

  // check if object was close to the image boundary so that noise has to be added
  if (xmin < 0 || ymin < 0 || xmax >= axsize0 || ymax >= axsize1) {
    if (!cutflag)
      text << "# Object close to image boundary: Possible cut-off. Extended area filled with noise." << std::endl;
    else
       text << "# Extended area filled with noise." << std::endl;
  }
  O.history.append(text);

  // define new object data set, find 1-sigma noise oscillations with more 
  // than 4 pixels and set their pixelmap flag to -2
  // in the end only object data into new vector of smaller size, the rest will
  // filled up with artificial noise
  NumVector<double>& objdata = O.accessData();
  objdata.resize((xmax-xmin+1)*(ymax-ymin+1));
  const NumVector<double>& data = FitsImage<double>::getData();

  for (int i =0; i < objdata.size(); i++) {
    // old coordinates derived from new pixel index i
    int axis0 = xmax-xmin+1;
    int x = i%axis0 + xmin;
    int y = i/axis0 + ymin;
    uint j = FitsImage<double>::getPixel(x,y);

    // if pixel is out of image region, fill noise from default values
    // since we fill same noise into data and into bgrms
    // the overall chi^2 (normalized by bg_rms) is unaffected by this region
    if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1) {
      objdata(i) = gsl_ran_gaussian (r, bg_rms);
      if (!subtractBG)
	objdata(i) += bg_mean;
    } 
    //now inside image region
    else {
      // mask other detected objects in the frame
      if (segMap(j) != objectList[nr].NUMBER && segMap(j) > 0) {
 	objdata(i) = gsl_ran_gaussian (r, bg_rms);
 	if (!subtractBG)
  	  objdata(i) += bg_mean;
      }
      else {
	objdata(i) = data(j);
	if (subtractBG) 
	  objdata(i) -= bg_mean;
      }
    }
  }
    
  // Grid will be changed but not shifted (all pixels stay at their position)
  O.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);
  
  // Fill other quantities into Object
  O.setNoiseMeanRMS(bg_mean,bg_rms);
  O.setNoiseModel("GAUSSIAN");
  O.computeFluxCentroid();
  O.setFlux(objectList[nr].FLUX_AUTO);
  Point2D centroid(objectList[nr].XWIN_IMAGE,objectList[nr].YWIN_IMAGE);
  O.setCentroid(centroid);

  O.setDetectionFlag(objectList[nr].FLAGS);
  O.setStarGalaxyProbability(objectList[nr].CLASS_STAR);
  O.setBlendingProbability(computeBlendingProbability(nr));
O.setBaseFilename(FitsImage<double>::getFilename());
}

void SExFrame::subtractBackground() {
  subtractBG = 1;
}

// fill in the column position into the SExCatFormat struct
void SExFrame::insertFormatField(std::string type, std::string columnnr) {
  unsigned int colnr = atoi(columnnr.c_str());
  if (type.compare("NUMBER")==0)
    sf.NUMBER = colnr;
  else if (type.compare("XMIN_IMAGE")==0)
    sf.XMIN_IMAGE = colnr;
  else if (type.compare("XMAX_IMAGE")==0)
    sf.XMAX_IMAGE = colnr;
  else if (type.compare("YMIN_IMAGE")==0)
    sf.YMIN_IMAGE = colnr;
  else if (type.compare("YMAX_IMAGE")==0)
    sf.YMAX_IMAGE = colnr;
  else if (type.compare("XWIN_IMAGE")==0)
    sf.XWIN_IMAGE = colnr;
  else if (type.compare("YWIN_IMAGE")==0)
    sf.YWIN_IMAGE = colnr;
  else if (type.compare("FLUX_AUTO")==0)
    sf.FLUX_AUTO = colnr;
  else if (type.compare("FLAGS")==0)
    sf.FLAGS = colnr;
  else if (type.compare("CLASS_STAR")==0)
    sf.CLASS_STAR = colnr;
}

// check if all necessary keywords are present
void SExFrame::checkFormat() {
  bool trouble = 0;
  if (sf.NUMBER == 0) {
    text << "SExFrame: NUMBER keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.XMIN_IMAGE == 0) {
    text << "SExFrame: XMIN_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.XMAX_IMAGE == 0) {
    text << "SExFrame: XMAX_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.YMIN_IMAGE == 0) {
    text << "SExFrame: YMIN_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.YMAX_IMAGE == 0) {
    text << "SExFrame: YMAX_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.XWIN_IMAGE == 0) {
    text << "SExFrame: XWIN_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.YWIN_IMAGE == 0) {
    text << "SExFrame: YWIN_IMAGE keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.FLUX_AUTO == 0) {
    text << "SExFrame: FLUX_AUTO keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.FLAGS == 0) {
    text << "SExFrame: FLAGS keyword not provided!" << std::endl;
    trouble = 1;
  }
  if (sf.CLASS_STAR == 0) {
    text << "SExFrame: CLASS_STAR keyword not provided!" << std::endl;
    trouble = 1;
  }
  history.append(text);

  if (trouble)
    terminate();
  else
    catChecked = 1;
}

// now extend to region around the object by
// typically objectsize/2, minimum 12 pixels
void SExFrame::addFrameBorder(double factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
  int xrange, yrange, xborder, yborder;
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  switch(xrange%2) {
  case 1: xmax++; break;
  }
  switch(yrange%2) {
  case 1: ymax++; break;
  }
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  // make the object frame square, because of const beta in both directions
  if (xrange < yrange) {
    yborder = GSL_MAX_INT((int)floor(yrange*factor), 12);
    xborder = yborder + (yrange - xrange)/2;
  } else {
    xborder = GSL_MAX_INT((int)floor(xrange*factor), 12);
    yborder = xborder + (xrange - yrange)/2;
  }
  xmin -= xborder;
  xmax += xborder-1;
  ymin -= yborder;
  ymax += yborder-1;
}

// compute probability of an object being blended with another one
// simple: if SExtractor FLAG contains 2, probability equals 1, else 0
double SExFrame::computeBlendingProbability(unsigned int nr) {
  unsigned short flag = objectList[nr].FLAGS;
  // if 2 is in FLAGS: object flaged as blended
  if ((flag >> 1)%2 == 1) 
    return 1;
  else 
    return 0;
}

const NumVector<int>& SExFrame::getObjectMap() {
  return segMap;
}


// estimate noise by iterative sigma clipping
void SExFrame::estimateNoise() {
  if (!segmapRead) {
    std::cout << "SExFrame: provide segmentation map before calling subtractNoise()!" << std::endl;
    std::terminate();
  }
   
  // for GSL sorting functions we have to copy NumVector to double*
  // only copy the pixel values where not object is found in the segmentation map
  int npixels = FitsImage<double>::getNumberOfPixels(), jmax = 0;
  double* D = (double *) malloc(npixels*sizeof(double));
  for (int i=0; i < npixels; i++) {
    if (segMap(i) == 0) {
      D[jmax] = (FitsImage<double>::getData())(i);
      jmax++;
    }
  }

  // first check left border of the image for pixel value variations
  // if its 0, its a (simulated) image with 0 noise
  // background_variance is set to the minmal value above 0
  if (gsl_stats_variance(D,1,axsize0-1) == 0) {
    bg_mean = 0;
    double min = gsl_stats_max(D,1,npixels);
    for (int i=0; i < npixels; i++)
      if (D[i] < min && D[i] > 0) min = D[i];
    bg_rms = sqrt(min*min);
  } 
  else {
    gsl_sort(D, 1, jmax-1);
    double sigma = gsl_stats_sd(D,1,jmax-1);
    double median = gsl_stats_median_from_sorted_data (D,1,jmax-1);

    // sigma clipping here: only using pixel not considered as object here
    int j;
    while (1) {
      j=0;
       for (int i = 0; i < jmax; i++ ) {
	 // if data is 3 sigma arround iterative median, keep it
	 // maybe additional constrain: D[i] > 0?
	 if (D[i] > median-3*sigma && D[i] < median+3*sigma)  {
	   D[j] = D[i];
	   j++;
	 }
       }
      // next time only work on the first jmax = j, all others are not sorted
      if (j >= 1) {
	gsl_sort(D, 1, j-1);
 	median = gsl_stats_median_from_sorted_data (D,1,j-1);
 	sigma = gsl_stats_sd(D,1,j-1);
	// no change in the selected pixels: converged
	// this is not as fast as 
	// if (fabs(newmedian-median) < median/npixels)
	if (jmax == j)
	  break;
	else {
	  jmax = j;
	}
      }
      else {
	history.append("# Sky background estimation did not converge!\n");
      }
    }
    // if distribution becomes definitely skewed to the brighter values
    // create histogram of pixel values between median+-3*sigma
    // if the pixels to the right become too bright (+1 error = sqrt(count))
    // set it to the values of the appropriate left bin
    // to symmetrize the histogram around the mean.
    if(gsl_stats_skew(D,1,j-1) > 0.1) {
      int nbins = j/100;
      if (nbins%2 ==1) nbins++;
      gsl_histogram * h = gsl_histogram_alloc (nbins);
      gsl_histogram_set_ranges_uniform (h,median-3*sigma,median+3*sigma);
      for(int i =0; i<j; i++)
	gsl_histogram_increment (h, D[i]);
      for(int binnr=nbins/2 +1; binnr<nbins; binnr++) {
	if (gsl_histogram_get(h,binnr) > gsl_histogram_get(h,nbins-binnr)+sqrt(gsl_histogram_get(h,nbins-binnr)))
	  h->bin[binnr] = gsl_histogram_get(h,nbins-binnr);
      }
      bg_mean = gsl_histogram_mean(h);
      bg_rms = gsl_histogram_sigma(h);
      gsl_histogram_free (h);
    }
    // no skewness, image contains enough noise for a solid noise estimation
    // using sigma-clipping only
    else {
      bg_rms = sigma;
      bg_mean = median;
    }
  }
  free(D);
  text << "# Background estimation: mean = " << bg_mean;
  text << ", sigma = " << bg_rms << std::endl;
  history.append(text);
  estimatedBG = 1;
}

double SExFrame::getNoiseMean() {
  if (!estimatedBG) estimateNoise();
  return bg_mean;
}

double SExFrame::getNoiseRMS() {
  if (!estimatedBG) estimateNoise();
  return bg_rms;
}

void SExFrame::setNoiseMeanRMS(double mean, double rms) {
  estimatedBG = 1;
  bg_mean = mean;
  bg_rms = rms;
  text << "# Noise estimates explicitly set: mean = " << mean << ", rms = " << rms << std::endl;
  history.append(text);
}

const History& SExFrame::getHistory () {
  return history;
}
