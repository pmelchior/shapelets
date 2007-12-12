#include <frame/SExFrame.h>
#include <NumVectorMasked.h>
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

SExFrame::SExFrame (std::string fitsfile) : Image<data_t>(fitsfile), weight() {
  catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of underlying Image copied since often used
  axsize0 = Image<data_t>::getSize(0);
  axsize1 = Image<data_t>::getSize(1);

  // for artificial noise
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
}

SExFrame::SExFrame (std::string fitsfile, std::string weightfile) : Image<data_t>(fitsfile), weight(weightfile) {
  catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of underlying Image copied since often used
  axsize0 = Image<data_t>::getSize(0);
  axsize1 = Image<data_t>::getSize(1);

  // for artificial noise
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
}

SExFrame::~SExFrame() {
  gsl_rng_free (r);
}

void SExFrame::readCatalog(std::string catfile) {
  catalog = Catalog(catfile);
  catRead = 1;
}

void SExFrame::readSegmentationMap(std::string segmapfile) {
  segMap = SegmentationMap(segmapfile);
  segmapRead = 1;
  // check if NOISE_MEAN and NOISE_RMS is given as header keyword in segmapfile
  fitsfile* fptr = openFITSFile(segmapfile);
  int status = readFITSKeyword(fptr,"NOISE_MEAN",bg_mean);
  status = readFITSKeyword(fptr,"NOISE_RMS",bg_rms);
  if (status == 0)
    estimatedBG = 1;
}

unsigned int SExFrame::getNumberOfObjects() {
  if (catRead)
    return catalog.size() - 1;
  else
    return 0;
}

const Catalog& SExFrame::getCatalog() {
  return catalog;
}

void SExFrame::fillObject(Object& O) {
  if (!estimatedBG) estimateNoise();
    
  unsigned int id = O.getID();
  int xmin, xmax, ymin, ymax;
  xmin = catalog[id].XMIN;
  xmax = catalog[id].XMAX;
  ymin = catalog[id].YMIN;
  ymax = catalog[id].YMAX;

  O.history << "# Extracting Object " << O.getID() << " (NUMBER = " <<  catalog[id].NUMBER << "), ";
  O.history << "found in the area (" << xmin << "/" << ymin << ") to (";
  O.history << xmax << "/" << ymax << ")" << std::endl;
  
  // check if outer sizes of the object are identical to the image
  // boundary, since then the objects is cutted 
  bool cutflag = 0;
  if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
    O.history << "# Object cut off at the image boundary!" << endl;
    cutflag = 1;
  }
  
  addFrameBorder(ShapeLensConfig::ADD_BORDER, xmin,xmax,ymin,ymax);
  O.history << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
  O.history << xmax << "/" << ymax << ")" << std::endl;

  // check if object was close to the image boundary so that noise has to be added
  if (xmin < 0 || ymin < 0 || xmax >= axsize0 || ymax >= axsize1) {
    if (!cutflag)
      O.history << "# Object close to image boundary: Possible cut-off. Extended area filled with noise." << std::endl;
    else
       O.history << "# Extended area filled with noise." << std::endl;
  }

  // define new object data set, find 1-sigma noise oscillations with more 
  // than 4 pixels and set their pixelmap flag to -2
  // in the end only object data into new vector of smaller size, the rest will
  // filled up with artificial noise
  const NumVector<data_t>& data = Image<data_t>::getData();
  NumVector<data_t>& objdata = O.accessData();
  objdata.resize((xmax-xmin+1)*(ymax-ymin+1));
  SegmentationMap& objSegMap = O.accessSegmentationMap();
  objSegMap.resize((xmax-xmin+1)*(ymax-ymin+1));
  NumVector<data_t>& objWeightMap = O.accessWeightMap();
  if (weight.size()!=0) 
    objWeightMap.resize((xmax-xmin+1)*(ymax-ymin+1));
  vector<uint> nearby_objects;

  for (int i =0; i < objdata.size(); i++) {
    // old coordinates derived from new pixel index i
    int axis0 = xmax-xmin+1;
    int x = i%axis0 + xmin;
    int y = i/axis0 + ymin;
    uint j = Image<data_t>::getGrid().getPixel(x,y);

    // if pixel is out of image region, fill noise from default values
    // since we fill same noise into data and into bgrms
    // the overall chi^2 (normalized by bg_rms) is unaffected by this region
    if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1) {
      objdata(i) = gsl_ran_gaussian (r, bg_rms);
      objSegMap(i) = 0;
      if (weight.size()!=0) 
	objWeightMap(i) = 1./gsl_pow_2(bg_rms);
      if (!subtractBG)
	objdata(i) += bg_mean;
    } 
    //now inside image region
    else {
      // mask other detected objects in the frame
      if ((segMap(j) > 0 && segMap(j) != catalog[id].NUMBER) || (segMap(j) < 0 && ShapeLensConfig::FILTER_SPURIOUS)) {
 	// if we have a weight map 
	if (weight.size()!=0)
	  objdata(i) = gsl_ran_gaussian (r, sqrt(1./weight(j)));
	else
	  objdata(i) = gsl_ran_gaussian (r, bg_rms);
	// this objects has to yet been found to be nearby
	if (std::find(nearby_objects.begin(),nearby_objects.end(),segMap(j)) == nearby_objects.end()) {
	  O.history << "# Object " << segMap(j) << " nearby, but not overlapping." << std::endl;
	  nearby_objects.push_back(segMap(j));
	}
 	if (!subtractBG)
  	  objdata(i) += bg_mean;
      }
      else {
	objdata(i) = data(j);
	if (subtractBG) 
	  objdata(i) -= bg_mean;
      }
      objSegMap(i) = segMap(j);
      if (weight.size()!=0) 
	objWeightMap(i) = weight(j);
    }
  }
    
  // Grid will be changed but not shifted (all pixels stay at their position)
  O.accessGrid() = objSegMap.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);

  // Fill other quantities into Object
  O.history << "# Segment:" << endl;
  O.computeFluxCentroid();
  O.setFlux(catalog[id].FLUX);
  Point2D centroid(catalog[id].XCENTROID,catalog[id].YCENTROID);
  O.setCentroid(centroid);
  O.setDetectionFlag(catalog[id].FLAGS);
  O.setStarGalaxyProbability(catalog[id].CLASSIFIER);
  O.setBlendingProbability(computeBlendingProbability(id));
  O.setBaseFilename(Image<data_t>::getFilename());
  O.setNumber(catalog[id].NUMBER);
  data_t obj_bg_mean = 0;
  if (!subtractBG)
    obj_bg_mean = bg_mean;
  O.setNoiseMeanRMS(obj_bg_mean,bg_rms);
}

void SExFrame::subtractBackground() {
  subtractBG = 1;
}

// now extend to region around the object by
// typically objectsize/2, minimum 12 pixels
void SExFrame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
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
data_t SExFrame::computeBlendingProbability(unsigned int nr) {
  unsigned short flag = catalog[nr].FLAGS;
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
    std::cerr << "SExFrame: provide segmentation map before calling subtractNoise()!" << std::endl;
    std::terminate();
  }
  // set noise estimates
  // since the position of the object is known from segMap, we can compute the
  // noise now on all pixels not in pixellist
  // 1) set mask(i)=1, when segMap(i) != 0
  // 2) create NumVectorMasked from objdata and mask
  // 3) compute std from that
  const NumVector<data_t>& data = Image<data_t>::getData();
  NumVector<bool> mask(data.size());
  for(int i =0; i < data.size(); i++)
    if (segMap(i) != 0) 
      mask(i) = 1;
  NumVectorMasked<data_t> masked(data,mask);
  masked.kappa_sigma_clip(bg_mean,bg_rms);
  estimatedBG = 1;
}

data_t SExFrame::getNoiseMean() {
  if (!estimatedBG) estimateNoise();
  return bg_mean;
}

data_t SExFrame::getNoiseRMS() {
  if (!estimatedBG) estimateNoise();
  return bg_rms;
}

void SExFrame::setNoiseMeanRMS(data_t mean, data_t rms) {
  estimatedBG = 1;
  bg_mean = mean;
  bg_rms = rms;
}
