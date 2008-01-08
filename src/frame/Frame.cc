#include <frame/Frame.h>
#include <NumVectorMasked.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <algorithm>
#include <bitset>

using namespace std;

typedef unsigned int uint;

Frame::Frame(string filename) : Image<data_t>(filename), weight(), segMap(*this), history(segMap.accessHistory()) {
  history << "# Reading FITS file " << filename << endl;
  history << "# Image properties: size = "<< Image<data_t>::getSize(0) << "/" << Image<data_t>::getSize(1) << endl; 
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

Frame::Frame(string datafile, string weightfile) : Image<data_t>(datafile), weight(weightfile), segMap(*this, weight), history(segMap.accessHistory()) {
  history << "# Reading data from " << datafile << " and weights from " << weightfile << endl;
  history << "# Image properties: size = "<< Image<data_t>::getSize(0) << "/" << Image<data_t>::getSize(1) << endl; 
  //weight = Image<data_t>(weightfile);
  if (weight.size() != (*this).size()) {
    history << "Frame: weight map has different layout than data!" << endl;
    terminate();
  } 
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

// estimate noise by iterative sigma clipping
void Frame::estimateNoise() {
  Image<data_t>::getData().kappa_sigma_clip(noise_mean,noise_rms);
  history << "# Background estimation: mean = " << noise_mean;
  history << ", sigma = " << noise_rms << std::endl;
  estimatedBG = 1;
}

data_t Frame::getNoiseMean() {
  if (!estimatedBG) estimateNoise();
  return noise_mean;
}

data_t Frame::getNoiseRMS() {
  if (!estimatedBG) estimateNoise();
  return noise_rms;
}

void Frame::setNoiseMeanRMS(data_t mean, data_t rms) {
  estimatedBG = 1;
  noise_mean = mean;
  noise_rms = rms;
  history << "# Noise estimates explicitly set: mean = " << mean << ", rms = " << rms << std::endl;
}
  
void Frame::subtractBackground() {
  if (!estimatedBG) estimateNoise();
  if (!subtractedBG) {
    for (int i=0; i < Image<data_t>::size(); i++) {
      (Image<data_t>::accessData())(i) -= noise_mean;
    }
    subtractedBG = 1;
    history << "# Background subtraction: noise level = " << noise_mean << std::endl;
    noise_mean = 0;
  }
}

// return threshold for findObject
// if weight map is provided, uses rms from this
data_t Frame::getThreshold(unsigned int pixel, data_t factor) {
  if (weight.size()==0) {
    if (!estimatedBG) estimateNoise();
    return noise_mean + factor*noise_rms;
  }
  else
    return noise_mean + factor*sqrt(1./weight(pixel));
}

void Frame::findObjects() {
  const NumVector<data_t>& data = Image<data_t>::getData();
  unsigned long counter = 0;
  unsigned int npixels = Image<data_t>::size();

  // clean stuff from previous runs of this function
  if (getNumberOfObjects() > 0) {
    objectsPixels.clear();
    catalog.clear();
    for (int i=0; i<npixels; i++)
      segMap(i) = 0;
  }  

  // set up pixellist and objectsPixels
  list<uint> pixellist;
  list<uint>::iterator iter;

  // reset catalog and create dummy CatObject
  CatObject co = {0,0,0,0,0,0,0,0,0};

  // look for all positive fluctuations with at least 1 pixel above
  // highThreshold and minPixels above significanceThreshold
  data_t highThreshold;
  for (int i =0; i < npixels; i++) {
    highThreshold = getThreshold(i,ShapeLensConfig::DETECT_THRESHOLD);
    if (data(i) > highThreshold && segMap(i) == 0) {
      counter++;
      segMap.linkPixelsSetMap (pixellist,i,counter,ShapeLensConfig::MIN_THRESHOLD,noise_mean,noise_rms,1);
      if (pixellist.size() >= ShapeLensConfig::MIN_PIXELS) {
	history << "# Object " << counter << " detected with " << pixellist.size() << " significant pixels at (" << i%(Image<data_t>::getSize(0)) << "/" << i/(Image<data_t>::getSize(0)) << ")"  << std::endl;
	objectsPixels[counter]= pixellist;
	catalog[counter] = co;
      }
      else {
	for(iter = pixellist.begin(); iter != pixellist.end(); iter++ )
	  segMap(*iter) = 0;
	counter--;
      }
    }
  }

//   // improve noise estimates (when weight map is not given)
//   // since the position of all object is known now, we can compute the
//   // noise now on all pixels not associated to an object
//   // 1) set mask(i)=1, when i is in one of the pixellists
//   // 2) create NumVectorMasked from data and mask
//   // 3) compute std from that
//   if (weight.size()==0) {
//     NumVector<bool> mask(data.size());
//     for (int i=0; i<numberofObjects; i++) {
//       list<uint>& pixellist = objectsPixels[i];
//       for(list<uint>::iterator iter = pixellist.begin(); iter != pixellist.end(); iter++ )
// 	mask(*iter) = 1;
//     }
//     NumVectorMasked<data_t> masked(data,mask);
//     masked.kappa_sigma_clip(noise_mean,noise_rms);
//     history << "# Improved background estimation (objects masked):";
//     history << " mean = " << noise_mean << ", sigma = " << noise_rms << endl;
//     history.append(history);
//   }
}

unsigned long Frame::getNumberOfObjects() {
  return catalog.size();
}

// cut the image to a small region around the object
// and set all pixels to zero than are not related to the image
void Frame::fillObject(Object& O) {
  // check if object is in the catalog
  // this will also provied us with the correct entry of catalog (if present)
  Catalog::iterator catiter = catalog.find(O.getID());
  if (catiter != catalog.end()) {
    // set up detection flags bitset
    std::bitset<8> flags(0);

    // for artificial noise for area contaminated by different object
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // define the points (xmin/ymin) and (xmax/ymax) 
    // that enclose the object
    list<uint>& pixellist = objectsPixels[(*catiter).first];
    int axsize0, axsize1,xmin, xmax, ymin, ymax;
    xmin = axsize0 = Image<data_t>::getSize(0);
    xmax = 0;
    ymin = axsize1 = Image<data_t>::getSize(1);
    ymax = 0;
    // loop over all pixels of current object
    for(list<uint>::iterator iter = pixellist.begin(); iter != pixellist.end(); iter++ ) {
	uint x = (*iter)%axsize0;
	uint y = (*iter)/axsize0;
	if (x < xmin) xmin = x;
	if (y < ymin) ymin = y;
	if (x > xmax) xmax = x;
	if (y > ymax) ymax = y;
    }

    O.history << "# Extracting Object " << (*catiter).first;
    O.history << " found in the area (" << xmin << "/" << ymin << ") to (";
    O.history << xmax << "/" << ymax << ")" << endl;
    
    // check if outer sizes of the object are identical to the image
    // boundary, since then the objects is cutted 
    if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
      flags[3] = 1;
      O.history << "# Object cut off at the image boundary!" << endl;
    }

    // add border around object for including the edges
    addFrameBorder(ShapeLensConfig::ADD_BORDER,xmin,xmax,ymin,ymax);
    O.history << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
    O.history << xmax << "/" << ymax << ")" << endl;

    // check if object was close to the image boundary so that noise has to be added
    if (xmin < 0 || ymin < 0 || xmax >= axsize0 || ymax >= axsize1) {
      flags[2] = 1;
      O.history << "# Object close to image boundary: Possible cut-off. Extending frame with noise." << endl;
    }

    // fill the object pixel data
    const NumVector<data_t>& data = Image<data_t>::getData();
    NumVector<data_t>& objdata = O.accessData();
    objdata.resize((xmax-xmin+1)*(ymax-ymin+1));
    SegmentationMap& objSegMap = O.accessSegmentationMap();
    History& objSegMapHistory = objSegMap.accessHistory();
    objSegMapHistory.setSilent();
    objSegMapHistory << history.str();
    objSegMapHistory.unsetSilent();
    objSegMap.resize((xmax-xmin+1)*(ymax-ymin+1));
    NumVector<data_t>& objWeightMap = O.accessWeightMap();
    if (weight.size()!=0) 
      objWeightMap.resize((xmax-xmin+1)*(ymax-ymin+1));
    vector<uint> nearby_objects;

    // lop over all object pixels
    for (int i =0; i < objdata.size(); i++) {
      // old coordinates derived from new pixel index i
      int axis0 = xmax-xmin+1;
      int x = i%axis0 + xmin;
      int y = i/axis0 + ymin;
      uint j = Image<data_t>::getGrid().getPixel(x,y);

      // if pixel is out of image region, fill noise
      if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1) {
	objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	objSegMap(i) = 0;
	if (weight.size()!=0) 
	  objWeightMap(i) = 1./gsl_pow_2(noise_rms);
      } 
      else {
	// filter other objects in the frame
	if ((segMap(j) > 0 && segMap(j) != (*catiter).first) || (segMap(j) < 0 && ShapeLensConfig::FILTER_SPURIOUS)) {
	  // if we have a weight map 
	  if (weight.size()!=0)
	    objdata(i) = noise_mean + gsl_ran_gaussian (r, sqrt(1./weight(j)));
	  else
	    objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	  flags[0] = 1;
	  // this objects has to yet been found to be nearby
	  if (std::find(nearby_objects.begin(),nearby_objects.end(),segMap(j)) == nearby_objects.end()) {
	    O.history << "# Object " << segMap(j) << " nearby, but not overlapping." << std::endl;
	    nearby_objects.push_back(segMap(j));
	  }
	} 
	// copy all other pixels into objdata
	else {
	  objdata(i) = data(j);
	}
	objSegMap(i) = segMap(j);
	if (weight.size()!=0) 
	  objWeightMap(i) = weight(j);
      }
    }

    gsl_rng_free (r);

    // Grid will be changed but not shifted (all pixels stay at their position)
    O.accessGrid() = objSegMap.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);

    // Fill other quantities into Object
    O.setDetectionFlags(flags);
    O.history << "# Segment:" << endl;
    O.computeFluxCentroid();

    // Update catalog with object values
    (*catiter).second.XMIN = xmin;
    (*catiter).second.XMAX = xmax;
    (*catiter).second.YMIN = ymin;
    (*catiter).second.YMAX = ymax;
    (*catiter).second.XCENTROID = O.getCentroid()(0);
    (*catiter).second.YCENTROID = O.getCentroid()(1);
    (*catiter).second.FLUX = O.getFlux();
    (*catiter).second.FLAGS = (unsigned char) flags.to_ulong();
    (*catiter).second.CLASSIFIER = 0;
  } 
  // this is the whole frame
  else if (O.getID()==0) {
    O.history.clear();
    O.history << "# Extracting Object 0 (whole Fits image)." << endl;
    O.accessData() = Image<data_t>::getData();
    O.accessGrid() = Image<data_t>::getGrid();
    O.accessSegmentationMap() = segMap;
    if (weight.size()!=0)
      O.accessWeightMap() = weight;
    O.computeFluxCentroid();
  } else {
    std::cerr << "# Frame: This Object does not exist!" << endl;
    terminate();
  }
  O.setNoiseMeanRMS(noise_mean,noise_rms);
  O.setBaseFilename(Image<data_t>::getFilename());
}

const SegmentationMap& Frame::getSegmentationMap() {
  return segMap;
}

// now extend to region around the object by
// typically objectsize/4, minimum 8 pixels
void Frame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
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

const History& Frame::getHistory () {
  return history;
}

const std::list<unsigned int>& Frame::getPixelList(unsigned long objectnr) {
  std::map< unsigned long, std::list<unsigned int> >::iterator iter = 
    objectsPixels.find(objectnr);
  if (iter != objectsPixels.end())
    return (*iter).second;
  else {
    std::cerr << "# Frame: This Object does not exist!" << std::endl;
    terminate();  
  }
}

const Catalog& Frame::getCatalog() {
  return catalog;
}
