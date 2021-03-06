#include "../../include/frame/Frame.h"
#include "../../include/utils/MathHelper.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <algorithm>
#include <bitset>
#include <list>

using namespace shapelens;
using namespace std;

typedef unsigned int uint;
typedef unsigned long ulong;

Frame::Frame() : Image<data_t>(), segmentation(), weight(), history(segmentation.history) {
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

Frame::Frame(string filename) : 
Image<data_t>(filename), segmentation(), weight(), history(segmentation.history) {
  history << "# Image properties: size = "<< Frame::getSize(0) << "/" << Frame::getSize(1) << endl;
  segmentation.resize(Frame::getSize(0)*Frame::getSize(1));
  segmentation.clear();
  segmentation.grid = Frame::grid;
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

Frame::Frame(string datafile, string weightfile) : 
Image<data_t>(datafile), weight(weightfile), segmentation(), history(segmentation.history) {
  history << "# Reading data from " << datafile << " and weights from " << weightfile << endl;
  history << "# Image properties: size = "<< Frame::getSize(0) << "/" << Frame::getSize(1) << endl; 
  segmentation.resize(Frame::getSize(0)*Frame::getSize(1));
  segmentation.clear();
  segmentation.grid = Frame::grid;
  if (weight.size() != (*this).size()) {
    history << "Frame: weight map has different layout than data!" << endl;
    cerr << "Frame: weight map has different layout than data!" << endl;
    terminate();
  } 
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

// estimate noise by iterative sigma clipping
void Frame::estimateNoise() {
  std::pair<data_t, data_t> mean_std = kappa_sigma_clip(*this);
  noise_mean = mean_std.first;
  noise_rms = mean_std.second;
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
    for (int i=0; i < Frame::size(); i++) {
      Frame::operator()(i) -= noise_mean;
    }
    subtractedBG = 1;
    history << "# Background subtraction: noise level = " << noise_mean << std::endl;
    noise_mean = 0;
  }
}

// return threshold for findObject
// if weight map is provided, uses rms from this
data_t Frame::getThreshold(unsigned long pixel, data_t factor) {
  if (weight.size()==0) {
    if (!estimatedBG) estimateNoise();
    // for noise-free images
    if (noise_rms == 0)
      return noise_mean + factor;
    else
      return noise_mean + factor*noise_rms;
  }
  else
    return noise_mean + factor*sqrt(1./weight(pixel));
}

void Frame::findObjects() {
  const NumVector<data_t>& data = *this;
  unsigned long counter = 0;
  unsigned int npixels = Frame::size();
  
  // clean stuff from previous runs of this function
  if (getNumberOfObjects() > 0) {
    objectsPixels.clear();
    catalog.clear();
    segmentation.clear();
  }  

  // if segmap is not initialized (default constructor):
  // allocate the required space
  if (segmentation.size() == 0) {
    segmentation.resize(Frame::getSize(0)*Frame::getSize(1));
    segmentation.clear();
    segmentation.grid = Frame::grid;
  }

  // set up pixellist and objectsPixels
  set<ulong> pixelset;
  set<ulong>::iterator iter;

  // reset catalog and create dummy CatObject
  CatObject co;

  // look for all positive fluctuations with at least 1 pixel above
  // highThreshold and minPixels above significanceThreshold
  data_t highThreshold, max, max_threshold;
  uint blending = 0;
  for (int i =0; i < npixels; i++) {
    highThreshold = getThreshold(i,ShapeLensConfig::DETECT_THRESHOLD);
    if (data(i) > highThreshold && segmentation(i) == 0) {
      counter++;
      linkPixels(pixelset,max, max_threshold, i);
      if (pixelset.size() >= ShapeLensConfig::MIN_PIXELS) {
	history << "# Object " << counter << " detected with " << pixelset.size() << " significant pixels at (" << i%(Frame::getSize(0)) << "/" << i/(Frame::getSize(0)) << ")"  << std::endl;
	for(iter = pixelset.begin(); iter != pixelset.end(); iter++ )
	  segmentation(*iter) = counter;
	objectsPixels[counter]= pixelset;
	catalog[counter] = co;
	catalog[counter].FLAGS = 2*blending;
      }
      else {
	counter--;
	for(iter = pixelset.begin(); iter != pixelset.end(); iter++ )
	  segmentation(*iter) = -1;
      }
    }
  }
}

// Find a list of connected pixels which have values above threshold.
// This is a Friend-of-Friend algorithm with a linking length of 1 pixel.
// It starts by putting startpixel into the pixelset if data > noise_rms * threshold and if
// it was not in the set yet
// Then it performs these steps for all neighbors of pixels in the set until no
// new pixels are found.
void Frame::linkPixels(std::set<unsigned long>& pixelset, data_t& max, data_t& max_threshold, unsigned long startpixel) {
  pixelset.clear();
  pixelset.insert(startpixel);
  // the list is necessary for the algorithm below (using an iterator)
  // the set will be used to check whether a pixels is already there
  list<ulong> pixellist;
  pixellist.push_back(startpixel);

  list<ulong>::iterator iter = pixellist.begin();
  uint pixelnumber = 0;
  max = noise_mean;
  max_threshold = noise_mean;
  const NumVector<data_t>& data = *this;
  while (pixelnumber < pixellist.size()) {
    ulong pixel = *iter;
    // loop over all direct neighbors and check if they are above/below threshold
    for (uint dir = 1; dir <= 8 ; dir++) {
      long neighbor = Frame::grid.getNeighborPixel(pixel,dir);
      if (neighbor != -1) {
	data_t threshold_neighbor = getThreshold(neighbor, ShapeLensConfig::MIN_THRESHOLD);
	if (data(neighbor) > threshold_neighbor) {
	  // if (segmentation(neighbor == 0) {
	  if (pixelset.find(neighbor) == pixelset.end()) {
	    pixellist.push_back(neighbor);
	    pixelset.insert(neighbor);
	    //segmentation(neighbor) = tag;
	    if (data(neighbor) > max)
	      max = data(neighbor);
	    if (threshold_neighbor > max_threshold)
	      max_threshold = threshold_neighbor;
	  }
	}
      }
    }
    iter++;
    pixelnumber++;
  }
}

data_t Frame::getTotalFlux(const set<ulong>& pixels) {
  data_t flux = 0;
  const NumVector<data_t>& data = *this;
  for(set<ulong>::const_iterator iter = pixels.begin(); iter != pixels.end(); iter++)
    flux += data(*iter);
  return flux;
}

unsigned long Frame::getNumberOfObjects() {
  return catalog.size();
}

// cut the image to a small region around the object
// and set all pixels to zero than are not related to the image
void Frame::fillObject(Object& O, Catalog::const_iterator& catiter) {
  if (catiter != catalog.end()) {
    O.id = catiter->first;
    // set up detection flags bitset, use existing ones as starting value
    std::bitset<8> flags(catiter->second.FLAGS);

    // for artificial noise for area contaminated by different object
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // define the points (xmin/ymin) and (xmax/ymax) 
    // that enclose the object
    set<ulong>& pixelset = objectsPixels[catiter->first];
    int axsize0, axsize1,xmin, xmax, ymin, ymax;
    xmin = axsize0 = Frame::getSize(0);
    xmax = 0;
    ymin = axsize1 = Frame::getSize(1);
    ymax = 0;
    // loop over all pixels of current object
    Point<int> P;
    for(set<ulong>::iterator iter = pixelset.begin(); iter != pixelset.end(); iter++ ) {
      P = Frame::grid.getCoords(*iter);
      if (P(0) < xmin) xmin = P(0);
      if (P(0) > xmax) xmax = P(0);
      if (P(1) < ymin) ymin = P(1);
      if (P(1) > ymax) ymax = P(1);
    }

    O.history << "# Extracting Object " << catiter->first;
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
    const NumVector<data_t>& data = *this;
    O.resize((xmax-xmin)*(ymax-ymin));
    // Grid will be changed but not shifted (all pixels stay at their position)
    O.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
    O.grid.setWCS(Image<data_t>::grid.getWCS());
    O.segmentation.history.setSilent();
    O.segmentation.history = history;
    O.segmentation.history.unsetSilent();
    O.segmentation.resize((xmax-xmin)*(ymax-ymin));
    O.segmentation.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
    if (weight.size()!=0) {
      O.weight.resize((xmax-xmin)*(ymax-ymin));
      O.weight.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
    }
    vector<uint> nearby_objects;

    // lop over all object pixels
    for (int i =0; i < O.size(); i++) {
      P = O.grid.getCoords(i);
      long j = Frame::grid.getPixel(P);

      // if pixel is out of image region, fill noise
      if (j == -1) {
	O(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	O.segmentation(i) = 0;
	if (weight.size()!=0) 
	  O.weight(i) = 1./gsl_pow_2(noise_rms);
      } 
      else {
	// filter other objects in the frame
	if ((segmentation(j) > 0 && segmentation(j) != catiter->first) || (segmentation(j) < 0 && ShapeLensConfig::FILTER_SPURIOUS)) {
	  // if we have a weight map 
	  if (weight.size()!=0)
	    O(i) = noise_mean + gsl_ran_gaussian (r, sqrt(1./weight(j)));
	  else
	    O(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	  flags[0] = 1;
	  // this objects has to yet been found to be nearby
	  if (std::find(nearby_objects.begin(),nearby_objects.end(),segmentation(j)) == nearby_objects.end()) {
	    O.history << "# Object " << segmentation(j) << " nearby, but not overlapping." << std::endl;
	    nearby_objects.push_back(segmentation(j));
	  }
	} 
	// copy all other pixels into objdata
	else {
	  O(i) = data(j);
	}
	O.segmentation(i) = segmentation(j);
	if (weight.size()!=0) 
	  O.weight(i) = weight(j);
      }
    }

    gsl_rng_free (r);

    // Fill other quantities into Object
    O.flags = flags;
    O.computeCentroid();

    // Update catalog with object values
    // therefore we need a iterator from the const_iterator
    Catalog::iterator write_iter(catalog.begin());
    advance(write_iter, distance<Catalog::const_iterator>(write_iter, catiter));

    write_iter->second.XMIN = xmin;
    write_iter->second.XMAX = xmax;
    write_iter->second.YMIN = ymin;
    write_iter->second.YMAX = ymax;
    write_iter->second.XCENTROID = O.centroid(0);
    write_iter->second.YCENTROID = O.centroid(1);
    write_iter->second.FLAGS = (unsigned char) flags.to_ulong();
  } 
  // this is the whole frame
  else if (O.id==0) {
    O.history.clear();
    O.history << "# Extracting Object 0 (whole Fits image)." << endl;
    O = *this;
    O.grid = Frame::grid;
    O.segmentation = segmentation;
    if (weight.size()!=0)
      O.weight = weight;
    O.computeCentroid();
  } else {
    std::cerr << "# Frame: This Object does not exist!" << endl;
    terminate();
  }
  if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
    O.noise_mean = noise_mean;
    O.noise_rms = noise_rms;
  }
  O.basefilename = getFilename();
}

const SegmentationMap& Frame::getSegmentationMap() {
  return segmentation;
}

// now extend to region around the object by
// typically objectsize/4, minimum 8 pixels
void Frame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
  if (factor > 0) {
    int xrange, yrange, xborder, yborder;
    xrange = xmax - xmin;
    yrange = ymax - ymin;
    if (xrange%2 == 1) {
      xmax++;
      xrange++;
    }
    if (yrange%2 == 1) {
      ymax++;
      yrange++;
    }
    // make the object frame square, because of const beta in both directions
    if (xrange < yrange) {
      yborder = GSL_MAX_INT((int)floor(yrange*factor), 6);
      xborder = yborder + (yrange - xrange)/2;
    } else {
      xborder = GSL_MAX_INT((int)floor(xrange*factor), 6);
      yborder = xborder + (xrange - yrange)/2;
    }
    xmin -= xborder;
    xmax += xborder;
    ymin -= yborder;
    ymax += yborder;
  }
}

const History& Frame::getHistory () {
  return history;
}

const std::set<unsigned long>& Frame::getPixelSet(unsigned long objectnr) {
  std::map< unsigned long, std::set<unsigned long> >::iterator iter = 
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

CorrelationFunction Frame::computeCorrelationFunction(data_t threshold) {
// findObject() has been called -> segmentation is meaningfull
  if(getNumberOfObjects())   
    return CorrelationFunction(*this,segmentation,threshold);
  else
    return CorrelationFunction(*this,threshold);
}
