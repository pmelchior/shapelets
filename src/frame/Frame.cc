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
#include <list>

using namespace shapelens;
using namespace std;

typedef unsigned int uint;
typedef unsigned long ulong;

Frame::Frame(string filename) : 
Image<data_t>(filename), segMap(), weight(), history(segMap.history) {
  history << "# Reading FITS file " << filename << endl;
  history << "# Image properties: size = "<< Frame::getSize(0) << "/" << Frame::getSize(1) << endl;
  segMap.resize(Frame::getSize(0)*Frame::getSize(1));
  segMap.clear();
  segMap.grid = Frame::grid;
  subtractedBG = estimatedBG = 0;
  noise_rms = noise_mean = 0;
  catalog.clear();
}

Frame::Frame(string datafile, string weightfile) : 
Image<data_t>(datafile), weight(weightfile), segMap(), history(segMap.history) {
  history << "# Reading data from " << datafile << " and weights from " << weightfile << endl;
  history << "# Image properties: size = "<< Frame::getSize(0) << "/" << Frame::getSize(1) << endl; 
  segMap.resize(Frame::getSize(0)*Frame::getSize(1));
  segMap.clear();
  segMap.grid = Frame::grid;
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
  Frame::kappa_sigma_clip(noise_mean,noise_rms);
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
    segMap.clear();
  }  

  // set up pixellist and objectsPixels
  set<ulong> pixelset;
  set<ulong>::iterator iter;

  // reset catalog and create dummy CatObject
  CatObject co = {0,0,0,0,0,0,0,0,0};

  // look for all positive fluctuations with at least 1 pixel above
  // highThreshold and minPixels above significanceThreshold
  data_t highThreshold, max, max_threshold;
  uint blending = 0;
  for (int i =0; i < npixels; i++) {
    highThreshold = getThreshold(i,ShapeLensConfig::DETECT_THRESHOLD);
    if (data(i) > highThreshold && segMap(i) == 0) {
      counter++;
      linkPixels(pixelset,max, max_threshold, i);
      if (ShapeLensConfig::BLENDING)
	blending = detectBlending(pixelset, max, max_threshold);
      if (pixelset.size() >= ShapeLensConfig::MIN_PIXELS) {
	history << "# Object " << counter << " detected with " << pixelset.size() << " significant pixels at (" << i%(Frame::getSize(0)) << "/" << i/(Frame::getSize(0)) << ")"  << std::endl;
	for(iter = pixelset.begin(); iter != pixelset.end(); iter++ )
	  segMap(*iter) = counter;
	objectsPixels[counter]= pixelset;
	catalog[counter] = co;
	catalog[counter].FLAGS = 2*blending;
      }
      else {
	counter--;
	for(iter = pixelset.begin(); iter != pixelset.end(); iter++ )
	  segMap(*iter) = -1;
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
  // in principle this could also be done with the segMap but does not have advantages
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
	  // if (segMap(neighbor == 0) {
	  if (pixelset.find(neighbor) == pixelset.end()) {
	    pixellist.push_back(neighbor);
	    pixelset.insert(neighbor);
	    //segMap(neighbor) = tag;
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

bool Frame::detectBlending(const set<ulong>& all, data_t max, data_t max_threshold) {
  // compute N thresholds between min_threshold an max
  uint N = ShapeLensConfig::BLEND_NTHRESH;
  NumVector<data_t> threshold(N);
  for (uint i=1; i <= N; i++)
    threshold(i-1) = i*(max - max_threshold)/(N+1);
  
  tree<set<ulong> > Tree;
  tree<set<ulong> >::iterator top;
  tree<set<ulong> >::fixed_depth_iterator diter;
  tree<set<ulong> >::sibling_iterator sib;
  
  // initialize tree with all
  top = Tree.begin();
  top = Tree.insert(top, all);
  // search for children of all
  insertNodesAboveThreshold(Tree,top, threshold(0));

  // go through all depth levels
  for (int i=1; i < N; i++) {
    if (Tree.max_depth(top) >= i) {
      diter = Tree.begin_fixed(top, i-1);
      for (uint j=0; j < diter.number_of_children(); j++) {
	sib = Tree.child(diter,j);
	insertNodesAboveThreshold(Tree, sib, threshold(i));
      }
    }
  }
  
  // now go thru all nodes from top to leafs
  // for nodes with more than a single child, check if
  // - the number of pixels are larger than ShapeLensConfig::MIN_PIXELS
  // - the flux of thos pixels is larger than a certain fraction of the total flux
  // for at least two children
  // otherwise remove all children apart from the biggest/brightest
  bool blending = 0;
  data_t allFlux = getTotalFlux(all);
  for(tree<set<ulong> >::breadth_first_iterator biter = Tree.begin_breadth_first(); biter != Tree.end_breadth_first(); biter++) {
    //sib = Tree.begin(biter);
    if (biter.number_of_children() > 1) {
      uint found = 0;
      for(sib = Tree.begin(biter); sib != Tree.end(biter); sib++)
	if ((*sib).size() > ShapeLensConfig::MIN_PIXELS)
	  if (getTotalFlux(*sib) > ShapeLensConfig::BLEND_MINCONT*allFlux)
	    found++;
      if (found >= 2)
	blending = 1;
    }
    if (blending)
      break;
  }
  return blending;
}

data_t Frame::getTotalFlux(const set<ulong>& pixels) {
  data_t flux = 0;
  const NumVector<data_t>& data = *this;
  for(set<ulong>::const_iterator iter = pixels.begin(); iter != pixels.end(); iter++)
    flux += data(*iter);
  return flux;
}

void Frame::insertNodesAboveThreshold(tree<set<ulong> >& Tree, tree<set<ulong> >::iterator_base& titer, data_t threshold) {
  set<ulong> above;
  set<ulong>& parent = *titer;
  const NumVector<data_t>& data = *this;
  // first find all pixels from parent that are above threshold
  for (set<ulong>::const_iterator iter = parent.begin(); iter != parent.end(); iter++)
    if (data(*iter) > threshold)
      above.insert(*iter);

  // iterate thru all above pixels, link the neighbors into first child
  // start a new child set when above pixel is not in a prior child
  set<ulong>::const_iterator aboveiter;
  tree<set<ulong> >::sibling_iterator sib;
  set<ulong> child;
  list<ulong> childlist;
  for (aboveiter = above.begin(); aboveiter != above.end(); aboveiter++) {
    ulong startpixel = *aboveiter;
    // check whether startpixel is already in one of the children
    bool found = 0;
    for (sib = Tree.begin(titer); sib != Tree.end(titer); sib++) {
      if ((*sib).find(startpixel) != (*sib).end()) {
	found = 1;
	break;
      }
   }
    // if it is not yet found, it belongs to a new child
    // thus, search all pixels of the new child
    if (!found) {
      child.clear();
      childlist.clear();
      childlist.push_back(startpixel);
      child.insert(startpixel);
      uint pixelnumber = 0;
      list<ulong>::iterator iter = childlist.begin();
      while (pixelnumber < childlist.size()) {
	ulong pixel = *iter;
	// loop over all direct neighbors
	for (uint dir = 1; dir <= 8 ; dir++) {
	  long neighbor = Frame::grid.getNeighborPixel(pixel,dir);
	  // neighbor is within the image
	  if (neighbor > 0)  {
	    // neighbor is not already in child
	    if (child.find(neighbor) == child.end()) {
	      // neighbor also in above set
	      if (above.find(neighbor) != above.end()) {
		childlist.push_back(neighbor);
		child.insert(neighbor);
	      }
	    }
	  }
	}
	iter++;
	pixelnumber++;	
      }
      // insert child into map of nodes and automatically increase key by using size()
      Tree.append_child(titer, child);
    }
  }
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
    for(set<ulong>::iterator iter = pixelset.begin(); iter != pixelset.end(); iter++ ) {
	uint x = (*iter)%axsize0;
	uint y = (*iter)/axsize0;
	if (x < xmin) xmin = x;
	if (y < ymin) ymin = y;
	if (x > xmax) xmax = x;
	if (y > ymax) ymax = y;
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
    O.segMap.history.setSilent();
    O.segMap.history = history;
    O.segMap.history.unsetSilent();
    O.segMap.resize((xmax-xmin)*(ymax-ymin));
    if (weight.size()!=0) 
      O.weight.resize((xmax-xmin)*(ymax-ymin));
    vector<uint> nearby_objects;

    // lop over all object pixels
    Point2D<int> P;
    for (int i =0; i < O.size(); i++) {
      // old coordinates derived from new pixel index i
      int axis0 = xmax-xmin;
      P(0) = i%axis0 + xmin;
      P(1) = i/axis0 + ymin;
      ulong j = Frame::grid.getPixel(P);

      // if pixel is out of image region, fill noise
      if (P(0) < 0 || P(1) < 0 || P(0) >= axsize0 || P(1) >= axsize1) {
	O(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	O.segMap(i) = 0;
	if (weight.size()!=0) 
	  O.weight(i) = 1./gsl_pow_2(noise_rms);
      } 
      else {
	// filter other objects in the frame
	if ((segMap(j) > 0 && segMap(j) != catiter->first) || (segMap(j) < 0 && ShapeLensConfig::FILTER_SPURIOUS)) {
	  // if we have a weight map 
	  if (weight.size()!=0)
	    O(i) = noise_mean + gsl_ran_gaussian (r, sqrt(1./weight(j)));
	  else
	    O(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	  flags[0] = 1;
	  // this objects has to yet been found to be nearby
	  if (std::find(nearby_objects.begin(),nearby_objects.end(),segMap(j)) == nearby_objects.end()) {
	    O.history << "# Object " << segMap(j) << " nearby, but not overlapping." << std::endl;
	    nearby_objects.push_back(segMap(j));
	  }
	} 
	// copy all other pixels into objdata
	else {
	  O(i) = data(j);
	}
	O.segMap(i) = segMap(j);
	if (weight.size()!=0) 
	  O.weight(i) = weight(j);
      }
    }

    gsl_rng_free (r);

    // Grid will be changed but not shifted (all pixels stay at their position)
    O.grid = O.weight.grid = O.segMap.grid = Grid(xmin,ymin,xmax-xmin,ymax-ymin);

    // Fill other quantities into Object
    O.flags = flags;
    O.computeFluxCentroid();

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
    write_iter->second.FLUX = O.flux;
    write_iter->second.FLAGS = (unsigned char) flags.to_ulong();
    write_iter->second.CLASSIFIER = 0;
  } 
  // this is the whole frame
  else if (O.id==0) {
    O.history.clear();
    O.history << "# Extracting Object 0 (whole Fits image)." << endl;
    O = *this;
    O.grid = Frame::grid;
    O.segMap = segMap;
    if (weight.size()!=0)
      O.weight = weight;
    O.computeFlux();
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
  return segMap;
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
// findObject() has been called -> segMap is meaningfull
  if(getNumberOfObjects())   
    return CorrelationFunction(*this,segMap,threshold);
  else
    return CorrelationFunction(*this,threshold);
}
