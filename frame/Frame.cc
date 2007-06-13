#include <Frame.h>
#include <NumVectorMasked.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
using namespace std;

typedef unsigned int uint;

Frame::Frame(string filename) : Image<double>(filename), segMap(*this) {
  text << "# Reading FITS file " << filename << endl;
  text << "# Image properties: size = "<< Image<double>::getSize(0) << "/" << Image<double>::getSize(1) << endl; 
  history.append(text);
  subtractedBG = estimatedBG = foundObjects= 0;
  noise_rms = noise_mean = 0;
  numberofObjects = 0;
}

// estimate noise by iterative sigma clipping
void Frame::estimateNoise() {
  Image<double>::getData().kappa_sigma_clip(noise_mean,noise_rms);
  text << "# Background estimation: mean = " << noise_mean;
  text << ", sigma = " << noise_rms << std::endl;
  history.append(text);
  estimatedBG = 1;
}

double Frame::getNoiseMean() {
  if (!estimatedBG) estimateNoise();
  return noise_mean;
}

double Frame::getNoiseRMS() {
  if (!estimatedBG) estimateNoise();
  return noise_rms;
}

void Frame::setNoiseMeanRMS(double mean, double rms) {
  estimatedBG = 1;
  noise_mean = mean;
  noise_rms = rms;
  text << "# Noise estimates explicitly set: mean = " << mean << ", rms = " << rms << std::endl;
  history.append(text);
}
  
void Frame::subtractBackground() {
  if (!estimatedBG) estimateNoise();
  if (!subtractedBG) {
    for (int i=0; i < Image<double>::size(); i++) {
      (Image<double>::accessData())(i) -= noise_mean;
    }
    subtractedBG = 1;
    text << "# Background subtraction: noise level = " << noise_mean << std::endl;
    history.append(text);
    noise_mean = 0;
  }
}

void Frame::findObjects(unsigned int minPixels, double significanceThresholdSigma, double detectionThresholdSigma) {
  if (!foundObjects) {
    const NumVector<double>& data = Image<double>::getData();
    int counter = 1;
    unsigned int npixels = Image<double>::size();

    // set up pixellist and objectsPixels
    list<uint> pixellist;
    list<uint>::iterator iter;
    objectsPixels.push_back(pixellist); // since object numbers start with 1, add a empty list

    // define threshold with global noise values and supplied threshold levels
    if (!estimatedBG) estimateNoise();
    double lowThreshold = noise_mean + significanceThresholdSigma*noise_rms;
    double highThreshold = noise_mean + detectionThresholdSigma*noise_rms;

    // find the maximum pixel in the field
    double maxvalue=0;
    uint maxindex = 0;
    for (uint i=0; i< npixels; i++) {
      if (data(i) > maxvalue) {
	maxvalue = data(i);
	maxindex = i;
      }
    }
    // there is a detected maximum (not a flat (noisy) image)
    if (data(maxindex) > highThreshold) {
      // search for all pixels connected to the maximum one
      segMap.linkPixelsSetMap (pixellist,maxindex,counter,lowThreshold,1);
      // if maximum has enough pixels define it as object 1
      if (pixellist.size() >= minPixels) {
	text << "# Maximum Object 1 detected with " << pixellist.size() << " significant pixels at (" << maxindex%(Image<double>::getSize(0)) << "/" << maxindex/(Image<double>::getSize(0)) << ")" << std::endl;
	history.append(text);
	objectsPixels.push_back(pixellist);
      // if not, mark the pixels int the active list with -1 (seen, but not big enough)
      } else {
	for(iter = pixellist.begin(); iter != pixellist.end(); iter++ )
	  segMap(*iter) = 0;
	counter--;
      }

      // now look for all other objects apart from the brightest one
      for (int i =0; i < npixels; i++) {
	if (data(i) > highThreshold && segMap(i) == 0) {
	  counter++;
	  segMap.linkPixelsSetMap(pixellist,i,counter,lowThreshold,1);
	  if (pixellist.size() >= minPixels) {
	    text << "# Object " << counter << " detected with " << pixellist.size() << " significant pixels at (" << i%(Image<double>::getSize(0)) << "/" << i/(Image<double>::getSize(0)) << ")"  << std::endl;
	    history.append(text);
	    objectsPixels.push_back(pixellist);
	  }
	  else {
	    for(iter = pixellist.begin(); iter != pixellist.end(); iter++ )
	      segMap(*iter) = 0;
	    counter--;
	  }
	}
      }
    } else {
      counter = 0;
    }

    numberofObjects = counter;
    if (numberofObjects >= 1)
      foundObjects = 1;
  }
}

unsigned int Frame::getNumberOfObjects() {
  return numberofObjects;
}

// cut the image to a small region around the object
// and set all pixels to zero than are not related to the image
void Frame::fillObject(Object& O) {
  if (O.getID() > 0 && O.getID() <= numberofObjects) {

    // for artificial noise for area contaminated by different object
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // define the points (xmin/ymin) and (xmax/ymax) 
    // that enclose the object
    list<uint>& pixellist = objectsPixels[O.getID()];
    int axsize0, axsize1,xmin, xmax, ymin, ymax;
    xmin = axsize0 = Image<double>::getSize(0);
    xmax = 0;
    ymin = axsize1 = Image<double>::getSize(1);
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

    O.history = history;
    text << "# Extracting Object " << O.getID() << ", ";
    text << "found in the area (" << xmin << "/" << ymin << ") to (";
    text << xmax << "/" << ymax << ")" << endl;
    O.history.append(text);
    
    // check if outer sizes of the object are identical to the image
    // boundary, since then the objects is cutted 
    if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
      O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),4));
      O.history.append("# Object cut off at the image boundary!\n");
    }

    // add border around object for including the edges
    addFrameBorder(xmin,xmax,ymin,ymax);
    text << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
    text << xmax << "/" << ymax << ")" << endl;
    O.history.append(text);

    // check if object was close to the image boundary so that noise has to be added
    if (xmin < 0 || ymin < 0 || xmax >= axsize0 || ymax >= axsize1) {
      O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),2));
      O.history.append("# Object close to image boundary: Possible cut-off. Extending frame with noise.\n");
    }

    // set noise estimates
    // since the position of the object is known now, we can compute the
    // noise now on all pixels not in pixellist
    // 1) set mask(i)=1, when i is pixellist
    // 2) create NumVectorMasked from objdata and mask
    // 3) compute std from that
    const NumVector<double>& data = Image<double>::getData();
    NumVector<bool> mask(data.size());
    for(list<uint>::iterator iter = pixellist.begin(); iter != pixellist.end(); iter++ )
      mask(*iter) = 1;
    NumVectorMasked<double> masked(data,mask);
    double masked_mean, masked_rms;
    masked.kappa_sigma_clip(masked_mean,masked_rms);
    text << "# Background estimation (object masked):";
    text << " mean = " << masked_mean << ", sigma = " << masked_rms << endl;
    history.append(text);
    O.setNoiseMeanRMS(masked_mean,masked_rms);

    // fill the object pixel data
    NumVector<double>& objdata = O.accessData();
    objdata.resize((xmax-xmin+1)*(ymax-ymin+1));
    SegmentationMap& objSegMap = O.accessSegmentationMap();
    objSegMap.resize((xmax-xmin+1)*(ymax-ymin+1));

    for (int i =0; i < objdata.size(); i++) {
      // old coordinates derived from new pixel index i
      int axis0 = xmax-xmin+1;
      int x = i%axis0 + xmin;
      int y = i/axis0 + ymin;
      uint j = Image<double>::getGrid().getPixel(x,y);

      // if pixel is out of image region, fill noise
      if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1) {
	objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	objSegMap(i) = 0;
      } 
      else {
	// filter other objects in the frame
	if (segMap(j) > 0 && segMap(j) != O.getID()) {
	  objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	  O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),1));
	  O.history.append("# Another object nearby, but not overlapping.\n");
	} 
	// copy all other pixels into objdata
	else {
	  objdata(i) = data(j);
	}
	objSegMap(i) = segMap(j);
      }
    }

    gsl_rng_free (r);

    // Grid will be changed but not shifted (all pixels stay at their position)
    O.accessGrid() = objSegMap.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);

    // Fill other quantities into Object
    O.setNoiseModel("GAUSSIAN");
    O.setBaseFilename(Image<double>::getFilename());
    // this calculates flux and centroid;
    O.history.append("# Segment:\n");
    O.computeFluxCentroid();

  } 
  // this is the whole frame
  else if (O.getID()==0) {
    O.history = history;
    O.history.append("# Extracting Object 0 (whole Fits image).\n");
    O.accessData() = Image<double>::getData();
    O.accessGrid() = Image<double>::getGrid();
    O.accessSegmentationMap() = segMap;
    O.setNoiseModel("GAUSSIAN");
    O.setNoiseMeanRMS(noise_mean,noise_rms);
    O.setBaseFilename(Image<double>::getFilename());
    O.computeFluxCentroid();
  } else {
    text.str("# Frame: This Object does not exist!\n");
    O.history.append(text);
    terminate();
  }
}

SegmentationMap& Frame::getSegmentationMap() {
  return segMap;
}

// now extend to region around the object by
// typically objectsize/4, minimum 8 pixels
void Frame::addFrameBorder(int& xmin, int& xmax, int& ymin, int& ymax) { 
  int xrange, yrange, xborder, yborder;
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  switch(xrange%4) {
  case 1: xmax--; break;
  case 2: xmin++; xmax++; break;
  case 3: xmax++; break;
  }
  switch(yrange%4) {
  case 1: ymax--; break;
  case 2: ymin++; ymax++; break;
  case 3: ymax++; break;
  }
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  // make the object frame square, because of const beta in both directions
  if (xrange < yrange) {
    yborder = GSL_MAX_INT(yrange/4, 10);
    xborder = yborder + (yrange - xrange)/2;
  } else {
    xborder = GSL_MAX_INT(xrange/4, 10);
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

const std::list<unsigned int>& Frame::getPixelList(unsigned int objectnr) {
  if (objectnr>0 && objectnr <= numberofObjects)
    return objectsPixels[objectnr];
  else {
    std::cout << "# Frame: This Object does not exist!" << std::endl;
    terminate();  
  }
}
std::list<unsigned int>& Frame::accessPixelList(unsigned int objectnr) {
  if (objectnr>0 && objectnr <= numberofObjects)
    return objectsPixels[objectnr];
  else {
    std::cout << "# Frame: This Object does not exist!" << std::endl;
    terminate();  
  }
}
