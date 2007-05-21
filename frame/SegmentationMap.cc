#include <SegmentationMap.h>
#include <IO.h>
#include <set>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <IO.h>

using namespace std;

typedef unsigned int uint;

SegmentationMap::SegmentationMap(const Image<double>& image) : Image<int>(), data(const_cast<NumVector<double>& >(image.getData())) {
  Image<int>::accessData().resize_clear(image.size());
  Image<int>::accessGrid() = image.getGrid();
}

SegmentationMap::SegmentationMap(string segMapFile, const Image<double>& image) : Image<int>(segMapFile), data(const_cast<NumVector<double>& >(image.getData())) {
}

SegmentationMap::SegmentationMap(const SegmentationMap& segIn) : Image<int>(), data(segIn.data) {
  Image<int>::operator=(segIn);
}

void SegmentationMap::operator=(const SegmentationMap& segMap) {
  Image<int>::operator=(segMap);
  data = segMap.data;
}

unsigned int SegmentationMap::getNumberOfObjects() {
  const NumVector<int>& segMap = Image<int>::getData();
  set<uint> objects;
  for (uint i=0; i<segMap.size(); i++)
    // if pixel belongs to object
    if (segMap(i) != 0)
      // inserts object number in set if it is not present yet
      objects.insert(segMap(i));
  return objects.size();
}

void SegmentationMap::linkPixelsSetMap(std::list<unsigned int>& pixellist, unsigned int startpixel, int tag, double threshold, bool positive) {
  NumVector<int>& segMap = Image<int>::accessData();
  pixellist.clear();
  pixellist.push_back(startpixel);
  segMap(startpixel) = tag;
  list<uint>::iterator theIterator;
  theIterator = pixellist.begin();
  uint& theEnd = pixellist.back();
  uint pixelnumber = 0;
  while (pixelnumber < pixellist.size()) {
    uint pixel = *theIterator;
    // loop over all direct neighbors and check if they are above/below threshold
    for (uint dir = 1; dir <= 8 ; dir++) {
      uint neighbor = Image<int>::getNeighborPixel(pixel,dir);
      if (neighbor != -1) {
	if ((data(neighbor) > threshold && positive) || 
	    (data(neighbor) < threshold && !positive)) {
	  if (segMap(neighbor) == 0) {
	    pixellist.push_back(neighbor);
	    segMap(neighbor) = tag;
	  }
	}
      }
    }
    theIterator++;
    pixelnumber++;
  }
}

void SegmentationMap::findHalo(uint nr, list<uint>& corelist, int& xmin, int& xmax, int& ymin, int& ymax, double correlationLength, double bg_mean, double bg_rms) {
  NumVector<int>& segMap = Image<int>::accessData();
  uint axsize0 = Image<int>::getSize(0), axsize1 = Image<int>::getSize(1);

  // save coordinate ranges for later usage
  int xmin0 = xmin, xmax0 = xmax, ymin0 = ymin, ymax0 = ymax;
  // use a rather large area for the search of fluctuations
  addFrameBorder(1.5, xmin,xmax,ymin,ymax);

  // clean this area of the segmap for already detected positive fluctuations
  // and halo pixels of other objects
  //cleanSegMapArea(xmin,xmax,ymin,ymax);

  // minimum nummber of pixels to be regarded as "group"
  uint minPixels = (uint) ceil(0.5*gsl_pow_2(correlationLength));
  // positive and negative thresholds
  double thresh_pos = bg_mean + bg_rms, thresh_neg = bg_mean - bg_rms;

  std::list<uint> pixellist;
  std::vector<std::list<uint> > groups;
  list<uint>::iterator iter;
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      uint j = Image<int>::getPixel(x,y);
      // if pixel is noise: check if its starting point of a noise oscillation 
      // if this is true, segMap(pixel) =  -1 (positive), -2 (negative)
      if (segMap(j) == 0) {
	bool positive = 0;
	// positive
	if (data(j) > thresh_pos) {
	  positive = 1;
	  linkPixelsSetMap(pixellist,j, -1, thresh_pos, positive);
	}
	// negative
	else if (data(j) < thresh_neg)
	  linkPixelsSetMap(pixellist,j, -2, thresh_neg, positive);

	// found no larger oscillation: restore pixelmap
	if (pixellist.size() <= minPixels)
	  for(iter = pixellist.begin(); iter != pixellist.end(); iter++)
	    segMap(*iter) = 0;
	// found bigger positive oscillation: put it into groups
	else if (positive)
	  groups.push_back(pixellist);
      }
    }
  }

  // remove inner pixels from list of core pixels
  std::list<uint> centrallist = corelist;
  removeObjectInnerPixels(centrallist, nr);

  // build histogram of distances
  // maximal distances are from the center to the top corners
  // since the size of the features we want tor resolve is correlationLength
  // we will use maxdist/correlationLength bins
  double maxdist = (double) GSL_MAX_INT(xmax-xmin,ymax-ymin)/M_SQRT2;
  gsl_histogram * h = gsl_histogram_alloc ((unsigned int) ceil(maxdist/correlationLength));
  gsl_histogram_set_ranges_uniform (h, 0, maxdist);

  // for each pixel in pixellist in groups compute minimal distance to object
  // which is to the boundary of the object
  double distance;
  for (int i=0; i < groups.size(); i++) {
    for(iter = groups[i].begin(); iter != groups[i].end(); iter++) {
      distance = distanceFromRim(centrallist,*iter);
      // store the distance temporarily in segmap
      segMap(*iter) = (int) ceil(distance);
      gsl_histogram_increment(h,distance);
    }
  }

  // find first minimum of histogram after the first peak (which should be our halo)
  // first find the peak
  int i=1;
  while (gsl_histogram_get(h,i) > gsl_histogram_get(h,i-1))
    i++;
  // now find the coming minimum
  while (gsl_histogram_get(h,i) < gsl_histogram_get(h,i-1))
    i++;
  int minbin = i-1;
  double lower, upper;
  gsl_histogram_get_range (h,minbin,&lower,&upper);

  // for all pixels in groups with smaller distance than lower:
  // set their segmap entry to nr-0.5 (probably halo)
  // also reset xmin, xmax etc. to minimum values.
  // and enlarge when halo pixels move outside this region
  xmin = xmin0;
  xmax = xmax0;
  ymin = ymin0;
  ymax = ymax0;
  for (int i=0; i < groups.size(); i++) {
    for(iter = groups[i].begin(); iter != groups[i].end(); iter++) {
      uint pixel = *iter;
      uint x, y;
      Image<int>::getCoords(pixel,x,y);
      // now set all found pixels (whose distance was stored in the segmap)
      // to nr if their distance is smaller than lower
      // also add them to coreList
      if (segMap(pixel) > 0 && segMap(pixel) <= lower) {
	segMap(pixel) = nr;
	corelist.push_back(pixel);
	// is new pixel outside the region of the object?
	// if yes set new xmin, xmax ...
	if (x < xmin)
	  xmin = x;
	else if (x > xmax)
	  xmax = x;
	if (y < ymin)
	  ymin = y;
	else if (y > ymax)
	  ymax = y;
      }
      // the other pixels: too far away to be halo
      else 
	segMap(pixel) = - 1;
    }
  }
}

// now extend to region around the object by
// typically objectsize/2, minimum 12 pixels
void SegmentationMap::addFrameBorder(double factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
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

// find list of pixels due to given object from segmentation map 
void SegmentationMap::findObjectPixels(std::list<uint>& pixellist, uint objectnr, int xmin, int xmax, int ymin, int ymax) {
  const NumVector<int>& segMap = Image<int>::getData();
  uint axsize0 = Image<int>::getSize(0), axsize1 = Image<int>::getSize(1);
  pixellist.clear();
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      uint j = Image<int>::getPixel(x,y);
      if (segMap(j) == objectnr)
	pixellist.push_back(j);
    }
  }
}

// removes all inner pixels of this object. Inner pixel are those with only neighbors 
// from the same object.
void SegmentationMap::removeObjectInnerPixels(std::list<uint>&  pixellist, uint nr) {
  const NumVector<int>& segMap = Image<int>::getData();
  std::list<uint> borderpixels;
  for(std::list<uint>::iterator iter = pixellist.begin(); iter != pixellist.end(); iter++) {
    int pixel = *iter;
    for (uint dir=0; dir<9; dir++) {
      int neighbor = Image<int>::getNeighborPixel(pixel,dir);
      if (neighbor != -1) {
	if (segMap(neighbor) != nr) {
	  borderpixels.push_back(neighbor);
	  break;
	}
      }
    }
  }
  pixellist = borderpixels;
}

// draws a rectangular frame with the given limits in the segmenation map
void SegmentationMap::setSegmentBorder(int tag, int xmin, int xmax, int ymin, int ymax) {
  NumVector<int>& segMap = Image<int>::accessData();
  uint axsize0 = Image<int>::getSize(0), axsize1 = Image<int>::getSize(1);
  // check if corners are within frame
  if (xmin<0) xmin=0;
  if (xmax>=axsize0) xmax=axsize0-1;
  if (ymin<0) ymin=0;
  if (ymax>=axsize1) xmax=axsize1-1;
  // low border
  for (int i=(xmin+ymin*axsize0); i<=(xmax+ymin*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // high border
  for (int i=(xmin+ymax*axsize0); i<=(xmax+ymax*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // left border
  for (int i=(xmin+(ymin+1)*axsize0);i<=(xmin+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // right border
  for (int i=(xmax+(ymin+1)*axsize0);i<=(xmax+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = tag;
}

// cleans the segmentation map within the given region from already detected
// positive fluctuations and halo pixels.
// The reason is that finding the fluctuations is faster than searching in the set
// of already found ones; and that pixels in the halo of other objects will be needed
// for finding the halo of this object (because of the rise after the first minimum).
void SegmentationMap::cleanSegMapArea(int xmin, int xmax, int ymin, int ymax) {
  NumVector<int>& segMap = Image<int>::accessData();
  uint axsize0 = Image<int>::getSize(0), axsize1 = Image<int>::getSize(1);
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      uint j = Image<int>::getPixel(x,y);
      // positive fluctuations: have to be cleaned to be refound again
      if (segMap(j) == -1)
	segMap(j) = 0;
    }
  }
}

// computes for all pixels in objectlist the euclidean distance to given pixel
double SegmentationMap::distanceFromRim(std::list<uint>& objectlist, uint pixel) {
  uint x0, y0;
  Image<int>::getCoords(pixel,x0,y0);
  uint x,y;
  double mindist=INFINITY;
  double distance;
  for(std::list<uint>::iterator iter = objectlist.begin(); iter != objectlist.end(); iter++) {
    Image<int>::getCoords(*iter,x,y);
    distance = sqrt(gsl_pow_2(x0-x) + gsl_pow_2(y0-y));
    if (distance < mindist)
      mindist = distance;
  }
  // the small offset if for discrimination with the 0 pixels
  return mindist + 0.01;
}


void SegmentationMap::save(std::string filename, std::map<std::string, std::string> keywords = std::map<std::string, std::string>()) {
  writeFITSFile(filename,Image<int>::getGrid(),Image<int>::getData(),keywords);
}


