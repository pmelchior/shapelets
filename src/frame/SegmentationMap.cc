#include <frame/SegmentationMap.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>

using namespace std;

typedef unsigned long ulong;

SegmentationMap::SegmentationMap() : 
Image<long>() {
}

SegmentationMap::SegmentationMap(string segMapFile) : 
Image<long>(segMapFile) {
}

unsigned long SegmentationMap::getNumberOfObjects() {
  const NumVector<long>& segMap = *this;
  set<long> objects;
  for (ulong i=0; i<segMap.size(); i++)
    // if pixel belongs to object
    if (segMap(i) != 0)
      // inserts object number in set if it is not present yet
      objects.insert(segMap(i));
  return objects.size();
}


// find list of pixels due to given object from segmentation map 
void SegmentationMap::findObjectPixels(std::set<ulong>& pixelset, ulong objectnr, long xmin, long xmax, long ymin, long ymax) {
  const NumVector<long>& segMap = *this;
  long axsize0 = Image<long>::grid.getSize(0), axsize1 = Image<long>::grid.getSize(1);
  pixelset.clear();
  for (long y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (long x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      ulong j = Image<long>::grid.getPixel(x,y);
      if (segMap(j) == objectnr)
	pixelset.insert(j);
    }
  }
}

// draws a rectangular frame with the given limits in the segmenation map
void SegmentationMap::setSegmentBorder(int tag, long xmin, long xmax, long ymin, long ymax) {
  NumVector<long>& segMap = *this;
  long axsize0 = Image<long>::grid.getSize(0), axsize1 = Image<long>::grid.getSize(1);
  // check if corners are within frame
  if (xmin<0) xmin=0;
  if (xmax>=axsize0) xmax=axsize0-1;
  if (ymin<0) ymin=0;
  if (ymax>=axsize1) xmax=axsize1-1;
  // low border
  for (long i=(xmin+ymin*axsize0); i<=(xmax+ymin*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // high border
  for (long i=(xmin+ymax*axsize0); i<=(xmax+ymax*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // left border
  for (long i=(xmin+(ymin+1)*axsize0);i<=(xmin+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = tag;
  // right border
  for (long i=(xmax+(ymin+1)*axsize0);i<=(xmax+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = tag;
}

void SegmentationMap::cleanSegMapArea(long xmin, long xmax, long ymin, long ymax) {
  NumVector<long>& segMap = *this;
  long axsize0 = SegmentationMap::getSize(0), axsize1 = SegmentationMap::getSize(1);
  for (long y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (long x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      ulong j = Image<long>::grid.getPixel(x,y);
      segMap(j) = 0;
    }
  }
}
