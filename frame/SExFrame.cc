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

SExFrame::SExFrame (std::string fitsfile) : FitsImage(fitsfile) {
  text << "# Reading FITS file " << fitsfile << endl;
  text << "# Image properties: size = "<< FitsImage::getSize(0) << "/" << FitsImage::getSize(1) << std::endl; 
  history.append(text);
  SExCatFormat empty = {0,0,0,0,0,0,0,0};
  sf = empty;
  catChecked = catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of unserlying FitsImage copied since often used
  axsize0 = FitsImage::getSize(0);
  axsize1 = FitsImage::getSize(1);
  // estimate background noise properties
  estimateNoise();

  // for artificial noise
  const gsl_rng_type * T;
  //gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  correlationLength = 4;
}

SExFrame::~SExFrame() {
  gsl_rng_free (r);
}

void SExFrame::readCatalog(std::string catfile) {
  // first inser empty object 0, since SExtractor starts with NUMBER 1
  objectList.clear();
  SExCatObject s0 = {0,0,0,0,0,0,0,0};
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
      so.BACKGROUND = (double) atof(column[sf.BACKGROUND-1].c_str());
      so.FLAGS = (unsigned int) atoi(column[sf.FLAGS-1].c_str());
      so.CLASS_STAR = (double) atof(column[sf.CLASS_STAR-1].c_str());
     // then push it on objectList
      objectList.push_back(so);
    }
  }
  catalog.close();
  catRead = 1;
}

void SExFrame::readSegmentationMap(std::string segmentfile) {
  FitsImage* seg = new FitsImage(segmentfile);
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

void SExFrame::setCharacteristicSize(uint size) {
  correlationLength = size;
}

void SExFrame::fillObject(Object& O) {
  unsigned int nr = O.getID();
  int xmin, xmax, ymin, ymax;
  xmin = objectList[nr].XMIN_IMAGE;
  xmax = objectList[nr].XMAX_IMAGE;
  ymin = objectList[nr].YMIN_IMAGE;
  ymax = objectList[nr].YMAX_IMAGE;

  O.history = history;
  text << "# Extracting Object " << O.getID() << ", ";
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
  
  // find the halo around the object and discard other fluctuations in the noise
  // therefore include quite a large border around
  addFrameBorder(1.5, xmin,xmax,ymin,ymax);
  findHalo(nr,xmin,xmax,ymin,ymax);
  text << "# Halo found: Enlarging object to (" << xmin << "/" << ymin << ") to (";
  text << xmax << "/" << ymax << ")" << std::endl;
  O.history.append(text);
  int innerxmin = xmin, innerxmax = xmax, innerymin = ymin, innerymax = ymax;

  addFrameBorder(0.25, xmin,xmax,ymin,ymax);
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
  // int the end only object data into new vector of smaller size, the rest will
  // filled up with artificial noise
  NumVector<double>& objdata = O.accessData();
  objdata = NumVector<double>((xmax-xmin+1)*(ymax-ymin+1));
  NumVector<double>& obj_bgmean = O.accessBackground();
  obj_bgmean = NumVector<double>((xmax-xmin+1)*(ymax-ymin+1));
  NumVector<double>& obj_bgrms = O.accessBackgroundRMS();
  obj_bgrms = NumVector<double>((xmax-xmin+1)*(ymax-ymin+1));
  
  makeLocalNoiseMap(nr,xmin,xmax,ymin,ymax,innerxmin,innerxmax,innerymin,innerymax,obj_bgmean,obj_bgrms);
  const NumVector<double>& data = FitsImage::getData();

  for (int i =0; i < objdata.size(); i++) {
    // old coordinates derived from new pixel index i
    int axis0 = xmax-xmin+1;
    int x = i%axis0 + xmin;
    int y = i/axis0 + ymin;
    uint j = FitsImage::getPixel(x,y);

    // if pixel is out of image region, fill noise from default values
    // since we fill same noise into data and into bgrms
    // the overall chi^2 (normalized by bg_rms) is unaffected by this region
    if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1) {
      objdata(i) = gsl_ran_gaussian (r, bg_rms);
      if (!subtractBG)
	objdata(i) += bg_mean;
      obj_bgmean(i) = bg_mean; 
      obj_bgrms(i) = bg_rms;
    } 
    //now inside image region
    else {
      // mask other objects in the frame and negative fluctuation in the inner region of 
      // the object (outside they are downweighted by their high variance)
      if ((segMap(j) != nr && segMap(j) != nr - 0.5 && segMap(j) > 0) ||
	  (x>innerxmin && x<innerxmax && y>innerymin && y<innerymax && segMap(j) == -3)){
 	objdata(i) = gsl_ran_gaussian (r, bg_rms);
 	if (!subtractBG)
  	  objdata(i) += bg_mean;
	obj_bgmean(i) = bg_mean; 
	obj_bgrms(i) = bg_rms;
      }
      else {
	objdata(i) = data(j);
	// bg mean not yet subtracted in subtractBackground()
	if (subtractBG) 
	  objdata(i) -= obj_bgmean(i);
      }
    }
  }
    
  // Grid will be changed but not shifted (all pixels stay at their position)
  O.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);
  
  // Fill other quantities into Object
  // take defaults for the noise values
  // if maps are given, the will be used instead anyway
  O.setNoiseMeanRMS(bg_mean,bg_rms);
  O.setNoiseModel("POISSONIAN");

  O.setDetectionFlag(objectList[nr].FLAGS);
  O.setStarGalaxyProbability(objectList[nr].CLASS_STAR);
  O.setBlendingProbability(computeBlendingProbability(nr));
  O.setBaseFilename(FitsImage::getFilename());
  // this calculates flux and centroid;
  O.getFlux();
}

// Searches for positive and negative oscillation largen than threshold (pixels)
// For each of the positive pixels it determines the minimum distance to the object 
// and fills distance into histogram. The halo is defined by all pixels with distances
// smaller than the distance to the first minimum in the histogram.
// The coordinates are updated to include the halo, the halo pixels are set to nr - 0.5 
// in the segmentation map, but can be overwritten in subsequent calls to this method.
void SExFrame::findHalo(uint nr, int& xmin, int& xmax, int& ymin, int& ymax) {
  const NumVector<double>& data = FitsImage::getData();
  // clean this area of the segmap for already detected positive fluctuations
  // and halo pixels of other objects
  cleanSegMapArea(nr,xmin,xmax,ymin,ymax);

  // set the minimal size of features to 1/3*correlationLength^2
  uint thresh = (uint)ceil(0.33*gsl_pow_2(correlationLength));
  std::list<uint> pixellist;
  std::vector<std::list<uint> > groups;
  list<uint>::iterator iter;
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      uint j = FitsImage::getPixel(x,y);
      // if pixel is noise: check if its starting point of a larger noise oscillation 
      // if this is true, segMap(j) = -1 (detected, but too small), 
      // -2 (positive), -3 (negative)
      if (segMap(j) == 0) {
	bool positive = 0;
	// positive
	if (data(j) > bg_mean + bg_rms) {
	  refinePixelMap(pixellist, j,1,xmin,xmax,ymin,ymax, bg_mean, bg_rms);
	  positive = 1;
	}
	// negative
	else if (data(j) < bg_mean - bg_rms)
	  refinePixelMap(pixellist, j,0,xmin,xmax,ymin,ymax, bg_mean, bg_rms);
	// found no larger oscillation: restore pixelmap
	if (pixellist.size() <= thresh)
	  for(iter = pixellist.begin(); iter != pixellist.end(); iter++)
	    segMap(*iter) = -1;
	// found bigger positive oscillation
	else if (positive)
	  groups.push_back(pixellist);
      }
    }
  }
  // find pixellist of the central object and remove inner pixel
  // such that only the border pixels are left over
  std::list<uint> centrallist;
  findObjectList(centrallist, nr, xmin, xmax, ymin, ymax);
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
      segMap(*iter) = distance;
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
  // also reset xmin, xmax etc. to values from objectList
  // and enlarge when halo pixels move outside this region
  xmin = objectList[nr].XMIN_IMAGE;
  xmax = objectList[nr].XMAX_IMAGE;
  ymin = objectList[nr].YMIN_IMAGE;
  ymax = objectList[nr].YMAX_IMAGE;
  for (int i=0; i < groups.size(); i++) {
    for(iter = groups[i].begin(); iter != groups[i].end(); iter++) {
      int pixel = *iter;
      int x, y;
      FitsImage::getCoords(pixel,x,y);
      // now set all found pixels (whose distance was stored in the segmap)
      // to nr - 0.5 if their distance is smaller than lower
      if (segMap(pixel) > 0 && segMap(pixel) <= lower) {
	segMap(pixel) = nr - 0.5;
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
	segMap(pixel) = - 2;
    }
  }
}

void SExFrame::subtractBackground() {
  subtractBG = 1;
}

// fill in the column position into the SExCatFormat struct
void SExFrame::insertFormatField(std::string type, std::string columnnr) {
  unsigned int colnr = atoi(columnnr.c_str());
  if (type.compare("NUMBER")==0)
    sf.NUMBER = colnr;
  if (type.compare("XMIN_IMAGE")==0)
    sf.XMIN_IMAGE = colnr;
  if (type.compare("XMAX_IMAGE")==0)
    sf.XMAX_IMAGE = colnr;
  if (type.compare("YMIN_IMAGE")==0)
    sf.YMIN_IMAGE = colnr;
  if (type.compare("YMAX_IMAGE")==0)
    sf.YMAX_IMAGE = colnr;
  if (type.compare("BACKGROUND")==0)
    sf.BACKGROUND = colnr;
  if (type.compare("FLAGS")==0)
    sf.FLAGS = colnr;
  if (type.compare("CLASS_STAR")==0)
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
  if (sf.BACKGROUND == 0) {
    text << "SExFrame: BACKGROUND keyword not provided!" << std::endl;
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

// based on the segmentation map and the pixel data the likelihood of blending
// of this object with other object(s) is estimated.
// Therefore the flux in average overlapping pixels (pixels of this object with a neighbor
// from another object) is compare the the average total flux of this object.
// This ratio (bounded by 1) is returned.
double SExFrame::computeBlendingProbability(unsigned int nr) {
  unsigned short flag = objectList[nr].FLAGS;
  // if 2 is in FLAGS: object flaged as blended
  if ((flag >> 1)%2 == 1) {
    int xmin, xmax, ymin, ymax;
    xmin = objectList[nr].XMIN_IMAGE;
    xmax = objectList[nr].XMAX_IMAGE;
    ymin = objectList[nr].YMIN_IMAGE;
    ymax = objectList[nr].YMAX_IMAGE;
    int total = 0, overlap = 0;
    double totalFlux = 0, overlapFlux = 0;
    const NumVector<double>& data = FitsImage::getData();

    // loop over all pixels within object borders
    set<unsigned int> overlapset;
    for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
      for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
	unsigned int pixel = FitsImage::getPixel(x,y);
	// pixel comes from selected object
	if (segMap(pixel) == nr) {
	  total++;
	  totalFlux += data(pixel) - bg_mean;
	  //if this pixel has not already been seen in overlap
	  // look at all neighbors
	  if (overlapset.find(pixel) == overlapset.end()) {
	    for (unsigned int dir = 1; dir < 9 ; dir++) {
	      int neighbor = neighborpixel(pixel,dir,axsize0,axsize1);
	      if (neighbor != -1) {
		// this pixel is due to another object
		if (segMap(neighbor) > 0 && segMap(neighbor) != nr && segMap(neighbor) != nr - 0.5) {
		  // store it first
		  overlapset.insert(pixel);
		  overlapFlux += data(pixel) - bg_mean;
		  break;
		}
	      }
	    }
	  }
	}
      }
    }
    overlap = overlapset.size();
    // idea: blend = average brightness in overlapping pixels in units of
    // average brightness in all pixels.
    // advantage: independet of the shape of the overlapping region
    double blend;
    if (overlap!=0) 
      blend = (overlapFlux*total)/(totalFlux*overlap);
    else
      blend = 0;
    // not bigger than 1
    blend = GSL_MIN(blend,1);
    return blend;
  }
  else return 0;
}

// returns the pixel index of the neighbor pixels or -1 if the pixel is outside the image
// area. The directions (starting from 1) go from top to right to bottom to left.
// Direction == 0 is the pixel itself.
int SExFrame::neighborpixel(int pixel,unsigned int direction, int axsize0, int axsize1) {
  int x,y,index;
  FitsImage::getCoords(pixel,x,y);
  switch(direction) {
  case 0: 
    // the pixel itself
    index = pixel;
    break;
  case 1: 
    if (y<axsize1-1) index = (y+1)*axsize0 + x ;  // top
    else index = -1;
    break;
  case 2:
    if (y<axsize1-1 && x<axsize0-1) index = (y+1)*axsize0 + x + 1;  // top right
    else index = -1;
    break;
  case 3:
    if (x<axsize0-1) index = y*axsize0 + x + 1;  // right neighbour
    else index = -1;
    break;
  case 4: 
    if (y>0 && x<axsize0-1) index = (y-1)*axsize0 + x + 1;  // bottom right
    else index = -1;
    break;  
  case 5: 
    if (y>0) index = (y-1)*axsize0 + x;  // bottom
    else index = -1;
    break;
  case 6: 
    if (y>0 && x>0) index = (y-1)*axsize0 + x - 1;  // bottom left
    else index = -1;
    break;   
  case 7: 
    if (x>0) index = y*axsize0 + x - 1; // left
    else index = -1;
    break;
  case 8: 
    if (y<axsize1-1 && x>0) index = (y+1)*axsize0 + x - 1;  // top left
    else index = -1;
    break;  
  }
  return index;
}

// find list of pixels due to given object from segmentation map 
void SExFrame::findObjectList(std::list<uint>& pixellist, int object, int xmin, int xmax, int ymin, int ymax) {
  pixellist.clear();
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      uint j = FitsImage::getPixel(x,y);
      if (segMap(j) == object)
	pixellist.push_back(j);
    }
  }
}

// removes all inner pixels of this object. Inner pixel are those with only neighbors 
// from the same object.
void SExFrame::removeObjectInnerPixels(std::list<uint>&  pixellist, uint nr) {
  std::list<uint> borderpixels;
  for(std::list<uint>::iterator iter = pixellist.begin(); iter != pixellist.end(); iter++) {
    int pixel = *iter;
    for (uint dir=0; dir<9; dir++) {
      int neighbor = neighborpixel(pixel,dir,axsize0,axsize1);
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

// this refines the pixelmap within the bounds of the smaller object frame
// it searches for noise oscillations above the 1 sigma level
void SExFrame::refinePixelMap(std::list<uint>& pixellist, uint startpixel, bool positive, int xmin, int xmax, int ymin, int ymax, double bg_mean, double bg_rms) {
  const NumVector<double>& data = FitsImage::getData();
  pixellist.clear();
  pixellist.push_back(startpixel);
  list<uint>::iterator theIterator;
  theIterator = pixellist.begin();
  int pixelnumber = 0;
  int x,y;
  while (pixelnumber < pixellist.size()) {
    int pixel = *theIterator;
    FitsImage::getCoords(pixel,x,y);
    if (x > xmin && x < xmax && y > ymin && y < ymax) {
      for (unsigned int dir = 1; dir < 9 ; dir++) {
	int neighbor = neighborpixel(pixel,dir,axsize0,axsize1);
	int xneighbor = neighbor%axsize0;
	int yneighbor = neighbor/axsize0;
	if (xneighbor >= xmin && xneighbor <= xmax && yneighbor >= ymin && yneighbor <= ymax) {
	  if (segMap(neighbor) == 0&& positive && data(neighbor) > bg_mean+bg_rms) {
	    pixellist.push_back(neighbor);
	    segMap(neighbor) = -2;
	  }
	  else if (segMap(neighbor) == 0 && !positive && data(neighbor) < bg_mean-bg_rms) {
	    pixellist.push_back(neighbor);
	    segMap(neighbor) = -3;
	  }
	}
      }
    }
    theIterator++;
    pixelnumber++;
  }
}

// Computes the statistics of the pixel values in boxes of 5x5 pixels around each pixel
// in the outer area without the inner area.
// Thus, we have a ring of local (noise) measurements which is then extrapolated into the
// inner region by comuting the four low-order polynomials on the border of the inner region
// (top, bottom, left, right). For finding the values inside the inner region, say for a
// point (x,y), we use the average of the values of the top and bottom polynomial at x and
// of the left and the right polynomial at y, weighted with their distance to the pixel.
void SExFrame::makeLocalNoiseMap(uint nr, int outerxmin, int outerxmax, int outerymin, int outerymax, int innerxmin, int innerxmax, int innerymin, int innerymax, NumVector<double>& mean, NumVector<double>& rms) {
  int cellsx = outerxmax-outerxmin+1, cellsy = outerymax - outerymin+1;
  int cellx=0, celly =0;
  // first define a ring of mean and rms values around the object
  for (int y = outerymin; y <= outerymax; y++) {
    cellx=0;
    for (int x=outerxmin; x<=outerxmax; x++) {
      if (!(x>innerxmin && x < innerxmax && y>innerymin && y<innerymax)) {
	int cell = cellx+celly*cellsx;
	// are we inside the image region?
	if (x>=0 && x < axsize0 && y >= 0 && y < axsize1)
	  getSubAreaStats(nr,FitsImage::getPixel(x,y),correlationLength,mean(cell),rms(cell));
	else {
	  mean(cell) = bg_mean;
	  rms(cell) = bg_rms;
	}
      }
      cellx++;
    }
    celly++;
  }
  // now do a LS fit to the top, bottom, left and right inner borders
  int innersx = innerxmax - innerxmin - 1, innersy = innerymax - innerymin - 1;
  NumVector<double> meantop(innersx), rmstop(innersx), meanbottom(innersx), rmsbottom(innersx), meanleft(innersy), rmsleft(innersy), meanright(innersy), rmsright(innersy);
  for (int i=0; i < innersx; i++) {
    meantop(i) = mean(i+innerxmin-outerxmin + (innerymax+1-outerymin)*cellsx);
    meanbottom(i) = mean(i+innerxmin-outerxmin + (innerymin-1-outerymin)*cellsx);
    rmstop(i) = rms(i+innerxmin-outerxmin + (innerymax+1-outerymin)*cellsx);
    rmsbottom(i) = rms(i+innerxmin-outerxmin + (innerymin-1-outerymin)*cellsx);
  }
  for (int j=0; j < innersy; j++) {
    meanleft(j) = mean(innerxmin-outerxmin-1 + (j+innerymin-outerymin)*cellsx);
    meanright(j) = mean(innerxmax-outerxmin+1 + (j+innerymin-outerymin)*cellsx);
    rmsleft(j) = rms(innerxmin-outerxmin-1 + (j+innerymin-outerymin)*cellsx);
    rmsright(j) = rms(innerxmax-outerxmin+1 + (j+innerymin-outerymin)*cellsx);
  }
  // build SVD of polynomial matrix
  int order = GSL_MIN_INT(4,GSL_MIN_INT(innersx,innersy)/10);
  NumMatrix<double> Atb(innersx,order+1);
  NumMatrix<double> Alr(innersy,order+1);
  for (int i=0; i<GSL_MAX_INT(innersx,innersy); i++) {
    for (int j=0;j<order+1; j++) {
      if (i<innersx)
	Atb(i,j) = gsl_pow_int(i,j);
      if (i<innersy)
	Alr(i,j) = gsl_pow_int(i,j);
    }
  }

  NumMatrix<double> Utb, Ulr, Stb, Slr, Vttb, Vtlr;
  Atb.svd(Utb,Stb,Vttb);
  Alr.svd(Ulr,Slr,Vtlr);
  for (int i=0; i< GSL_MIN_INT(Stb.getRows(),Stb.getColumns()); i++)
    if (Stb(i,i) != 0)
      Stb(i,i) = 1./Stb(i,i);
  for (int i=0; i< GSL_MIN_INT(Slr.getRows(),Slr.getColumns()); i++)
    if (Slr(i,i) != 0)
      Slr(i,i) = 1./Slr(i,i);

  Stb = Atb*Vttb.transpose()*Stb.transpose()*Utb.transpose();
  Slr = Alr*Vtlr.transpose()*Slr.transpose()*Ulr.transpose();
  
  // replace border lines by their polynomial approximation
  meantop = Stb*meantop;
  meanbottom = Stb*meanbottom;
  rmstop = Stb*rmstop;
  rmsbottom = Stb*rmsbottom;
  meanleft = Slr*meanleft;
  meanright = Slr*meanright;
  rmsleft = Slr*rmsleft;
  rmsright = Slr*rmsright;
  
  // now build weighted means of these value in the inner region
  // weigting is done by the distance to the appropriate line
  for (int i=0; i < innersx; i++) {
    for (int j=0; j < innersy; j++) {
      int cell = i+innerxmin+1 -outerxmin+ (j+innerymin+1-outerymin)*cellsx;
      mean(cell) = 0.5*(meantop(i)*j/(innersy-1) + meanbottom(i)*(innersy-1-j)/(innersy-1) +
      			meanleft(j)*(innersx-1-i)/(innersx-1) + meanright(j)*i/(innersx-1));
      rms(cell) = 0.5*(rmstop(i)*j/(innersy-1) + rmsbottom(i)*(innersy-1-j)/(innersy-1) +
		       rmsleft(j)*(innersx-1-i)/(innersx-1) + rmsright(j)*i/(innersx-1));
    }
  }
}

// computes the statistics (mean and rms w.r.t. the mean within squares of size sidelength+1;
// the given startpixel is in the center of the box.
// If we are at the boundary of the image, standard noise is added.
// In the inner region, we neglect pixels that are due to the object.
void SExFrame::getSubAreaStats(uint nr, uint startpixel, int sidelength, double& mean, double& rms) {
  const NumVector<double>& data = FitsImage::getData();
  int npixels = (sidelength+1)*(sidelength+1);
  double subarea[npixels];
  int index = 0;
  int x0,y0;
  FitsImage::getCoords(startpixel,x0,y0);
  uint pixel;
  for (int j=-sidelength/2; j<=sidelength/2; j++) {
    for (int i=-sidelength/2; i<=sidelength/2; i++) {
      int x = x0 + i, y = y0 + j;
      pixel = FitsImage::getPixel(x,y);
      // are we inside the image region?
      if (x>=0 && x < axsize0 && y >= 0 && y < axsize1) {
	// discard object pixels
	if (segMap(pixel) != nr && segMap(pixel) != nr - 0.5) {
	  subarea[index] = data(pixel);
	  index++;
	}
      }
      else {
	subarea[index] = bg_mean + gsl_ran_gaussian (r, bg_rms);
	index++;
      }
    }
  }
  mean = gsl_stats_mean(subarea,1,index);
  rms = gsl_stats_sd_m(subarea,1,index,mean);
}

// draws a rectangular frame with the given limits in the segmenation map
void SExFrame::paintSegmentBorder(int xmin, int xmax, int ymin, int ymax) {
  // check if corners are within frame
  if (xmin<0) xmin=0;
  if (xmax>=axsize0) xmax=axsize0-1;
  if (ymin<0) ymin=0;
  if (ymax>=axsize1) xmax=axsize1-1;
  // low border
  for (int i=(xmin+ymin*axsize0); i<=(xmax+ymin*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = -5;
  // high border
  for (int i=(xmin+ymax*axsize0); i<=(xmax+ymax*axsize0); i++)
    if (segMap(i) == 0)
      segMap(i) = -5;
  // left border
  for (int i=(xmin+(ymin+1)*axsize0);i<=(xmin+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = -5;
  // right border
  for (int i=(xmax+(ymin+1)*axsize0);i<=(xmax+(ymax-1)*axsize0); i+=axsize0)
    if (segMap(i) == 0)
      segMap(i) = -5;
}

// cleans the segmentation map within the given region from already detected
// positive fluctuations and halo pixels.
// The reason is that finding the fluctuations is faster than searching in the set
// of already found ones; and that pixels in the halo of other objects will be needed
// for finding the halo of this object (because of the rise after the first minimum).
void SExFrame::cleanSegMapArea(int nr, int xmin, int xmax, int ymin, int ymax) {
  int runningnr;
  for (int y = GSL_MAX_INT(ymin,0); y <= GSL_MIN_INT(ymax,axsize1-1); y++) {
    for (int x = GSL_MAX_INT(xmin,0); x <= GSL_MIN_INT(xmax,axsize0-1); x++) {
      runningnr = nr;
      uint j = FitsImage::getPixel(x,y);
      // positive fluctuations: have to be cleaned to be refound again
      if (segMap(j) == -2)
	segMap(j) = 0;
      // now pixels from halos of already detected objects
      else {
	while (runningnr>1) {
	  runningnr--;
	  if (segMap(j) == runningnr-0.5) {
	    segMap(j) = 0;
	    break;
	  }
	}
      }
    }
  }
}

const NumVector<double>& SExFrame::getObjectMap() {
  return segMap;
}


// estimate noise by iterative sigma clipping
void SExFrame::estimateNoise() {
  // for GSL sorting functions we have to copy NumVector to double*
  int npixels = FitsImage::getNumberOfPixels();
  double* D = (double *) malloc(npixels*sizeof(double));
  for (int i=0; i < npixels; i++)
    D[i] = (FitsImage::getData())(i);

  // first check left border of the image for pixel value variations
  // if its 0, its a (simulated) image with 0 noise
  // FIXME: background_variance is set to the minmal value above 0
  // Is this corect? 
  if (gsl_stats_variance(D,1,axsize0-1) == 0) {
    bg_mean = 0;
    double min = gsl_stats_max(D,1,npixels);
    for (int i=0; i < npixels; i++)
      if (D[i] < min && D[i] > 0) min = D[i];
    bg_rms = sqrt(min*min);
  } 
  else {
    gsl_sort(D, 1, npixels);
    double sigma = gsl_stats_sd(D,1,npixels);
    double median = gsl_stats_median_from_sorted_data (D,1,npixels);

    // sigma clipping here
    int j, jmax = npixels;
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

// computes for all pixels in objectlist the euclidean distance to given pixel
double SExFrame::distanceFromRim(std::list<uint>& objectlist, uint pixel) {
  int x0, y0;
  FitsImage::getCoords(pixel,x0,y0);
  int x,y;
  double mindist=axsize0;
  double distance;
  for(std::list<uint>::iterator iter = objectlist.begin(); iter != objectlist.end(); iter++) {
    FitsImage::getCoords(*iter,x,y);
    distance = sqrt(gsl_pow_2(x0-x) + gsl_pow_2(y0-y));
    if (distance < mindist)
      mindist = distance;
  }
  // the small offset if for discrimination with the 0 pixels
  return mindist + 0.01;
}

const History& SExFrame::getHistory () {
  return history;
}
