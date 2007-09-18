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

SExFrame::SExFrame (std::string fitsfile) : Image<data_t>(fitsfile) {
  SExCatFormat empty = {0,0,0,0,0,0,0,0,0,0};
  sf = empty;
  catChecked = catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of underlying Image copied since often used
  axsize0 = Image<data_t>::getSize(0);
  axsize1 = Image<data_t>::getSize(1);

  text << "# Reading FITS file " << fitsfile << endl;
  text << "# Image properties: size = "<< axsize0 << "/" << axsize1 << std::endl; 
  history.append(text);

  // for artificial noise
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
}

SExFrame::SExFrame (std::string fitsfile, std::string weightfile) : Image<data_t>(fitsfile) {
  weight = Image<data_t>(weightfile);
  SExCatFormat empty = {0,0,0,0,0,0,0,0,0,0};
  sf = empty;
  catChecked = catRead = segmapRead = subtractBG = estimatedBG = 0;
  bg_mean = bg_rms = 0;
  // axsizes of underlying Image copied since often used
  axsize0 = Image<data_t>::getSize(0);
  axsize1 = Image<data_t>::getSize(1);

  text << "# Reading FITS file " << fitsfile << endl;
  text << "# Image properties: size = "<< axsize0 << "/" << axsize1 << std::endl; 
  history.append(text);

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
      so.CLASS_STAR = (data_t) atof(column[sf.CLASS_STAR-1].c_str());
     // then push it on objectList
      objectList.push_back(so);
    }
  }
  catalog.close();
  catRead = 1;
}

void SExFrame::readSegmentationMap(std::string segmentfile) {
  Image<int>* seg = new Image<int>(segmentfile);
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
  unsigned int id = O.getID();
  int xmin, xmax, ymin, ymax;
  xmin = objectList[id].XMIN_IMAGE;
  xmax = objectList[id].XMAX_IMAGE;
  ymin = objectList[id].YMIN_IMAGE;
  ymax = objectList[id].YMAX_IMAGE;

  O.history = history;
  text << "# Extracting Object " << O.getID() << " (NUMBER = " <<  objectList[id].NUMBER << "), ";
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
  
  addFrameBorder(ShapeLensConfig::ADD_BORDER, xmin,xmax,ymin,ymax);
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
      if ((segMap(j) > 0 && segMap(j) != objectList[id].NUMBER) || (segMap(j) < 0 && ShapeLensConfig::FILTER_SPURIOUS)) {
 	// if we have a weight map 
	if (weight.size()!=0)
	  objdata(i) = gsl_ran_gaussian (r, sqrt(1./weight(j)));
	else
	  objdata(i) = gsl_ran_gaussian (r, bg_rms);
	// this objects has to yet been found to be nearby
	if (std::find(nearby_objects.begin(),nearby_objects.end(),segMap(j)) == nearby_objects.end()) {
	  text << "# Object " << segMap(j) << " nearby, but not overlapping." << std::endl;
	  O.history.append(text);
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
  if (weight.size()==0) {
    O.setNoiseMeanRMS(bg_mean,bg_rms);
    O.setNoiseModel("GAUSSIAN");
  }
  else {
    O.setNoiseMeanRMS(-1,-1);
    O.setNoiseModel("WEIGHT");
  }
  O.history.append("# Segment:\n");
  O.computeFluxCentroid();
  O.setFlux(objectList[id].FLUX_AUTO);
  Point2D centroid(objectList[id].XWIN_IMAGE,objectList[id].YWIN_IMAGE);
  O.setCentroid(centroid);

  O.setDetectionFlag(objectList[id].FLAGS);
  O.setStarGalaxyProbability(objectList[id].CLASS_STAR);
  O.setBlendingProbability(computeBlendingProbability(id));
  O.setBaseFilename(Image<data_t>::getFilename());
  O.setNumber(objectList[id].NUMBER);
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
  text << "# Background estimation (object masked):";
  text << " mean = " << bg_mean << ", sigma = " << bg_rms << endl;
  history.append(text);
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
  text << "# Noise estimates explicitly set: mean = " << mean << ", rms = " << rms << std::endl;
  history.append(text);
}

const History& SExFrame::getHistory () {
  return history;
}
