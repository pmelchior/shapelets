#include "../../include/frame/SExFrame.h"
#include "../../include/utils/MathHelper.h"
#include <gsl/gsl_randist.h>
#include <set>
#include <vector>
#include <bitset>

using namespace shapelens;
using namespace std;

typedef unsigned int uint;

SExFrame::SExFrame (std::string datafile, std::string catfile, std::string segmapfile, std::string weightfile) : 
  catalog(catfile), basefilename(datafile) {
  subtractBG = false;
  bg_mean = bg_rms = 0;

  // open datafile
  fptr = IO::openFITSFile(datafile);
  long naxis[2];
  int status = 0;
  int dimensions;
  fits_get_img_dim (fptr, &dimensions, &status);
  if (dimensions == 2) {
    fits_get_img_size(fptr, dimensions, naxis, &status);
    axsize0 = naxis[0];
    axsize1 = naxis[1];
    grid.setSize(0,0,axsize0, axsize1);
    if (ShapeLensConfig::USE_WCS) {
#ifdef HAS_WCSLIB
      grid.setWCS(WCSTransformation(fptr));
#else
      throw std::runtime_error("IO: WCS usage requested, but HAS_WCSLIB not specified");
#endif
    }
  } else
    throw std::invalid_argument("SExFrame: FITS file does not provide valid image!");
  
  // open weight and segmentation maps if specified
  if (weightfile != "")
    fptr_w = IO::openFITSFile(weightfile);
  else
    fptr_w = NULL;
  if (segmapfile != "")
    fptr_s = IO::openFITSFile(segmapfile);
  else
    fptr_s = NULL;

  // for artificial noise
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  // check if NOISE_MEAN and NOISE_RMS is given as header keyword 
  // in datafile or segmapfile
  try {
    IO::readFITSKeyword(fptr,"BG_MEAN",bg_mean);
    IO::readFITSKeyword(fptr,"BG_RMS",bg_rms);
  } catch (std::exception) {
    if (fptr_s != NULL) {
      try {
	IO::readFITSKeyword(fptr_s,"BG_MEAN",bg_mean);
	IO::readFITSKeyword(fptr_s,"BG_RMS",bg_rms);
      } catch (std::exception) {}
    }
  }
}

SExFrame::~SExFrame() {
  gsl_rng_free (r);
  IO::closeFITSFile(fptr);
  IO::closeFITSFile(fptr_w);
  IO::closeFITSFile(fptr_s);
}

unsigned long SExFrame::getNumberOfObjects() {
  return catalog.size();
}

const Catalog& SExFrame::getCatalog() {
  return catalog;
}

void SExFrame::fillObject(Object& O, Catalog::const_iterator& catiter) {
  if (catiter != catalog.end()) {
    O.id = catiter->first;
    int xmin, xmax, ymin, ymax;
    xmin = catiter->second.XMIN;
    xmax = catiter->second.XMAX;
    ymin = catiter->second.YMIN;
    ymax = catiter->second.YMAX;

    O.history << "# Extracting Object " << (*catiter).first;
    O.history << " found in the area (" << xmin << "/" << ymin << ") to (";
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

    // prepare containers
    O.resize((xmax-xmin)*(ymax-ymin));
    // Grid will be changed but not shifted (all pixels stay at their position)
    O.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
    O.centroid = Point<data_t>(catiter->second.XCENTROID,
			       catiter->second.YCENTROID);
    if (ShapeLensConfig::USE_WCS)
      O.grid.setWCS(grid.getWCS());

    if (fptr_w != NULL) {
      O.weight.resize((xmax-xmin)*(ymax-ymin));
      O.weight.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
      if (ShapeLensConfig::USE_WCS)
      O.weight.grid.setWCS(grid.getWCS());
    }
    if (fptr_s != NULL) {
      O.segMap.resize((xmax-xmin)*(ymax-ymin));
      O.segMap.grid.setSize(xmin,ymin,xmax-xmin,ymax-ymin);
      if (ShapeLensConfig::USE_WCS)
      O.segMap.grid.setWCS(grid.getWCS());
    }

    // copy pixel data
    int status = 0;
    data_t nullval = 0;
    int anynull = 0;
    long firstpix[2] = {xmin+1,ymin+1}, lastpix[2] = {xmax, ymax}, inc[2] = {1,1};
    fits_read_subset(fptr, IO::getFITSDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.c_array(), &anynull, &status);
    if (fptr_w != NULL) 
      fits_read_subset(fptr_w, IO::getFITSDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.weight.c_array(), &anynull, &status);
    if (fptr_s != NULL) 
      fits_read_subset(fptr_s, IO::getFITSDataType(long(0)), firstpix, lastpix, inc, &nullval, O.segMap.c_array(), &anynull, &status);

    // check image pixels: replace pixels outside the original frame
    // or those belonging to another object by background noise
    if (subtractBG)
      O -= bg_mean;
    if (ShapeLensConfig::CHECK_OBJECT && (anynull != 0 || fptr_s != NULL)) {
      vector<uint> nearby_objects;
      for (int i =0; i < O.size(); i++) {
	bool fill = false;
	Point<int> P = O.grid.getCoords(i);
	// outside
	if (P(0) < 0 || P(0) >= axsize0 || P(1) < 0 || P(1) >= axsize1)
	  fill = true;

	// check segmap
	if (fptr_s != NULL) {
	  if ((O.segMap(i) > 0 && O.segMap(i) != catiter->first) || (O.segMap(i) < 0 && ShapeLensConfig::FILTER_SPURIOUS))
	    fill = true;
	  // this objects has to yet been found to be nearby
	  if (std::find(nearby_objects.begin(),nearby_objects.end(),O.segMap(i)) == nearby_objects.end()) {
	    O.history << "# Object " << O.segMap(i) << " nearby, but not overlapping." << std::endl;
	    nearby_objects.push_back(O.segMap(i));
	  }
	}

	if (fill) {
	  if (fptr_w != NULL) {
	    if (O.weight(i) != 0)
	      O(i) = gsl_ran_gaussian (r, 1./sqrt(O.weight(i)));
	  }
	  else
	    O(i) = gsl_ran_gaussian (r, bg_rms);
	  if (!subtractBG)
	    O(i) += bg_mean;
	}
      }
    }

    // Fill other quantities into Object
    O.flags = std::bitset<8>(catiter->second.FLAGS);
    O.basefilename = basefilename;
    if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
      O.noise_rms = bg_rms;
      if (subtractBG)
	O.noise_mean = 0;
      else
	O.noise_mean = bg_mean;
    }
  }
  else {
    std::ostringstream mess;
    mess << "SExFrame: Object " << O.id << " does not exist!";
    throw std::invalid_argument(mess.str());
  }
}

void SExFrame::subtractBackground() {
  subtractBG = 1;
}

// now extend to region around the object by
// typically objectsize/2, minimum 16 pixels
void SExFrame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
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
      yborder = std::max((int)floor(yrange*factor), 16);
      xborder = yborder + (yrange - xrange)/2;
    } else {
      xborder = std::max((int)floor(xrange*factor), 16);
      yborder = xborder + (xrange - yrange)/2;
    }
    xmin -= xborder;
    xmax += xborder;
    ymin -= yborder;
    ymax += yborder;
  }
}

data_t SExFrame::getNoiseMean() {
  return bg_mean;
}

data_t SExFrame::getNoiseRMS() {
  return bg_rms;
}

void SExFrame::setNoiseMeanRMS(data_t mean, data_t rms) {
  bg_mean = mean;
  bg_rms = rms;
}
