#include "../../include/frame/HugeFrame.h"
#include "../../include/ShapeLensConfig.h"
#include <numla/NumVectorMasked.h>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <set>
#include <vector>
#include <bitset>

using namespace shapelens;
using namespace std;
using namespace boost;

typedef unsigned int uint;

HugeFrame::HugeFrame (std::string datafile, std::string catfile) :
  catalog(catfile) {
  bg_mean = bg_rms = subtractBG = 0;

  fptr = IO::openFITSFile(datafile);
  fptr_w = NULL;

  // axsizes of underlying image (otherwise unknown)
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
    throw std::invalid_argument("HugeFrame: FITS file does not provide valid image!");
}

HugeFrame::HugeFrame (std::string datafile, std::string weightfile, std::string catfile) : 
  catalog(catfile) {
  bg_mean = bg_rms = subtractBG = 0;

  fptr = IO::openFITSFile(datafile);
  fptr_w = IO::openFITSFile(weightfile);

  // axsizes of underlying image (otherwise unknown)
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
      throw std::runtime_error("IO: WCS usage requested, but HAS_WCSToolsLib not specified");
#endif
    }
  } else
    throw std::invalid_argument("HugeFrame: FITS file does not provide valid image!");
}

HugeFrame::~HugeFrame() {
  IO::closeFITSFile(fptr);
  if (fptr_w != NULL) 
    IO::closeFITSFile(fptr_w);
}

unsigned long HugeFrame::getNumberOfObjects() {
  return catalog.size();
}

const Catalog& HugeFrame::getCatalog() {
  return catalog;
}

void HugeFrame::fillObject(Object& O, Catalog::const_iterator& catiter) {
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
  
    // add several pixel to original box
    addFrameBorder(ShapeLensConfig::ADD_BORDER,xmin,xmax,ymin,ymax);
    
    // check if outer sizes of the object are identical to the image
    // boundary, since then the objects is cutted 
    bool cutflag = 0;
    if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
      O.history << "# Object cut off at the image boundary!" << endl;
      cutflag = 1;
    }
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
    }

    // copy pixel data
    int status = 0;
    data_t nullval = 0;
    int anynull = 0;
    long firstpix[2] = {xmin+1,ymin+1}, lastpix[2] = {xmax, ymax}, inc[2] = {1,1};
    fits_read_subset(fptr, IO::getFITSDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.c_array(), &anynull, &status);
    if (subtractBG)
      O -= bg_mean;

    if (fptr_w != NULL) 
      fits_read_subset(fptr_w, IO::getFITSDataType(data_t(0)), firstpix, lastpix, inc, &nullval, O.weight.c_array(), &anynull, &status);

    // Fill other quantities into Object
    O.flags = std::bitset<8>(catiter->second.FLAGS);
    O.basefilename = basefilename;
    O.noise_rms = bg_rms;
    O.noise_mean = bg_mean;
  }
  else {
    std::ostringstream mess;
    mess << "HugeFrame: Object " << O.id << " does not exist!";
    throw std::invalid_argument(mess.str());
  }
}

data_t HugeFrame::getNoiseMean() {
  return bg_mean;
}

data_t HugeFrame::getNoiseRMS() {
  return bg_rms;
}

void HugeFrame::setNoiseMeanRMS(data_t mean, data_t rms) {
  bg_mean = mean;
  bg_rms = rms;
}

void HugeFrame::subtractBackground() {
  subtractBG = 1;
}

void HugeFrame::addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax) { 
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
      yborder = GSL_MAX_INT((int)floor(yrange*factor), 16);
      xborder = yborder + (yrange - xrange)/2;
    } else {
      xborder = GSL_MAX_INT((int)floor(xrange*factor), 16);
      yborder = xborder + (xrange - yrange)/2;
    }
    xmin -= xborder;
    xmax += xborder;
    ymin -= yborder;
    ymax += yborder;

    // ensure cutout is inside the whole frame
    xmin = std::max(0,xmin);
    ymin = std::max(0,ymin);
    xmax = std::min(axsize0-1,long(xmax));
    ymax = std::min(axsize1-1,long(ymax));

  }
}
