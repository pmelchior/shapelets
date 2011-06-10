#ifdef HAS_WCSLIB
#include "../../include/frame/WCSTransformation.h"
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>
#include <wcslib/wcs.h>

namespace shapelens {

  WCSTransformation::WCSTransformation(fitsfile* fptr, bool intermediate_) : 
    intermediate(intermediate_)  {
    int status = 0, nkeyrec, nreject, nwcs;
    char *header;
    
    // Read in the FITS header, excluding COMMENT and HISTORY keyrecords.
    if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
      status = 0;
      fits_file_name(fptr, header, &status);
      throw std::runtime_error("WCSTransformation: Cannot read header of file at " + std::string(header));
    }

    // Interpret the WCS keywords
    struct wcsprm* wcss; // use pointer to read multiple wcs structs if necessary
    status = wcspih(header, nkeyrec, WCSHDR_all, 0, &nreject, &nwcs, &wcss);
    free(header);
    if (status)
      throw std::runtime_error("WCSTransformation: Cannot read WCS header keywords (" + std::string(wcshdr_errmsg[status]) + ")");

    // check there is one (and one only) WCS with 2 coordinate axes
    if (wcss == NULL)
      throw std::runtime_error("WCSTransformation: No world coordinate systems found.");
    else if (nwcs > 1) {
      wcsvfree(&nwcs, &wcss);
      throw std::runtime_error("WCSTransformation: More than one world coordinate systems found.");
    }
    else if (wcss->naxis != 2) {
      wcsvfree(&nwcs, &wcss);
      throw std::runtime_error("WCSTransformation: WCS does not have 2 axes");
    }

    // initialize this wcs structure and copy it from first (and only)
    // entry of wcss
    wcs.flag = -1;           // wcslib implementation detail
    wcsini(1, 2, &wcs);      // 1: allocate memory for 2: axes
    wcscopy(0, wcss, &wcs);
    status = wcsset(&wcs);   // set remaining wcs fields from read information
    wcsvfree(&nwcs, &wcss);  // free the read-in structure

    if (status)
      throw std::runtime_error("WCSTransformation: wcsset error (" + std::string(wcs_errmsg[status]) + ")");

    // initialize coordinate containers
    world  = (double*) realloc(NULL,  2 * sizeof(double));
    imgcrd = (double*) realloc(NULL, 2 * sizeof(double));
    pixcrd = (double*) realloc(NULL, 2 * sizeof(double));
    stat   = (int*) realloc(NULL,   2 * sizeof(int));
  }
  
  // explicit definition since we have to allocate containers
  // and perform a deep copy of wcs
  WCSTransformation::WCSTransformation(const WCSTransformation& W) {
    world  = (double*) realloc(NULL,  2 * sizeof(double));
    imgcrd = (double*) realloc(NULL, 2 * sizeof(double));
    pixcrd = (double*) realloc(NULL, 2 * sizeof(double));
    stat   = (int*) realloc(NULL,   2 * sizeof(int));
    wcs.flag = -1;
    wcsini(1, 2, &wcs);
    wcscopy(0, &(W.wcs), &wcs);
    wcsset(&wcs);
    intermediate = W.intermediate;
  }
  
  // explicit definition to deallocate all structures
  WCSTransformation::~WCSTransformation() {
    int nwcs = 1;
    wcsfree(&wcs);
    free(world);
    free(imgcrd);
    free(pixcrd);
    free(stat);
  }

  void WCSTransformation::transform(Point<data_t>& P) const {
    *pixcrd = P(0);    
    *(pixcrd+1) = P(1);
    // use intermediate world coordinates
    // rather then celestial
    if (intermediate) {
      linp2x(const_cast<linprm*>(&(wcs.lin)), 1, 2, pixcrd, imgcrd);
      P(0) = *imgcrd;
      P(1) = *(imgcrd+1);
    } else {
      double phi, theta;
      wcsp2s(const_cast<wcsprm*>(&wcs), 1, 2, pixcrd, imgcrd, &phi, &theta, world, stat);
      P(0) = *world;
      P(1) = *(world+1);
    }
    stack_transform(P);
  }
  void WCSTransformation::inverse_transform(Point<data_t>& P) const {
    // inverse: this trafo comes last
    stack_inverse_transform(P);
    // use intermediate world coordinates (as input)
    // rather than celestial
    if (intermediate) {
      *imgcrd = P(0);    
      *(imgcrd+1) = P(1);
      linx2p(const_cast<linprm*>(&(wcs.lin)), 1, 2, imgcrd, pixcrd);
    }
    else {
      *world = P(0);    
      *(world+1) = P(1);
      double phi, theta;
      wcss2p (const_cast<wcsprm*>(&wcs), 1, 2, world, &phi, &theta, imgcrd, pixcrd, stat);
    }
    P(0) = (*pixcrd);
    P(1) = *(pixcrd+1);
  }

  boost::shared_ptr<CoordinateTransformation> WCSTransformation::clone() const {
    return boost::shared_ptr<CoordinateTransformation>(new WCSTransformation(*this));
  }
} // end namespace
#endif // HAS_WCSLIB
