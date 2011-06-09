#ifndef SHAPELENS_WCSTRANSFORMATION_H
#define SHAPELENS_WCSTRANSFORMATION_H

#include "CoordinateTransformation.h"
#include <fitsio.h>
#ifdef HAS_WCSLIB
  #include <wcslib/wcs.h>
#endif

namespace shapelens {
  /// Class for WCS coordinate transformations.
  /// This class depends on the presence of libwcs (and its headers).
  /// To make use of it, define \p -DHAS_WCSLIB in the environment 
  /// variable \p SPECIALFLAGS. 
  class WCSTransformation : public CoordinateTransformation {
  public:
    /// Constructor.
    WCSTransformation(fitsfile* fptr, bool intermediate = true);
    /// Copy constructor.
    WCSTransformation(const WCSTransformation& wcs);
    ///Destructor.
    ~WCSTransformation(); 
    /// Apply transformation to \p P.
    virtual void transform(Point<data_t>& P) const;
    /// Apply inverse transformation to \p P.
    virtual void inverse_transform(Point<data_t>& P) const;
    /// Get a copy of \p this.
    virtual boost::shared_ptr<CoordinateTransformation> clone() const;
  private:
    bool intermediate;
#ifdef HAS_WCSLIB  
    struct wcsprm wcs;
    double *imgcrd, *pixcrd, *world;
    int *stat;
#endif
  };
} // end namespace shapelens
 
#endif // SHAPELENS_WCSTRANSFORMATION_H

