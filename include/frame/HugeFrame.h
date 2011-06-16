#ifndef SHAPELENS_HUGEFRAME_H
#define SHAPELENS_HUGEFRAME_H

#include <vector>
#include <list>
#include <gsl/gsl_rng.h>
#include "../Typedef.h"
#include "../ShapeLensConfig.h"
#include "Grid.h"
#include "Image.h"
#include "Object.h"
#include "Catalog.h"

namespace shapelens {

  /// Frame class for huge images.
  /// Insted of keeping the entire image in memory this class will only 
  /// keep pointers to the image files.\n
  /// fillObject() then copies only the regions which are enclosed by the 
  /// corner points <tt>XMIN..XMAX,YMIN..YMAX</tt> as given in the catalog.
  class HugeFrame {
  public:
    /// Constructor with image file and catalog file.
    HugeFrame(std::string datafile, std::string catfile);
    /// Constructor with image, weight map, and catalog.
    HugeFrame(std::string datafile, std::string weightfile, std::string catfile);
    /// Destruktor.
    ~HugeFrame();
    /// Return noise mean \f$\mu_n\f$.
    data_t getNoiseMean();
    /// Return noise RMS \f$\sigma_n\f$.
    data_t getNoiseRMS();
    /// Set noise features explicitly.
    void setNoiseMeanRMS(data_t mean, data_t rms);
    /// Fill all necessary information into an Object.
    /// The object is selected by the iterator \p catiter which must must
    /// be constructed from the Catalog obtained by calling getCatalog().\n
    /// The image will be cut to a region around the selected object. 
    /// If another objects is inside this region, its pixel will be replaced 
    /// by noise according to the noise statistics derived in estimateBackground().\n
    /// If subtractBackground() was called before, the background noise level
    /// will be subtracted.
    void fillObject(Object& obj, Catalog::const_iterator& catiter);
    /// Return number of objects found by SExtractor.
    /// If it returns 0, the catalog file is either empty or its format is wrongly specified.
    unsigned long getNumberOfObjects();
    /// Return Catalog read by readCatalog().
    const Catalog& getCatalog();
    /// Grid.
    Grid grid;

  private:
    void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
    Catalog catalog;
    fitsfile *fptr, *fptr_w;
    data_t bg_mean, bg_rms;
    long axsize0, axsize1;
    std::string basefilename;
  };
} // end namespace
#endif
