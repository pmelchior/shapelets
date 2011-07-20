#ifndef SHAPELENS_SEXFRAME_H
#define SHAPELENS_SEXFRAME_H

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
  /// Wrapper class for SExtractor.
  /// Provides segmentation of a FitsFile into various Object entities employing 
  /// SExtractor's output files and some additional processing.\n
  /// It is assumed that for each input FitsFile there's a SExtractor catalog file
  /// (in <tt>ASCII_HEAD</tt> format) and a segmentation map (which is by SExtractor 
  /// generated if <tt>CHECKIMAGE_TYPE SEGMENTATION</tt> is specified).\n
  /// The required catalog
  /// parameters are 
  /// - <tt>NUMBER,</tt>
  /// - <tt>XMIN_IMAGE, XMAX_IMAGE, YMIN_IMAGE, YMAX_IMAGE,</tt>
  /// - <tt>XWIN_IMAGE,YWIN_IMAGE,</tt>
  /// - <tt>FLUX_AUTO,</tt>
  /// - <tt>FLAGS</tt> and
  /// - <tt>CLASS_STAR</tt>;
  ///
  /// the ordering is not important.\n
  /// The segmentation map is expected to have <tt>TINT</tt> FITS datatype (which is not provided by
  /// SExtractor). To convert (and compress) it one can use the <tt>fitscopy</tt> tool from 
  /// the FITS toolbox:
  /// \code
  /// fitscopy 'input.fits[pixi X]' '!output.fits'
  /// gzip 'output.fits'
  /// \endcode
  ///
  /// The noise characteristics are estimated globally by \f$\sigma\f$-clipping.
  ///
  /// The Object entities will have these features:
  /// - The frames is quadratic such that the area within the object was located 
  /// by SExtractor (<tt>XMIN..XMAX,YMIN..YMAX</tt>) plus additional border area is included.
  /// - Overlapping objects are masked with noise.
  /// - If subtractBackground() has been called before, the global background mean is
  /// subtracted from the pixel data.
  /// - If a weight map is not given, <tt>noiseModel="GAUSSIAN"</tt>, otherwise 
  ///   <tt>noisemodel="WEIGHT"</tt>
  /// - \p detectionFlag and \p StarGalaxyProbability of Object are filled with SExtractor's
  /// \p FLAGS and \p CLASS_STAR, respectively.
  /// - \p BlendingProbability \f$b\f$ is set to 1 if SExtractor's flag contains 2,
  /// otherwise it is set to 0.
  ///
  /// Example:
  /// \code
  /// SExFrame* f = new SExFrame(fitsfile, catfile);
  /// f->subtractBackground();
  /// const Catalog& cat = f->getCatalog();
  /// Catalog::const_iterator iter;
  /// for(iter = cat.begin(); iter != cat.end(); iter++) {
  ///   unsigned long id = (*iter).first;
  ///   Object obj;
  ///   f->fillObject(obj,iter);
  ///   if ((*iter).second.FLAGS == 0) {
  ///     // do something useful with obj ...
  ///   }
  /// }
  ///\endcode
  ///
  ///
  /// This Fits file can be produced by SExtractor by setting 
  /// <tt>CHECKIMAGE_TYPE  SEGMENTATION</tt>\n
  /// The noise statistics are read from the header keywords <tt>NOISE_MEAN</tt>
  /// and <tt>NOISE_RMS</tt> in either the data file or the segmentation file.

  class SExFrame {
  public:
    /// Argumented constructor for including a weight map image.
    /// \p datafile, \p weightfile and \p segmapfile are the names of the data,
    /// weight map and segmenation map FITS files 
    /// (Extensions or other selections can be passed in the standard cfitsio way:
    /// <tt>filename[extension]</tt>).\n
    /// \p catfile is the name of the SExtractor catalog file.
    /// \p weightfile and \p segmapfile can be empty string <tt>""</tt> then
    /// they will be ignored.
    SExFrame(std::string datafile, std::string catfile,  std::string segmapfile = "", std::string weightfile = "");
    /// Destructor.
    ~SExFrame();
    /// Return noise mean \f$\mu_n\f$.
    data_t getNoiseMean();
    /// Return noise RMS \f$\sigma_n\f$.
    data_t getNoiseRMS();
    /// Set noise features explicitly.
    void setNoiseMeanRMS(data_t mean, data_t rms);
    /// Subtract constant background noise level.
    /// The noise level is obtained from the BACKGROUND keyword in the catalog file 
    /// for object independently.
    void subtractBackground();
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
    fitsfile *fptr, *fptr_w, *fptr_s;
    History history;
    std::string basefilename;
    void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
    Catalog catalog;
    bool subtractBG;
    data_t bg_mean, bg_rms;
    long axsize0, axsize1;
    gsl_rng* r;
  };
} // end namespace
#endif
