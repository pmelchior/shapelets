#ifndef SEXFRAME_H
#define SEXFRAME_H

#include <vector>
#include <list>
#include <gsl/gsl_rng.h>
#include <NumMatrix.h>
#include <Typedef.h>
#include <History.h>
#include <ShapeLensConfig.h>
#include <frame/Grid.h>
#include <frame/Image.h>
#include <frame/Object.h>
#include <frame/Catalog.h>

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
/// SExFrame* f = new SExFrame(fitsfile);
/// f->readCatalog(catfile);
/// f->readSegmentationMap(segfile);
/// f->subtractBackground();
/// for(int n = 1; n <= f->getNumberOfObjects(); n++) {
///   Object* obj = new Object(i);
///   f->fillObject(obj);
///   NumVector<data_t>& data = obj->getData();
///   Grid& grid = obj->getGrid();
///   ...
/// }
///\endcode

class SExFrame : public Image<data_t> {
 public:
  /// Argumented constructor.
  /// <tt>filename</tt> is the name of the Fits file.\n Extensions or other selections can be 
  // passed in the standard cfitsio way: <tt>filename[extension]</tt>.
  SExFrame(std::string filename);
  /// Argemented constructor for including a weight map image.
  /// <tt>data_file</tt> is the name of the Fits file containing the data, 
  /// <tt>weight_file</tt> the one of the corresponding weight map.\n
  /// Extensions or other selections can be 
  /// passed in the standard cfitsio way: <tt>data_file[extension]</tt>.
  SExFrame(std::string data_file, std::string weight_file);
  /// Destructor.
  ~SExFrame();
  /// Read the SExtractor catalog file.
  /// The format of this file has to be specified by setting the SExCatFormat properly.
  void readCatalog(std::string catfile);
  /// Read the segmentation map.
  /// This Fits file can be produced by SExtractor by setting 
  /// <tt>CHECKIMAGE_TYPE  SEGMENTATION</tt>\n
  /// The noise estimation step can be bypassed if the header keywords <tt>NOISE_MEAN</tt>
  /// and <tt>NOISE_RMS</tt> are set in this FITS file.
  void readSegmentationMap(std::string segmentfile);
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
  /// The Object(id) will then have a segmented version of the area around
  /// the SExtractor object specified by NUMBER.
  /// If subtractBackground() was called before, the background noise level
  /// will be subtracted.
  void fillObject(Object& obj);
  /// Return number of objects found by SExtractor.
  /// If it returns 0, the catalog file is either empty or its format is wrongly specified.
  unsigned int getNumberOfObjects();
  /// Get the map of found objects.
  /// The object map defines for each pixel, if it is part of an object.
  /// The convention is:
  /// - -1: non-significant pixel group, meaning above noise, but too small or faint
  /// - -2: positive small-scale noise oscillation bigger than 8 pixels close to detected object
  /// - -3: negative small-scale noise oscillation bigger than 8 pixels close to detected object
  /// - -4: borderline of the segments
  /// - 0: noise
  /// - 1..N: objectID of an identified significant pixel group
  const NumVector<int>& getObjectMap();
  /// Return Catalog read by readCatalog().
  const Catalog& getCatalog();

 private:
  void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
  data_t computeBlendingProbability(unsigned int objectNr);
  void estimateNoise();
  Catalog catalog;
  bool segmapRead, catRead, subtractBG, estimatedBG;
  NumVector<int> segMap;
  data_t bg_mean, bg_rms;
  unsigned int axsize0, axsize1;
  Image<data_t> weight;
  gsl_rng* r;
};



#endif
