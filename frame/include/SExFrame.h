#ifndef SEXFRAME_H
#define SEXFRAME_H

#include <Grid.h>
#include <NumMatrix.h>
#include <FitsImage.h>
#include <Object.h>
#include <History.h>
#include <vector>
#include <list>
#include <gsl/gsl_rng.h>

/// Wrapper class for SExtractor.
/// Provides segmentation of an FitsFile into various Object entities employing 
/// SExtractor's output files and various other methods.\n
/// It is assumed that for each input FitsFile there's a SExtractor catalog file
/// (in <tt>ASCII_HEAD</tt> format) and a segmentation map (which is by SExtractor 
/// generated if <tt>CHECKIMAGE_TYPE SEGMENTATION</tt> is specified). The required catalog
/// parameters are <tt>NUMBER, XMIN_IMAGE, XMAX_IMAGE, YMIN_IMAGE, YMAX_IMAGE, BACKGROUND,
/// FLAGS, CLASS_STAR</tt>; the ordering is not important.
///
/// The noise characteristics are first estimated globally by \f$\sigma\f$-clipping.
/// Based on these values, in the regions of interest all sufficiently large fluctuations
/// (based on the FWHM of features in the image, set in setCharacteristicSize())
/// are identified. The positive fluctuations in the direct
/// neighborhood of the object are likely to be related to the object (halo) and thus
/// added to the definition of the object. The halo of the object (which is usually left 
/// out by SExtractor) is found using a
/// histogramming technique: For each pixel which is member of a positive fluctuation its
/// minimal distance to the object is computed and filled into a histogram. The halo should
/// create a high number of pixels with short distances. We therefore define the halo as
/// those pixels with smaller distances than the distance of the minimum after the first 
/// peak in the histogram.
/// \image html histogram_halo_small.png
///
/// In the region outside the object (incl. halo) the mean and variance w.r.t. the mean
/// is computed in boxes of sidelength given by setCharacteristicSize(). 
/// These maps are then extrapolated to the
/// region of the object. These maps provide the necessarily local information on the background, which are used for subtractBackground() and the definition of the weight map. 
///
/// The Object entities will have these features:
/// - The frames is quadratic such that the area within the object was located 
/// by SExtractor (<tt>XMIN..XMAX,YMIN..YMAX</tt>), extended by the halo, is completely 
/// included and additional border area is included.
/// - Overlapping objects are masked with noise.
/// - If subtractBackground() has been called before, the local background mean map is
/// subtracted from the pixel data.
/// - \p noiseModel is "POISSONIAN".
/// - \p detectionFlag and \p StarGalaxyProbability of Object are filled with SExtractor's
/// \p FLAGS and \p CLASS_STAR, respectively.
/// - \p BlendingProbability \f$b\f$ is obtained from SExtractor's segementation map: 
/// \f$ b= \min(\frac{F_{overlap}}{F_{total}} \frac{total}{overlap},1)\f$, where \f$total\f$
/// and \f$overlap\f$ are the number of total and overlapping pixels and \f$F\f$ is the flux
/// integrated within these pixels. A pixel is regarded as overlapping if it has a direct
/// neighbor pixel which is due to another object. With <tt>DEBLEND_MINCONT  0.010</tt> a
/// reasonable value for discarding a object because of blending is \f$b \approx 0.5..0.6\f$.
/// - When the background map and/or the background RMS map are provided, these maps are stored
/// in the region of the image segement.
///
/// Example:
/// \code
/// SExFrame* f = new SExFrame(fitsfile,extension);
/// f->readCatalog(catfile, sformat);
/// f->readSegmentationMap(segfile,0);
/// f->subtractBackground();
/// f->setCharacteristicSize(4);
/// for(int n = 1; n <= f->getNumberOfObjects(); n++) {
///   Object* obj = new Object(i);
///   f->fillObject(obj);
///   NumVector<double>& data = obj->getData();
///   Grid& grid = obj->getGrid();
///   ...
/// }
///\endcode

class SExFrame : public FitsImage {
 public:
  /// Argumented constructor.
  /// The name of the Fits file and to appropriate extension have to be given.
  SExFrame(std::string fitsfile, unsigned int extension);
  /// Destructor.
  ~SExFrame();
  /// Read the SExtractor catalog file.
  /// The format of this file has to be specified by setting the SExCatFormat properly.
  void readCatalog(std::string catfile);
  /// Read the segmentation map.
  /// This Fits file can be produced by SExtractor by setting 
  /// <tt>CHECKIMAGE_TYPE  SEGMENTATION</tt>
  void readSegmentationMap(std::string segmentfile, unsigned int extension);
  /// Set the characteristic size (FWHM) of features in the image (in pixels).
  /// It is used for estimating the correlation length of pixels in the noise and for
  /// boxsize of the local noise estimates.
  /// Size must be given as multiple of 2.
  void setCharacteristicSize(unsigned int size);
  /// Return noise mean \f$\mu_n\f$.
  double getNoiseMean();
  /// Return noise RMS \f$\sigma_n\f$.
  double getNoiseRMS();
  /// Set noise features explicitly.
  void setNoiseMeanRMS(double mean, double rms);
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
  const NumVector<double>& getObjectMap();
   /// Get the segmentation history of this image.
  const History& getHistory ();

 private:
  void insertFormatField(std::string type, std::string columnnr);
  void checkFormat();
  void addFrameBorder(double factor, int& xmin, int& xmax, int& ymin, int& ymax);
  double computeBlendingProbability(unsigned int objectNr);
  void makeLocalNoiseMap(unsigned int nr, int outerxmin, int outerxmax, int outerymin, int outerymax, int innerxmin, int innerxmax, int innerymin, int innerymax, NumVector<double>& mean, NumVector<double>& rms);
  void getSubAreaStats(unsigned int nr, unsigned int startindex, int sidelength, double& mean, double& var);
  int neighborpixel(int pixel,unsigned int direction, int axsize0, int axsize1);
  void findObjectList(std::list<unsigned int>& pixellist, int object, int xmin, int xmax, int ymin, int ymax);
  void removeObjectInnerPixels(std::list<unsigned int>&  pixellist, unsigned int nr);
  void refinePixelMap(std::list<unsigned int>& pixellist, unsigned int startpixel, bool positive, int xmin, int xmax, int ymin, int ymax, double bg_mean, double bg_rms);
  void findHalo(unsigned int id, int& xmin, int& xmax, int& ymin, int& ymax);
  void paintSegmentBorder(int xmin, int xmax, int ymin, int ymax);
  void cleanSegMapArea(int type, int xmin, int xmax, int ymin, int ymax);
  void estimateNoise();
  double distanceFromRim(std::list<unsigned int>& objectlist, unsigned int pixel);
  // Define the format of the SExtractor catalog list.
  // The entries of this struct have all to be filled with the number 
  // of the column of the respective keywords in the SExtractor catalog file.
  struct SExCatFormat {
    unsigned short NUMBER;
    unsigned short XMIN_IMAGE;
    unsigned short XMAX_IMAGE;
    unsigned short YMIN_IMAGE;
    unsigned short YMAX_IMAGE;
    unsigned short BACKGROUND;
    unsigned short FLAGS;
    unsigned short CLASS_STAR;
  };
  // Store the entries of the SExtractor catalog for each object
  struct SExCatObject {
    unsigned short NUMBER;
    int XMIN_IMAGE;
    int XMAX_IMAGE;
    int YMIN_IMAGE;
    int YMAX_IMAGE;
    double BACKGROUND;
    unsigned short FLAGS;
    double CLASS_STAR;
  };
  SExCatFormat sf;
  std::vector<SExCatObject> objectList;
  bool segmapRead, catRead, catChecked, subtractBG, estimatedBG;
  NumVector<double> segMap;
  double bg_mean, bg_rms;
  unsigned int axsize0, axsize1;
  gsl_rng* r;
  unsigned int correlationLength;
  History history;
  std::ostringstream text;
};



#endif
