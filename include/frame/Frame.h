#ifndef FRAME_H
#define FRAME_H

#include <vector>
#include <list>
#include <NumVector.h>
#include <Typedef.h>
#include <History.h>
#include <ShapeLensConfig.h>
#include <frame/Grid.h>
#include <frame/Image.h>
#include <frame/Object.h>
#include <frame/SegmentationMap.h>
#include <frame/Catalog.h>

/// Image preprocessing class.
/// The purpose of this class is to read in a Image and segment it into individual
/// Object entities.\n
/// The image processing includes:
/// - background noise estimation
/// - detection of objects in the frame
/// - segmentation of the frame into segments with one object each
/// - cleaning of the segments from overlapping objects.
///
/// It thus allows the automatic processing of all significant objects in the image.
/// Two typical scenarios are presented here:
/// \code
/// Frame* f = new Frame(fitsfile);
/// f->subtractBackground();
/// NumVector<data_t>& data = f->getData();
/// Grid& grid = f->getGrid();
/// \endcode
/// This reads the extension of a FITS file, estimates the noise statistics
/// subtracts the global background noise level. The data of the whole frame
/// and an adequate grid are returned.\n\n
/// \code
/// Frame* f = new Frame(fitsfile);
/// f->subtractBackground();
/// f->findObjects(50,1.5,5);
/// for(int n = 1; n <= f->getNumberOfObjects(); n++) {
///   Object* obj = new Object(i);
///   f->fillObject(obj);
///   if (obj.getDetectionFlag() < 3) {
///     NumVector<data_t>& data = obj->getData();
///     Grid& grid = obj->getGrid();
///   }
/// }
///\endcode
/// This reads the extension from the FITS file, subtracts the global noise background level
/// and searches for signifcant objects inside the image.\n
/// Then it iterates through the object list, filling the necessary information in the Object
/// \p obj; \p data and \p grid now contain only the segmented frame of the selected object.\n\n
///
/// The Object entities will have these features:
/// - The frames will be quadratic such that the area within the object was located 
/// plus addition border area is included.
/// - Overlapping objects will be masked with noise.
/// - If subtractBackground() has been called, the global background level
/// is subtracted from the object pixels.
/// - \p detectionFlag is set according to the following scheme 
/// (higher values supersed lower ones):
///  - 0: OK
///  - 1: another object nearby, but not overlapping
///  - 2: object close to the image boundary, frame extended with noise, possible cut-off
///  - 4: object cut-off at the image boundary
/// - If a weight map is not given, <tt>noiseModel="GAUSSIAN"</tt>, otherwise 
///   <tt>noisemodel="WEIGHT"</tt>
/// - \p StarGalaxyProbability and \p BlendingProbability are not set.
/// - Its segmentation map is an appropriate cutout of the full segmentation map.
/// - If provided, its weight map is an appropriate cutout of the full weight map.
///
/// The full SegmentationMap of the Image can be obtained from 
/// getSegmentationMap() after calling findObjects().\n\n
/// If detection of blended objects or stars is neccessary, use SExFrame instead.
/// 

class Frame : public Image<data_t> {
 public:
  /// Argumented constructor.
  /// <tt>filename</tt> is the name of the Fits file.\n Extensions or other selections can be 
  // passed in the standard cfitsio way: <tt>filename[extension]</tt>.
  Frame(std::string filename);
  /// Argemented constructor for including a weight map image.
  /// <tt>data_file</tt> is the name of the Fits file containing the data, 
  /// <tt>weight_file</tt> the one of the corresponding weight map.\n
  /// Extensions or other selections can be 
  /// passed in the standard cfitsio way: <tt>data_file[extension]</tt>.
  Frame(std::string data_file, std::string weight_file);
  /// Return noise mean \f$\mu_n\f$.
  /// The mean is obtained by \f$\kappa-\sigma\f$-clipping of the noisy pixel data.
  /// see NumVector::kappa_sigma_clip() for details.
  data_t getNoiseMean();
  /// Return noise RMS \f$\sigma_n\f$.
  /// The RMS is obtained by \f$\kappa-\sigma\f$-clipping of the noisy pixel data.
  /// see NumVector::kappa_sigma_clip() for details.
  data_t getNoiseRMS();
  /// Set noise features explicitly.
  void setNoiseMeanRMS(data_t mean, data_t rms);
  /// Subtract background brightness from image.
  /// The background mean \f$\mu_n\f$ will be taken from the value set by calling
  /// setNoiseMeanRMS() or estimated by running getNoiseMean().
  void subtractBackground();
  /// Identify significant objects in the image.
  /// Searches for connected pixel groups with more than ShapeLensConfig::MIN_PIXELS 
  /// above a lower significance threshold ShapeLensConfig::MIN_THRESHOLD
  /// \f$ \tau_{low} = \mu_n + F_{signific}\cdot\sigma_n\f$
  /// and at least 1 pixel above the higher detection threshold 
  /// ShapeLensConfig::DETECT_TRESHOLD
  /// \f$ \tau_{high} = \mu_n + F_{detect}\cdot\sigma_n\f$.\n
  /// Note that objects cannot be separated if there is a path from the one to the other
  /// that is entirely above \f$\tau_{low}\f$ and that the runtime of the algorithm depends
  /// linearly on the number of detected objects, thus on \f$\tau_{high}\f$.\n
  /// The noise measures \f$\mu_n\f$ and \f$\sigma_n\f$ are the ones defined 
  /// in estimateNoise().\n
  void findObjects();
  /// Return number of objects found during findObjects().
  /// If it returns 0, no objects have been found.
  unsigned int getNumberOfObjects();
  /// Select the object to use for further analysis.
  /// Object 0 is the whole image, Objects 1 to N have to be found in findObjects().
  /// Object 1 is the brightest, all others are essentially unordered.\n
  /// The image will be cut to a region around the selected object. If another objects
  /// is inside this region, its pixel will be replaced by noise according to the
  /// noise statistics derived in estimateBackground().\n
  void fillObject(Object& obj);
  /// Get the SegmentationMap.
  /// The segmentation map defines for each pixel, if it is part of an object.
  /// The convention is:
  /// - 0: noise
  /// - 1..N: objectID of an identified significant pixel group
  SegmentationMap& getSegmentationMap();
  /// Get list of pixels for the given <tt>objectnr</tt>.
  const std::list<unsigned int>& getPixelList(unsigned int objectnr);
  /// Access list of pixels for the given <tt>objectnr</tt>.
  std::list<unsigned int>& accessPixelList(unsigned int objectnr);
  /// Get Catalog of detected objects.
  /// After calling findObject(), only the CatObject.NUMBER is set; if the method
  /// fillObject(Object& obj) is called, all values of the catalog entry <tt>obj.getID()</tt>
  /// are set. Thus, if you want to save the catalog detected by Frame, 
  /// run through all detected objects before.
  const Catalog& getCatalog();
  /// Get the history object of this image.
  const History& getHistory ();

 private:
  data_t noise_mean, noise_rms;
  Image<data_t> weight;
  void estimateNoise();
  void reset();
  void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
  data_t getThreshold(unsigned int pixel, data_t factor);
  SegmentationMap segMap;
  History& history;
  std::vector< std::list<unsigned int> > objectsPixels;
  bool subtractedBG, estimatedBG;
  Catalog catalog;
};

#endif
