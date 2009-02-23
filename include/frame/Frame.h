#ifndef FRAME_H
#define FRAME_H

#include <vector>
#include <set>
#include <NumVector.h>
#include <Typedef.h>
#include <ShapeLensConfig.h>
#include <utils/tree.hh>
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
/// - detection of blended substructure within objects
/// - segmentation of the frame into segments with one object each
/// - cleaning of the segments from overlapping objects.
///
/// It thus allows the automatic processing of all significant objects in the image.
/// A typical scenario is presented here:
/// \code
/// Frame* f = new Frame(fitsfile);
/// f->subtractBackground();
/// f->findObjects();
/// Catalog& cat = f->getCatalog();
/// Catalog::iterator iter;
/// for(iter = cat.begin(); iter != cat.end(); iter++) {
///   unsigned long id = (*iter).first;
///   Object obj;
///   f->fillObject(obj,iter);
///   if ((*iter).second.FLAGS == 0) {
///     // do something useful with obj ...
///   }
/// }
///\endcode
/// This reads the extension from the FITS file, subtracts the global noise background level
/// and searches for signifcant objects inside the image. The minum number of pixels or the 
/// detection threshold etc. are obtained from ShapeLensConfig settings.\n
/// Then it iterates through the Catalog entries. If no problems with this object have been found 
/// (<tt>FLAGS == 0</tt>), it fills the necessary information in the Object container \p obj.
///
/// The Object entities will have these features:
/// - The frames will be quadratic such that the area within the object was located 
/// plus addition border area is included.
/// - Overlapping objects will be masked with noise.
/// - If subtractBackground() has been called, the global background level
/// is subtracted from the object pixels.
/// - \p detectionFlags(i) are set according to the following scheme, which is similar to 
/// the one of SExtractor.
///  - <tt>i = 0</tt>: another object nearby, but not overlapping
///  - <tt>i = 1</tt>: object blended
///  - <tt>i = 2</tt>: object close to the image boundary, frame extended with noise, possible cut-off
///  - <tt>i = 3</tt>: object cut-off at the image boundary
/// - If a weight map is not given, <tt>noiseModel="GAUSSIAN"</tt>, otherwise 
///   <tt>noisemodel="WEIGHT"</tt>
/// - \p StarGalaxyProbability is not set.
/// - Its segmentation map is an appropriate cutout of the full segmentation map.
/// - If provided, its weight map is an appropriate cutout of the full weight map.
///
/// The full SegmentationMap of the Image can be obtained from getSegmentationMap() after 
/// calling findObjects(); a Catalog of the objects therein is also produced at that time, but
/// quantities like centroid positions are only valid after calling fillObject() for the appropriate
/// Object.

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
  unsigned long getNumberOfObjects();
  /// Select the object to use for further analysis.
  /// The object is selected by the iterator \p catiter which must must
  /// be constructed from the Catalog obtained by calling getCatalog().\n
  /// The image will be cut to a region around the selected object. 
  /// If another objects is inside this region, its pixel will be replaced 
  /// by noise according to the noise statistics derived in estimateBackground().
  void fillObject(Object& obj, Catalog::iterator& catiter);
  /// Get the SegmentationMap.
  /// The segmentation map defines for each pixel, if it is part of an object.
  /// The convention is:
  /// - 0: noise
  /// - 1..N: objectID of an identified significant pixel group
  const SegmentationMap& getSegmentationMap();
  /// Get list of pixels for the given <tt>objectnr</tt>.
  const std::set<unsigned long>& getPixelSet(unsigned long objectnr);
  /// Get Catalog of detected objects.
  /// After calling findObject(), only the CatObject.NUMBER is set; if the method
  /// fillObject(Object& obj) is called, all values of the catalog entry <tt>obj.getID()</tt>
  /// are set. Thus, if you want to save the catalog detected by Frame, 
  /// run through all detected objects before.
  Catalog& getCatalog();
  /// Get the history object of this image.
  const History& getHistory ();
  /// Compute correlation function from the pixel data of the entire Frame.
  /// \p threshold is the minimum significance of the correlation to
  /// to be considered (reasonable value are around 2).
  CorrelationFunction computeCorrelationFunction(data_t threshold);

 private:
  data_t noise_mean, noise_rms;
  Image<data_t> weight;
  void estimateNoise();
  void reset();
  void addFrameBorder(data_t factor, int& xmin, int& xmax, int& ymin, int& ymax);
  data_t getThreshold(unsigned long pixel, data_t factor);
  void linkPixels(std::set<unsigned long>& pixelset, data_t& max, data_t& max_threshold, unsigned long startpixel);
  bool detectBlending(const std::set<unsigned long>& all, data_t max, data_t max_threshold);
  void insertNodesAboveThreshold(tree<std::set<unsigned long> >& tree, tree<std::set<unsigned long> >::iterator_base& diter, data_t threshold);
  data_t getTotalFlux(const std::set<unsigned long>& pixels);
  SegmentationMap segMap;
  History& history;
  std::map< unsigned long, std::set<unsigned long> > objectsPixels;
  bool subtractedBG, estimatedBG;
  Catalog catalog;
};

#endif
