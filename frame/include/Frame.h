#ifndef FRAME_H
#define FRAME_H

#include <list>
#include <Grid.h>
#include <NumVector.h>
#include <FitsImage.h>
#include <Object.h>
#include <History.h>

/// Image preprocessing class.
/// The purpose of this class is to read in a FitsImage and segment it into individual
/// Object entities.\n
/// The image processing includes (features in brackets are to come):
/// - background noise estimation
/// - detection of objects in the frame
/// - segmentation of the frame into segments with one object each
/// - cleaning of the segments from disturbing noise or overlapping objects.
/// - (identification and deblending of blended objects)
/// - (galaxy/star discrimination)
///
/// It thus allows the automatic processing of all significant objects in the image.
/// Two typical scenarios are presented here:
/// \code
/// Frame* f = new Frame(fitsfile);
/// f->subtractBackground();
/// NumVector<double>& data = f->getData();
/// Grid& grid = f->getGrid();
/// \endcode
/// This simply reads the extension of a FITS file, estimates the noise statistics
/// globally and subtracts the global background noise level. The data of the whole frame
/// and the adequate grid are returned.\n\n
/// \code
/// Frame* f = new Frame(fitsfile);
/// f->subtractBackground();
/// f->findObjects(50,1.5,5);
/// for(int n = 1; n <= f->getNumberOfObjects(); n++) {
///   Object* obj = new Object(i);
///   f->fillObject(obj);
///   if (obj.getDetectionFlag() < 3) {
///     NumVector<double>& data = obj->getData();
///     Grid& grid = obj->getGrid();
///   }
/// }
///\endcode
/// This reads the extension from the FITS file, subtracts the global noise background level
/// and searches for signifcant objects inside the image.\n
/// Then it iterates through the object list, filling the necessary information in a Object.
/// data and grid now contain only the segmented frame of the selected object.\n\n
///
/// The Object entities will have these features:
/// - The frames will be quadratic such that the area within the object was located 
/// is completely included and addition border areas is included.
/// - Overlapping objects will be masked with noise.
/// - If subtractBackground() has been called the global background level from estimateNoise()
/// is subtracted from the object pixels.
/// - \p detectionFlag is set according to the output of getProcessingFlag().
/// - \p noiseModel is "POISSONIAN".
/// - \p StarGalaxyProbability and \p BlendingProbability are not yet set.
/// 
/// In addition, the segmentation map of the FitsImage can be obtained from getObjectMap()
/// after calling findObjects().\n
/// If detection of blended objects or stars is neccessary, use SExFrame instead.
/// 
/// \todo - blending, galaxy/star discrimination


class Frame : public FitsImage<double> {
 public:
  /// Argumented constructor.
  /// The filename of the Fits file and extension have to given
  /// in the standard cfitsio way: fitsfile[extension].
  Frame(std::string filename);
  /// Estimate the statistical noise features.
  /// The method used is \f$\sigma\f$-clipping (as in SExtractor). This iteratively selects
  /// only those pixels with values in a 3\f$\sigma\f$ interval around the median, until
  /// convergence is reached (no further pixels discarded).\n
  /// If the selected sample is not strongly skewed (skewness < 0.1), the median and the
  /// standard deviation of this sample will be used as noise mean \f$\mu_n\f$ and sigma
  /// \f$\sigma_n\f$, respectively.\n
  /// Otherwise a histogram of the selection will be constructed (with 100 pixels per bin,
  /// on average). Assuming that the left wing of the histogram is unbiased by the brighter
  /// objects pixels, the right wing will be set to the appropriate value on the left side,
  /// whenever it becomes significantly too bright. Then the mean and the standard deviation
  /// of the histogram will be used as \f$\mu_n\f$ and \f$\sigma_n\f$.\n\n
  /// If the method is called for an objectID > 0, the noise features will be estimated
  /// locally inside each objects frame. Be aware, that local noise estimation is much more
  /// subject to small scale variations and thus easily overestimates the noise measures.
  void estimateNoise();
  /// Subtract background brightness from image.
  /// The subtracted background level is set to the noise mean \f$\mu_n\f$, estimated in
  /// estimateNoise(), which is implicitly called before.\n
  /// If the background should be subtracted locally around a particular object, this
  /// objects has to be selected and estimateNoise() has to be called for this
  /// object.\n
  /// If the background has been subtracted globally,
  /// it cannot be subtracted locally anymore.
  void subtractBackground();
  /// Return noise mean \f$\mu_n\f$.
  double getNoiseMean();
  /// Return noise RMS \f$\sigma_n\f$.
  double getNoiseRMS();
  /// Set noise features explicitly.
  void setNoiseMeanRMS(double mean, double rms);

  /// Identify significant objects in the image.
  /// Searches for connected pixel groups with more than minPixels 
  /// above a lower significance threshold
  /// \f$ \tau_{low} = \mu_n + F_{signific}\cdot\sigma_n\f$
  /// and at least 1 pixel above a higher detection threshold
  /// \f$ \tau_{high} = \mu_n + F_{detect}\cdot\sigma_n\f$.\n
  /// Note that objects can not be separated if there is a path from the one to the other
  /// that is entirely above \f$\tau_{low}\f$ and that the runtime of the algorithm depends
  /// linearly on the number of detected objects, thus on \f$\tau_{high}\f$.\n
  /// The noise measures \f$\mu_n\f$ and \f$\sigma_n\f$ are the ones defined 
  /// in estimateNoise().\n
  /// If \f$F_{signific}<1\f$ or \f$F_{detect}<1\f$, the default values 1.5 and 3 are used.
  /// If minPixels == 0, 50 pixels are used instead.
  void findObjects(unsigned int minPixels, double F_signific, double F_detect);
  /// Return number of objects found during findObjects().
  /// If it returns 0, no objects have been found.
  unsigned int getNumberOfObjects();
  /// Select the object to use for further analysis.
  /// Object 0 is the whole image, Objects 1 to N have to be found in findObjects().
  /// Object 1 is the brightest, all others are unordered.\n
  /// The image will be cut to a region around the selected object. If another objects
  /// is inside this region, its pixel will be replaced by noise according to the
  /// noise statistics derived in estimateBackground().
  void fillObject(Object& obj);
  /// Get the map of found objects.
  /// The object map defines for each pixel, if it is part of an object.
  /// The convention is:
  /// - -1: non-significant pixel group, meaning above noise, but too small or faint
  /// - -2: positive small-scale noise oscillation bigger than 8 pixels close to detected object
  /// - -3: negative small-scale noise oscillation bigger than 8 pixels close to detected object
  /// - -4: borderline of the segments
  /// - 0: noise
  /// - 1..N: objectID of an identified significant pixel group
  NumVector<int>& getObjectMap();
  // Get list of pixels for the actual object.
  //std::list<int>& getPixelList();
  /// Get flag from the various image processing steps
  /// Higher values supersede lower ones.
  /// - 0: OK
  /// - 1: Image without noise 
  /// - 2: Noise estimation corrected by using symmetrized histogram
  /// - 3: Noise estimate did not converge 
  unsigned short getProcessingFlag();
  /// Get the history object of this image.
  const History& getHistory ();

 private:
  double noise_mean, noise_rms;
  void reset();
  int neighborpixel(int pixel,unsigned int direction);
  void definePixelMap(int startpixel, int counter, double threshold);
  void refinePixelMap(int startpixel, bool positive, int xmin, int xmax, int ymin, int ymax);
  void addFrameBorder(int& xmin, int& xmax, int& ymin, int& ymax);
  void paintSegmentBorder(int xmin, int xmax, int ymin, int ymax);
  NumVector<int> pixelmap;
  std::list<int> pixellist;
  bool subtractedBG, estimatedBG, foundObjects;
  unsigned short flag;
  unsigned int numberofObjects;
  History history;
  std::ostringstream text;
};

#endif
