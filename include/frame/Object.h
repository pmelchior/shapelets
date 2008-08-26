#ifndef OBJECT_H
#define OBJECT_H

#include <History.h>
#include <Typedef.h>
#include <IO.h>
#include <frame/Image.h>
#include <frame/SegmentationMap.h>
#include <frame/PixelCovarianceMatrix.h>
#include <frame/CorrelationFunction.h>
#include <bitset>

/// Central object representing class.
/// The purpose of this class is to faciliate the exchange of object related information
/// between different codes. By our definition, an object is a significant peak of the 
/// brightness distribution, indentified by some image processing code (e.g. Frame).\n\n
/// A respresentant of this class describes an object in real (pixel) space.
/// It consists of the pixel values, the Grid and additional properties (flux, 
/// centroid, ellipticity).
/// In addition, it can store information from the preceding image processing steps:
/// - noise characteristics (including the noise model)
/// - pixel covariance matrix
/// - detection flags
/// - processing history
/// - arbitrary classifier variable
/// - segmentation map
/// 
///\todo Save PixelCovarianceMatrix and CorrelationFunction to Fits file during save()

class Object : public Image<data_t> {
 public:
  /// Argumented constructor.
  /// The \p id is  determined during the segmentation process.
  Object(unsigned long id);
  /// Argumented constructor for loading an object from a Fits file.
  /// The Fits file shold have been created by Object::save().
  Object (std::string fitsfile);
  /// Copy constructor from base type.
  Object (const Image<data_t>& base);
  /// Set the \p id.
  void setID(unsigned long id);
  /// Return the \p id.
  unsigned long getID() const;
  /// Get object classifier.
  /// This can be an arbitrary floating point number, e.g. <tt>CLASS_STAR</tt> from
  /// SExFrame.\n
  /// When the objects come from reading a Catalog (e.g. in SExFrame), CatObject.CLASSIFIER
  /// will be used to set this variable.
  data_t getClassifier() const;
  /// Set the object classifier.
  void setClassifier(data_t classifier);
  /// Return the weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  const NumVector<data_t>& getWeightMap() const;
  /// Access the weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  NumVector<data_t>& accessWeightMap();
  /// Return flux of object.
  /// To get the flux from pixel data, you have to call computeFluxCentroid() before.
  data_t getFlux() const;
  /// Set flux of object.
  void setFlux(data_t F);
  /// Return the position of the object's centroid.
  /// To get the flux from pixel data, you have to call computeFluxCentroid() before.
  const Point2D& getCentroid() const;
  /// Set the position of the object's centroid.
  void setCentroid(const Point2D& xc);
  /// Return 2nd brightness moments.
  /// The 2nd brightness moments are defined relative to the centroid computed
  /// with getCentroid().
  NumMatrix<data_t> get2ndBrightnessMoments();
  /// Return ellipticity of the object.
  /// Ellipticity is calculated from the 2nd brightness moments \f$Q_{ij}\f$ as defined
  /// in Bartelmann & Schneider (2001):\n
  /// \f$ \epsilon = (Q_{11}-Q_{22}+2iQ_{12})/(Q_{11}+Q_{22}+2\sqrt{Q_{11}Q_{22}-Q_{12}^2})\f$
  complex<data_t> getEllipticity();
  /// Return the detection flags.
  /// They indicates problems during the various procedures and are filled by a frameing class.
  /// Thus, look at Frame or SExFrame for the meaning of the individual bits.
  std::bitset<8> getDetectionFlags() const;
  /// Set detection flag (from Frame).
  void setDetectionFlags(const std::bitset<8>& flags);
  /// Return noise mean \f$\mu_n\f$.
  data_t getNoiseMean() const;
  /// Return noise RMS \f$\sigma_n\f$.
  data_t getNoiseRMS() const;
  /// Set noise features from external measurements (e.g. in Frame).
  void setNoiseMeanRMS(data_t mean, data_t rms);
  /// Get the segmentation map.
  const SegmentationMap& getSegmentationMap() const;
  /// Access the segmentation map.
  SegmentationMap& accessSegmentationMap();
  /// Get the pixel covariance matrix.
  const PixelCovarianceMatrix& getPixelCovarianceMatrix() const;
  /// Access the pixel covariance matrix.
  PixelCovarianceMatrix& accessPixelCovarianceMatrix();
  /// Get the correlation function.
  const CorrelationFunction& getCorrelationFunction() const;
  /// Access the correlation function.
  CorrelationFunction& accessCorrelationFunction();
  /// Set the filename from which this object is derived.
  /// This will be stored in the FITS header when using save().
  /// The string should not exceed 80 characters, otherwise it becomes truncated 
  /// in the FITS header.
  void setBaseFilename(std::string filename);
  /// Get the filename from which this object is derived.
  /// see setBaseFilename() for details.
  std::string getBaseFilename() const;
  /// Save the object information in a Fits file.
  /// Data and SegmentationMap will go to pHDU and 1st extHDU, respectively. 
  /// All other information goes to the pHDU header.
  /// If a weight map is provided, these will be stored in the 2nd extHDU.
  void save(std::string fitsfile);
  /// Get History of the object.
  std::string getHistory() const;
  /// Access History of the object.
  History& accessHistory();
  /// Computes the flux and the position of the centroid of the object from pixel data.
  void computeFluxCentroid();

  friend class Frame;
  friend class SExFrame;
  
 private:
  unsigned long id;
  NumVector<data_t> weight;
  SegmentationMap segMap;
  PixelCovarianceMatrix cov;
  CorrelationFunction xi;
  Point2D centroid;
  data_t flux, noise_mean, noise_rms, classifier;
  std::bitset<8> flag;
  std::string basefilename;
  History history;
};
#endif
