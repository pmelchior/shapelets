#ifndef OBJECT_H
#define OBJECT_H

#include <frame/History.h>
#include <frame/Image.h>
#include <frame/SegmentationMap.h>
#include <frame/PixelCovarianceMatrix.h>
#include <frame/CorrelationFunction.h>
#include <frame/IO.h>

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
/// - blending information
/// - galaxy/star discrimination index
/// - segmentation map
/// 
///\todo Save PixelCovarianceMatrix and CorrelationFunction to Fits file during save()

class Object : public Image<double> {
 public:
  /// Argumented constructor.
  /// The ID is the object ID determined during the segmentation process
  Object(unsigned int id);
  /// Argumented constructor for loading an object from a Fits file.
  /// The Fits file shold have been created by Object::save().
  Object (std::string fitsfile);
  /// Return the ID.
  unsigned int getID() const;
  /// Set the object number.
  /// This can but does not have to be identical to the <tt>ID</tt> of this object,
  /// this number can thus carry an additional identifier.
  void setNumber(unsigned int num);
  /// Get object number.
  unsigned int getNumber() const;
  /// Return the weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  const NumVector<double>& getWeightMap() const;
  /// Access the weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  NumVector<double>& accessWeightMap();
  /// Return flux of object.
  /// To get the flux from pixel data, you have to call computeFluxCentroid() before.
  double getFlux() const;
  /// Set flux of object.
  void setFlux(double F);
  /// Return the position of the object's centroid.
  /// To get the flux from pixel data, you have to call computeFluxCentroid() before.
  const Point2D& getCentroid() const;
  /// Set the position of the object's centroid.
  void setCentroid(const Point2D& xc);
  /// Return 2nd brightness moments.
  /// The 2nd brightness moments are defined relative to the centroid computed
  /// with getCentroid().
  NumMatrix<double> get2ndBrightnessMoments();
  /// Return ellipticity of the object.
  /// Ellipticity is calculated from the 2nd brightness moments \f$Q_{ij}\f$ as defined
  /// in Bartelmann & Schneider (2001):\n
  /// \f$ \epsilon = (Q_{11}-Q_{22}+2iQ_{12})/(Q_{11}+Q_{22}+2\sqrt{Q_{11}Q_{22}-Q_{12}^2})\f$
  complex<double> getEllipticity();
  /// Return the detection flag.
  /// It indicates problems during the various procedures. 
  /// Higher numbers supersede lower ones.
  /// - 0: OK
  /// - 1: Object nearby, but not overlapping
  /// - 2: Object close to the image boundary, frame extended with noise, possible cut-off
  /// - 3: Object overlapped by other object: BLENDED = 1
  /// - 4: Object cut-off at the image boundary
  /// 
  /// Images with flag > 2 should not be used for further analysis.
  unsigned short getDetectionFlag() const;
  /// Set detection flag (from Frame).
  void setDetectionFlag(unsigned short flag);
  /// Return noise mean \f$\mu_n\f$.
  double getNoiseMean() const;
  /// Return noise RMS \f$\sigma_n\f$.
  double getNoiseRMS() const;
  /// Set noise features from external measurements (e.g. in Frame).
  void setNoiseMeanRMS(double mean, double rms);
  /// Set the noise model.
  void setNoiseModel(std::string noisemodel);
  /// Get the noise model.
  std::string getNoiseModel() const;
  /// Get the blending probability.
  /// If the probability is small, the object is considered to be unblended.
  double getBlendingProbability() const;
  /// Set the blending probability.
  void setBlendingProbability(double blend);
  /// Get star/galaxy discrimination index.
  /// The number is normalized, high values are given for objects, which are very 
  /// likely stars, low values for galaxies.
  double getStarGalaxyProbability() const;
  /// Set the star/galaxy discrimination index.
  void setStarGalaxyProbability(double s_g);
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
  /// Save the object information in a Fits file.
  /// Data and SegmentationMap will go to pHDU and 1st extHDU, respectively. 
  /// All other information goes to the pHDU header.
  /// If a weight map is provided, these will be stored in the 2nd extHDU.
  void save(std::string fitsfile);
  /// History of the object
  History history;
  /// Computes the flux and the position of the centroid of the object from pixel data.
  void computeFluxCentroid();

  
 private:
  unsigned int id, number;
  NumVector<double> weight;
  SegmentationMap segMap;
  PixelCovarianceMatrix cov;
  CorrelationFunction xi;
  Point2D centroid;
  double flux, noise_mean, noise_rms, s_g, blend;
  unsigned short flag;
  std::string noisemodel, basefilename;
};

#endif
