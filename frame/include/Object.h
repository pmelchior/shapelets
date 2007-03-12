#ifndef OBJECT_H
#define OBJECT_H

#include <History.h>
#include <Grid.h>
#include <NumVector.h>
#include <NumMatrix.h>
#include <IO.h>

/// Central object representing class.
/// The purpose of this class is to faciliate the exchange of object related information
/// between different codes. By our definition, an object is a significant peak of the 
/// brightness distribution, indentified by some image processing code (e.g. Frame).\n\n
/// A respresentant of this class describes an object in real (pixel) space.
/// It consists of the pixel values, the Grid and additional properties (flux, 
/// centroid, ellipticity).
/// In addition, it provides information from the preceding image processing steps:
/// - noise characteristics (including the noise model)
/// - detection flags
/// - processing history
/// - blending information
/// - galaxy/star discrimination index

class Object {
 public:
  /// Constructor.
  /// The ID is the object ID determined during the segmentation process
  Object(unsigned int id);
  /// Return the ID.
  unsigned int getID() const;
  /// Return the pixel data of this object.
  const NumVector<double>& getData() const;
  /// Access the pixel data of this object.
  NumVector<double>& accessData();
  /// Return the background mean values in the region of this object.
  const NumVector<double>& getBackground() const;
  /// Access the background mean values in the region of this object.
  /// This information is optional.
  NumVector<double>& accessBackground();
  /// Return the background rms values in the region of this object.
  const NumVector<double>& getBackgroundRMS() const;
  /// Access the background rms values in the region of this object.
  /// This information is optional.
  NumVector<double>& accessBackgroundRMS();
  /// Return the grid, on which the object is defined.
  const Grid& getGrid() const;
  /// Access the grid, on which the object is defined.
  Grid& accessGrid();
  /// Return size in units of the grid stepsize.
  int getSize(bool direction) const;
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
  /// The noise model has to be either
  /// - "GAUSSIAN" or
  /// - "POISSONIAN"
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
  /// Set the filename from which this object is derived.
  /// This will be stored in the FITS header when using save().
  /// The string should not exceed 80 characters, otherwise it becomes truncated 
  /// in the FITS header.
  void setBaseFilename(std::string filename);
  /// Save the object information in a Fits file.
  /// Data and Grid will go into the data unit, all other information into the header.
  /// If background maps are provided, these will be stored in FITS extensions.
  void save(std::string fitsfile);
  /// History of the object
  History history;
  /// Computes the flux and the position of the centroid of the object from pixel data.
  void computeFluxCentroid();

  
 private:
  unsigned int id;
  NumVector<double> data, bg_mean, bg_rms;
  Grid grid;
  Point2D centroid;
  double scaleSize, flux, noise_mean, noise_rms, s_g, blend;
  unsigned short flag, blended;
  std::string noisemodel, basefilename;
};

#endif
