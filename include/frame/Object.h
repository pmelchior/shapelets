#ifndef SHAPELENS_OBJECT_H
#define SHAPELENS_OBJECT_H

#include "../Typedef.h"
#include "Image.h"
#include "SegmentationMap.h"
#include "PixelCovarianceMatrix.h"
#include "CorrelationFunction.h"
#include "Moments.h"
#include "../utils/FFT.h"
#include <bitset>

namespace shapelens {

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

class Object : public Image<data_t> {
 public:
  /// Constructor.
  Object();
  /// Argumented constructor for loading an object from a Fits file.
  /// The Fits file shold have been created by Object::save().
  Object (std::string fitsfile);
  /// Copy constructor from base class.
  Object (const Image<data_t>& base);
  /// Copy operator from base class.
  /// It copies only the Image part, but does not replace any other
  /// data members of Object.
  void operator=(const Image<data_t>& base);
  /// Save the object information in a Fits file.
  /// Data and SegmentationMap will go to pHDU and 1st extHDU, respectively. 
  /// All other information goes to the pHDU header.
  /// If a weight map is provided, these will be stored in the extension \p WEIGTH
  /// If a correlation function is provided, these will be stored in the extension \p CORRELATION.
  void save(std::string fitsfile);

  /// Computes Object::flux from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.
  void computeFlux();
  /// Computes Object::centroid from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.\n
  /// \b CAUTION: The results depends on the value of Object::flux.
  void computeCentroid();
  /// Computes Object::flux and Object::centroid from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.
  void computeFluxCentroid();
  /// Computes the quadrupole moment Object::Q from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.\n
  /// \b CAUTION: The results depends on the value of Object::flux and 
  /// Object::centroid.
  void computeQuadrupole();
  /// Computes the octupole moment Object::O from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.\n
  /// \b CAUTION: The results depends on the value of Object::flux and 
  /// Object::centroid.
  void computeOctupole();
  /// Computes the hexadecupole moment Object::H from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.\n
  /// \b CAUTION: The results depends on the value of Object::flux and 
  /// Object::centroid.
  void computeHexadecupole();
  /// Short-hand for computing Object::Q, Object::O, and Object::H.
  /// If Object::weight is non-empty, the pixel weights are considered.\n
  /// \b CAUTION: The results depends on the value of Object::flux and
  /// Object::centroid.
  void computeMoments();
  /// Compute correlation function Object::xi from pixel data.
  /// \p threshold is the minimum significance of the correlation to
  /// to be considered (reasonable value are around 2).
  void computeCorrelationFunction(data_t threshold);
  /// Compute Fourier transform Object::fourier from pixel data.
  void computeFFT();
  /// Convolve with given \p kernel.
  /// For this to work, \p kernel has to have to same size as this object.
  /// If this is not the case, it will be resized accordingly.
  /// Also, kernel::fourier will be computed if the present.
  void convolve(Object& kernel);


  /// The \p id of the Object.
  unsigned long id;
  /// The object classifier.
  /// This can be an arbitrary floating point number, e.g. <tt>CLASS_STAR</tt> from
  /// SExFrame.\n
  /// When the objects come from reading a Catalog (e.g. in SExFrame), CatObject.CLASSIFIER
  /// will be used to set this variable.
  data_t classifier;
  /// The weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  Image<data_t> weight;
  /// The flux of this Object.
  data_t flux;
  /// The position of the object's centroid.
  Point2D<data_t> centroid;
  /// The 2nd brightness moments, defined relative to centroid.
  Quadrupole Q;
  /// The 3rd brightness moments, defined relative to centroid.
  Octupole O;
  /// The 3rd brightness moments, defined relative to centroid.
  Hexadecupole H;
  /// The detection flags.
  /// They indicates problems during the various procedures and are filled by a frameing class.
  /// Thus, look at Frame or SExFrame for the meaning of the individual bits.
  std::bitset<8> flags;
  /// The noise mean \f$\mu_n\f$.
  data_t noise_mean;
  /// The noise RMS \f$\sigma_n\f$.
  data_t noise_rms;
  /// The segmentation map.
  SegmentationMap segMap;
  /// The correlation function.
  CorrelationFunction xi;
  /// The Fourier transform of the pixel data.
  FourierTransform2D fourier;
};
} // end namespace
#endif
