#ifndef SHAPELENS_OBJECT_H
#define SHAPELENS_OBJECT_H

#include "../Typedef.h"
#include "Image.h"
#include "SegmentationMap.h"
#include "PixelCovarianceMatrix.h"
#include "CorrelationFunction.h"
#include "../utils/FFT.h"
#include <bitset>

namespace shapelens {

/// Central object representing class.
/// The purpose of this class is to faciliate the exchange of object related information
/// between different codes. By our definition, an object is a significant peak of the 
/// brightness distribution, indentified by some image processing code (e.g. Frame).\n\n
/// A respresentant of this class describes an object in real (pixel) space.
/// It consists of the pixel values, the Grid, and the centroid.
/// centroid, ellipticity).
/// In addition, it can store information from the preceding image processing steps:
/// - noise characteristics (including the noise model)
/// - pixel covariance matrix
/// - detection flags
/// - processing history
/// - segmentation map
/// - Fourier transform

class Object : public Image<data_t> {
 public:
  /// Constructor.
  Object();
  /// Argumented constructor for loading an object from a Fits file.
  /// The Fits file shold have been created by Object::save().
  Object (std::string fitsfile);
  /// Copy constructor from base class.
  Object (const Image<data_t>& base);
  /// Save the object information in a Fits file.
  /// Data and SegmentationMap will go to pHDU and 1st extHDU, respectively. 
  /// All other information goes to the pHDU header.
  /// If a weight map is provided, these will be stored in the extension \p WEIGTH
  /// If a correlation function is provided, these will be stored in the extension \p CORRELATION.
  void save(std::string fitsfile);

  /// Computes Object::flux and Object::centroid from pixel data.
  /// If Object::weight is non-empty, the pixel weights are considered.
  void computeCentroid();
  /// Compute correlation function Object::xi from pixel data.
  /// \p threshold is the minimum significance of the correlation to
  /// to be considered (reasonable value are around 2).
  void computeCorrelationFunction(data_t threshold);
  /// Compute Fourier transform Object::fourier from pixel data.
  void computeFFT();
  /// Convolve with given \p kernel.
  /// fourier and kernel::fourier will be used if present.
  void convolve(const Object& kernel);


  /// The \p id of the Object.
  unsigned long id;
  /// The weight (inverse variance) map in the region of this object.
  /// This map is employed when <tt>noisemodel==WEIGHT</tt>.
  Image<data_t> weight;
  /// The position of the object's centroid.
  Point<data_t> centroid;
  /// The detection flags.
  /// They indicates problems during the various procedures and are filled by a frameing class.
  /// Thus, look at Frame or SExFrame for the meaning of the individual bits.
  std::bitset<8> flags;
  /// The noise mean \f$\mu_n\f$.
  data_t noise_mean;
  /// The noise RMS \f$\sigma_n\f$.
  data_t noise_rms;
  /// The segmentation map.
  SegmentationMap segmentation;
  /// The correlation function.
  CorrelationFunction xi;
  /// The Fourier transform of the pixel data.
  FourierTransform2D fourier;
};
} // end namespace
#endif
