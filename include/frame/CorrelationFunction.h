#ifndef CORRELATIONFUNCTION_H
#define CORRELATIONFUNCTION_H

#include <map>
#include <NumVector.h>
#include <Typedef.h>
#include <frame/Image.h>
#include <frame/SegmentationMap.h>
#include <frame/Grid.h>

/// Class for computing and storing the pixel correlation function.
/// The pixel correlation function \f$\xi\f$ is defined as
/// \f[\xi_{\vec\Delta} \equiv N_{\vec\Delta}^{-1}\sum_{i,j}[I(\vec{x}_i) I(\vec{x}_j)]\,\f]
/// where \f$N_{\vec\Delta}\f$ denotes the number of pixels pairs 
/// considered for the given displacement vector \f$\vec{\Delta}\equiv \vec{x}_i - \vec{x}_j\f$.\n
/// This is also applicable for non-isotropic noise correlations.\n\n
/// \b Example
/// \code
/// Frame* f = new Frame(filename);
/// f->subtractBackground();
/// f->findObjects();
///
/// CorrelationFunction corr(*f,f->getSegmentationMap(),5,true);
/// corr.applyThreshold(2);
/// const std::map<Point2D<grid_t>, data_t>& xi = corr.getCorrelationFunction();
/// for (std::map<Point2D<grid_t>, data_t>::const_iterator iter = xi.begin(); iter != xi.end(); iter++)
///   std::cout << iter->first << "\t" << iter->second << std::endl;
/// ...
///  delete f;
/// \endcode
/// This would open a file, detect all objects inside and use the pixel data
/// plus the segmentation to compute the correlation function only pixels
/// which are not considered as object pixels.

class CorrelationFunction {
 public:
  /// Default constructor.
  CorrelationFunction ();
  /// Constructor for computing the correlation function from an Image.
  /// Only pixel correlation with a significance higher than \p threshold 
  /// are considered. The computation is performed on increasingly bigger boxes,
  /// until some correlation become insignificant.\n
  /// If <tt>limit > 0</tt>, the computation is performed in boxes of maximum
  /// sidelength of <tt>2*maxLength+1</tt>.
  CorrelationFunction (const Image<data_t>& data, data_t threshold, int limit = 0);
  /// Constructor for computing the correlation function from an Image.
  /// The SegmentationMap is used to mask pixels with <tt>segMap(i)!=0</tt>.\n
  /// Only pixel correlation with a significance higher than \p threshold 
  /// are considered. The computation is performed on increasingly bigger boxes,
  /// until some correlation become insignificant.\n
  /// If <tt>limit > 0</tt>, the computation is performed in boxes of maximum
  /// sidelength of <tt>2*maxLength+1</tt>.
  CorrelationFunction (const Image<data_t>& im, const SegmentationMap& segMap, data_t threshold, int limit = 0);
  /// Constructor from a correlation matrix;
  /// \p corr must have the layout used by getCorrelationMatrix().
  CorrelationFunction (const NumMatrix<data_t>& corr);
  /// Get correlation function \f$\xi\f$ as matrix.
  /// Returns a <tt>(2*maxLength+1)^2</tt> matrix, where the indices encode
  /// the spatial information of \f$\vec\Delta\f$.
  /// The pixel with \f$\xi_{\vec 0}\f$ is in the center.
  NumMatrix<data_t> getCorrelationMatrix() const;
  /// Get correlation function \f$\xi\f$.
  const std::map<Point2D<grid_t>, data_t>& getCorrelationFunction() const;
  /// Get the error \f$\sigma(\xi)\f$ of the correlation function.
  /// The error provided here is the standard deviation of the mean of the correlation
  /// at the considered displacements, 
  /// \f$\sigma(\xi)_{\vec\Delta} = \bigl[\sqrt{N(N-1)}\bigr]^{-1}\bigl[\sum_{i,j}{[I(\vec{x}_i) I(\vec{x}_j)-\xi_{\vec\Delta}]^2}\bigr]^{\frac{1}{2}}\f$ 
  /// with  \f$\vec\Delta\f$ and \f$N_{\vec\Delta}\f$ as defined above.
  const std::map<Point2D<grid_t>, data_t>& getCorrelationError() const;
  /// Apply significance threshold to entries of \f$\xi\f$.
  /// Keep only entries, where \f$\xi_i \geq \tau \sigma(\xi)_i\f$.
  void applyThreshold(data_t tau);
  /// Get maximal length used for constructing \f$\xi\f$.
  unsigned int getMaxLength() const;
  /// Get amount of pixels relevant for \f$\xi\f$.
  unsigned int getSize() const;
 private:
  void compute(const Image<data_t>& data);
  void compute(const Image<data_t>& im, const SegmentationMap& segMap);
  std::map<Point2D<grid_t>, data_t> xi, sigma;
  std::map<Point2D<grid_t>, unsigned int> num;
  int maxLength;
  void setPoints();
};

#endif
