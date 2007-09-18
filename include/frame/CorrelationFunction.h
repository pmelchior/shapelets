#ifndef CORRELATIONFUNCTION_H
#define CORRELATIONFUNCTION_H

#include <map>
#include <NumVector.h>
#include <Typedef.h>
#include <frame/Image.h>
#include <frame/SegmentationMap.h>
#include <frame/Grid.h>

/// Class for computing and storing the pixel correlation function.
/// For computing the correlation function \f$\xi_{ij} = \langle I(x_i) I(x_j)\rangle\f$
/// it is assumed that the pixel correlation is isotropic, thus \f$\xi_{ij} = \xi(r)\f$ with
/// \f$r = |x_i-x_j|\f$.

class CorrelationFunction {
 public:
  /// Default constructor.
  CorrelationFunction ();
  /// Constructor for computing the correlation function from an Image.
  /// The SegmentationMap is used to mask pixels with <tt>segMap(i)!=0</tt>,
  /// unless <tt>mask=0</tt>.\n
  /// Since correlation is expected only on small scales, the averging is done in a
  /// box of area <tt>(2*size+1)^2</tt> centered on each pixel.
  CorrelationFunction (const Image<data_t>& im, const SegmentationMap& segMap, unsigned int size = 2, bool mask=1);
  /// Constructor for computing the correlation function from an NumVector.
  /// The correlation function is computed on all pixels of <tt>data</tt>.\n
  /// Since correlation is expected only on small scales, the averging is done in a
  /// box of area <tt>(2*size+1)^2</tt> centered on each pixel.
  CorrelationFunction (const NumVector<data_t>& data, const Grid& grid, unsigned int size = 2);
  /// Copy constructor.
  void operator= (const CorrelationFunction& xi);
  /// Get correlation function \f$\xi(r)\f$.
  /// This is the mean of the product of pixel values,
  /// \f$\xi(r) = N^{-1}\sum_{i,j}[I(x_i) I(x_j)]\f$ with \f$|x_i-x_j|=r\f$; \f$N(r)\f$
  /// denotes the number of pixels pairs considered.\n
  /// The vector of distances \f$r\f$ is provided by getDistances().
  const NumVector<data_t>& getCorrelationFunction() const;
  /// Access correlation function \f$\xi(r)\f$.
  NumVector<data_t>& accessCorrelationFunction();
  /// Get distances \f$r\f$.
  const NumVector<data_t>& getDistances() const;
  /// Access distances \f$r\f$.
  NumVector<data_t>& accessDistances();
  /// Get the error of the correlation function.
  /// The error provided here is the standard deviation of the mean of the correlation
  /// at the considered distances, 
  /// \f$\Delta\xi(r) = \bigl[\sqrt{N(r)(N(r)-1)}\bigr]^{-1}\bigl[\sum_{i,j}{[I(x_i) I(x_j)-\xi(r)]^2}\bigr]^{\frac{1}{2}}\f$ 
  /// with \f$|x_i-x_j|=r\f$; \f$N(r)\f$ denotes the number of pixels pairs considered.
  const NumVector<data_t>& getCorrelationError() const;
  /// Access the error of the correlation function.
  NumVector<data_t>& accessCorrelationError();
  
 private:
  NumVector<data_t> xi,sigma,dist;
  unsigned int size;
  void makeIndexSetDistances(std::map<unsigned int, unsigned int>& index);
};

#endif
