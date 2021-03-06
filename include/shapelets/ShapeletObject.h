#ifndef SHAPELENS_SHAPELETOBJECT_H
#define SHAPELENS_SHAPELETOBJECT_H

#include <fstream>
#include <string>
#include <bitset>
#include <numla/NumMatrix.h>
#include <numla/NumVector.h>
#include "../Typedef.h"
#include "../ShapeLensConfig.h"
#include "../utils/Property.h"
#include "../utils/SourceModel.h"
#include "../frame/Object.h"
#include "Composite.h"
#include "OptimalDecomposite.h"
#include "PolarTransformation.h"
#include "ImageTransformation.h"

namespace shapelens {

/// Central class for 2D shapelet objects.
/// Provides all functionalities related to the work with 2D shapelet objects.\n
/// The shapelet objects can be defined by giving a SIFFile (Shapelet Image Format) name, 
/// cartesian or polar coefficients or a Object entity, which is then decomposed into shapelets.\n\n
/// With the coefficients one can then compose shapelet models (Composite entities) 
/// and employ transformations in shapelet space.\n
/// For efficient work with shapelet images the active set of image parameters can be saved
/// to a SIFFile and later be loaded from that file again.

class ShapeletObject : public Composite {
 public:
  /// Default constructor.
  ShapeletObject();
  /// Copy constructor.
  ShapeletObject(const ShapeletObject&);
  /// Constructor for reading a SIF file.
  /// If <tt>preserve_config==1</tt>, the global ShapeLensConfig parameters
  /// are preserved during the loading of <tt>sifFile</tt>, otherwise
  /// they are overwritten by the values given in <tt>sifFile</tt>.
  ShapeletObject(std::string sifFile, bool preserve_config=1);
  /// Constructor, using cartesian coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on 
  /// given grid.
  ShapeletObject(const CoefficientVector<data_t>& cartesianCoeffs, data_t beta, const Point<data_t>& xcentroid = Point<data_t>(0,0), const Grid& grid = Grid(-25,-25,50,50));
  ///  Constructor, using polar coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on
  /// given grid.
  ShapeletObject(const CoefficientVector<std::complex<data_t> >& polarCoeffs, data_t inbeta, const Point<data_t>& xcentroid = Point<data_t>(0,0), const Grid& grid = Grid(-25,-25,50,50));
  /// Constructor for decomposing an Object.
  /// The only thing necessary is a properly filled Object.
  /// The decomposition will find the optimal shapelet parameters automatically,
  /// within the bounds given by the ShapeLensConfig.\n\n
  /// \b ATTENTION: If you want to change the configuration settings, you have to 
  /// change them \e BEFORE calling this constructor.
  ShapeletObject(const Object& obj);
  ///Copy operator.
  ShapeletObject& operator=(const ShapeletObject&);

  /// Set cartesian coefficients.
  void setCoeffs(const CoefficientVector<data_t>& cartesianCoeffs);
  /// Set polar coefficients.
  void setPolarCoeffs(const CoefficientVector<std::complex<data_t> >& polarCoeffs);
  /// Set errors of the cartesian coefficients.
  /// The given coefficient errors will populate the diagonal elements of 
  /// the coefficient covariance matrix.
  void setErrors(const CoefficientVector<data_t>& errors);

  /// Get polar coefficients.
  CoefficientVector<std::complex<data_t> > getPolarCoeffs();
  /// Get polar coefficients.
  CoefficientVector<std::complex<data_t> > getPolarCoeffs() const;

  // methods depending on the decomposition
  /// Return best fit \f$\chi^2\f$ from decomposition.
  /// It will return 0, if ShapeletObject is not constructed from a Fits file.
  data_t getChiSquare() const;
  /// Get the cartesian coefficient errors.
  /// The coefficient errors are computed from the covariance matrix 
  /// (cf. getCovarianceMatrix()) by quadratic summation:
  /// \f[ \Delta I_{n_1,n_2} = \sqrt{\sum_{n_1,n_2} cov(n_1,n_2)}\f]
  CoefficientVector<data_t> getErrors() const;

  // methods depending on ImageTransformations
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(data_t rho);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(data_t kappa);
  /// Shear the image by \f$\gamma_0, \gamma_1\f$.
  /// The dimension of the coefficient matrix will be increased by 2.
  void shear(std::complex<data_t> gamma);
  /// Apply flexion to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void flex(std::complex<data_t> F, std::complex<data_t> G);
  /// Apply lensing operations converge, shear and flex to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void lens(data_t kappa, std::complex<data_t> gamma, std::complex<data_t> F, std::complex<data_t> G);
  /// Translate the image by \f$dx0, dx1\f$.
  /// This is only valid for small translations (less than 1 pixel). 
  /// For larger translations, adjust the centroid by calling setCentroid().\n
  /// \f$dx0, dx1\f$ are assumed to be in pixel scale.
  void translate(data_t dx0, data_t dx1);
  /// Circularize the image.
  void circularize();
  /// Flip the image along the X axis
  void flipX();
  /// Brighten the image by the given factor.
  void brighten(data_t factor);
  /// Convolve the image with another image.
  /// The convolution kernel is given by it's cartesian coefficients and it's beta.
  void convolve(const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel);
  /// Deconvolve the image from another image.
  /// The convolution kernel is given by it's cartesian coefficients and it's beta.
  void deconvolve(const CoefficientVector<data_t>& kernelCoeffs, data_t beta_kernel);
  /// Rescale the image.
  /// This changes the coefficients such, that they show the same object with
  /// the new scale size.\n
  /// This is only correct for small changes of \f$\beta\f$.
  /// Larger changes of \f$\beta\f$ can be achieved by setBeta().
  void rescale(data_t newBeta);

  /// Save active set of image parameters to given SIFFile.
  void save(std::string sifFile) const;
  /// Load active set of image parameters from SIFFile.
  /// If <tt>preserve_config==1</tt>, the global ShapeLensConfig parameters
  /// are preserved during the loading of <tt>sifFile</tt>, otherwise
  /// they are overwritten by the values given in <tt>sifFile</tt>.
  void load(std::string sifFile, bool preserve_config = 1);

  /// Return history string of image.
  /// This contains all procedure parameters of decomposition, transformations etc.
  std::string getHistory() const;
  /// Set the history string of the image to an arbitrary string.
  /// This can be used for erasing the history by using a empty string.
  void setHistory(std::string comment);

  /// Get the filename from which this object originated.
  std::string getBaseFilename() const;
  /// Get the object id assigned to this object.
  /// See Object::getID() for details.
  unsigned long getObjectID() const;
  /// Get the object extraction and decomposition flags.
  /// The extraction flags populate the lower 8 bits, the decomposition
  /// flags upper ones.\n
  /// See OptimalDecomposite::getDecompositionFlag() and Object::getDetectionFlag()
  /// for details.
  const std::bitset<16>& getFlags() const;
  /// Flexible storage container.
  Property prop;
  
  // legacy functions:
  void setName(std::string name);
  std::string getName() const;
  void setTag(data_t tag);
  void setID(unsigned long id);
  data_t getTag() const;
  data_t getObjectClassifier() const;
  
  friend class SIFFile;
  friend class ShapeletObjectList;

  private:
  CoefficientVector<data_t>& coeffs;
  NumMatrix<data_t>& cov;
  CoefficientVector<std::complex<data_t> > polarCoeffs;
  NumMatrix<std::complex<data_t> > polarCov;
  PolarTransformation c2p;
  ImageTransformation trafo;
  data_t chisquare;
  bool fits, updatePolar;
  History history;
  std::bitset<16> flags;
  unsigned long id;
  std::string basefilename, name;
  void correctCovarianceMatrix();

  // Transform arbitrary matrix into triangular matrix of appropraite dimension.
  // Copies entries from input matrix into lower left corner of the triangular matrix
  // as long as there are entries unequal to 0 on the diagonal. 
  template <class T>
    CoefficientVector<T> triangularizeCoeffs(const NumMatrix<T>& input, data_t cutoff=0) {
    int nmax = getMaxOrder(input,cutoff);
    CoefficientVector<T> cv(nmax);
    unsigned int n1, n2;
    const IndexVector& nVector = cv.getIndexVector();
    for (unsigned int i=0; i < cv.size(); i++) {
      n1 = nVector.getState1(i);
      n2 = nVector.getState2(i);
      cv(i) = input(n1,n2);
    }
    return cv;
  }

  // Get maximal order to keep for the triangular matrix.
  // Reads the matrix from top right to bottom left along the diagonals and returns
  // the number of the diagonal where not all entries are 0.
  template <class T>
    int getMaxOrder(const NumMatrix<T>& matrix, data_t cutoff) {
    int result;
    // start in lower right corner and work back to top left corner
    for (int diag = (matrix.getRows() + matrix.getColumns() -2); diag>=0; diag-- ) {
      result = diag;
      bool isempty = 1;
      int startj, stopj;
      if (diag >= matrix.getRows()-1) {
	startj = diag - matrix.getRows() + 1;
	stopj = matrix.getColumns() - 1;
      } else {
	startj = 0;
	stopj = diag;
      }
      for (int j = startj; j <= stopj; j++) {
	int i = diag -j;
	if (abs(matrix(i,j)) > cutoff) {
	  isempty = 0;
	  break;
	}
      }
      if (!isempty) break;
    }
    return result;
  }
};

/// Model from ShapeletObject.
/// The class provides convenient sampling from a ShapeletObject.
class ShapeletModel : public SourceModel {
public:
  /// Constructor.
  /// By giving \p centroid, one moves \p sobj to this centroid and adjusts
  /// the support without altering \p sobj itself. Similarly, \p flux 
  /// automatically rescales \p sobj to the desired total flux.\n
   ShapeletModel(const boost::shared_ptr<ShapeletObject>& sobj, data_t flux, const CoordinateTransformation* CT = NULL);
  /// Sample model at \p P.
  virtual data_t getValue(const Point<data_t>& P) const;
  /// Get total flux of model.
  virtual data_t getFlux() const;
  /// Get type of model.
  virtual char getModelType() const;
private:
  boost::shared_ptr<ShapeletObject> sobj;
  const Point<data_t>& scentroid;
  data_t flux_scale;
};

} // end namespace
#endif
