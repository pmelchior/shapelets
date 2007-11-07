#ifndef SHAPELETOBJECT_H
#define SHAPELETOBJECT_H

#include <fstream>
#include <string>
#include <bitset>
#include <NumMatrix.h>
#include <NumVector.h>
#include <Typedef.h>
#include <History.h>
#include <IO.h>
#include <ShapeLensConfig.h>
#include <frame/Point2D.h>
#include <frame/Grid.h>
#include <frame/Image.h>
#include <frame/Object.h>
#include <frame/Profile.h>
#include <shapelets/Composite2D.h>
#include <shapelets/OptimalDecomposite2D.h>
#include <shapelets/PolarTransformation.h>
#include <shapelets/ImageTransformation.h>
#include <shapelets/MatrixManipulations.h>


/// Central class for 2D shapelet objects.
/// Provides all functionalities related to the work with 2D shapelet objects.\n
/// The shapelet objects can be defined by giving a binary SIF (Shapelet Image Format) file, 
/// cartesian or polar coefficients or a FITS image, that is then decomposed into shapelets.\n\n
/// With the coefficients one can then compose images on arbitrary grids for evaluation
/// and transform the image in shapelet space.\n
/// For efficient work with shapelet images the active set of image parameters can be saved
/// in a SIF file and later be loaded from that file again.

class ShapeletObject : public Composite2D {
 public:
  /// Default constructor.
  ShapeletObject();
  /// Constructor for reading a SIF file.
  /// If <tt>preserve_config==1</tt>, the global ShapeLensConfig parameters
  /// are preserved during the loading of <tt>sifFile</tt>, otherwise
  /// they are overwritten by the values given in <tt>sifFile</tt>.
  ShapeletObject(std::string sifFile, bool preserve_config=1);
  /// Copy constructor.
  //ShapeletObject(ShapeletObject& sobj);
  /// Constructor, using cartesian coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on 
  /// given grid.
  ShapeletObject(const NumMatrix<data_t>& cartesianCoeffs, data_t beta, const Point2D& xcentroid);
  ///  Constructor, using polar coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on
  /// given grid.
  ShapeletObject(const NumMatrix<complex<data_t> >& polarCoeffs, data_t beta, const Point2D& xcentroid);
  /// Constructor for decomposing an Object.
  /// The only thing necessary is a properly filled Object.
  /// The decomposition will find the optimal shapelet parameters automatically.\n
  /// As default, it uses extremely loose bounds on \f$n_{max}\f$ and \f$\beta\f$
  /// (see DEFAULTS); a regularization (to avoid regions of negative flux) is not
  /// employed.\n
  /// If you want to change these settings, change them BEFORE calling this
  /// constructor.
  ShapeletObject(const Object& obj);
  
  /// Set new cartesian coefficients.
  void setCartesianCoeffs(const NumMatrix<data_t>& cartesianCoeffs);
  /// Set cartesian coefficient errors.
  void setCartesianCoeffErrors(const NumMatrix<data_t>& errors);
  /// Set new polar coeficients.
  void setPolarCoeffs(const NumMatrix<complex<data_t> >& polarCoeffs);
  /// Return active cartesian coefficients.
  const NumMatrix<data_t>& getCartesianCoeffs() const;
  /// Return active polar coeficients.
  const NumMatrix<complex<data_t> >& getPolarCoeffs() const;

  // methods depending on the decomposition
  /// Return best fit \f$\chi^2\f$ from decomposition.
  /// It will return 0, if ShapeletObject is not constructed from a Fits file.
  data_t getDecompositionChiSquare() const;
  /// Return error matrix of the cartesian shapelet coefficients.
  /// It will return empty matrix if ShapeletObject is not constructed from a Fits file.
  const NumMatrix<data_t>& getDecompositionErrors() const;

  // methods depending on ImageTransformations
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(data_t rho);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(data_t kappa);
  /// Shear the image by \f$\gamma_0, \gamma_1\f$.
  /// The dimension of the coefficient matrix will be increased by 2.
  void shear(complex<data_t> gamma);
  /// Apply flexion to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void flex(const NumMatrix<data_t>& Dgamma);
  /// Apply lensing operations converge, shear and flex to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void lens(data_t kappa, complex<data_t> gamma, complex<data_t> F, complex<data_t> G);
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
  void convolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel);
  /// Deconvolve the image from another image.
  /// The convolution kernel is given by it's cartesian coefficients and it's beta.
  void deconvolve(const NumMatrix<data_t>& KernelCoeffs, data_t beta_kernel);
  /// Rescale the image.
  /// This changes the coefficients such, that they show the same object with
  /// the new scale size.\n
  /// This is only correct for small changes of \f$\beta\f$.
  /// Larger changes of \f$\beta\f$ can be achieved by setBeta().
  void rescale(data_t newBeta);

  /// Save active set of image parameters to given file.
  void save(std::string filename) const;
  /// Load active set of image parameters from file.
  /// If <tt>preserve_config==1</tt>, the global ShapeLensConfig parameters
  /// are preserved during the loading of <tt>sifFile</tt>, otherwise
  /// they are overwritten by the values given in <tt>sifFile</tt>.
  void load(std::string filename, bool preserve_config = 1);

  /// Return history string of image.
  /// This contains all procedure parameters of decomposition, transformations etc.
  std::string getHistory() const;
  /// Set the history string of the image to an arbitrary string.
  /// This can be used for erasing the history by using a empty string.
  void setHistory(std::string comment);

  /// Get the value of the regularization parameter \f$R\f$.
  /// See OptimalDecomposite2D::regularize() for details.
  data_t getRegularizationR() const;
  /// Get mean of pixel noise \f$\mu_n\f$ obtained from Frame or SExFrame.
  data_t getNoiseMean() const;
  /// Get RMS of pixel noise \f$\sigma_n\f$ obtained from Frame or SExFrame.
  data_t getNoiseRMS() const;
  /// Get the filename from which this object originated.
  std::string getBaseFilename() const;
  /// Get the object id assigned to this object.
  /// See Object::getID() for details.
  unsigned int getObjectID() const;
  /// Get the object number assigned to this object.
  /// See Object::getNumber() for details.
  unsigned int getObjectNumber() const;
  /// Get the object extraction and decomposition flags.
  /// The extraction flags populate the lower 8 bits, the decomposition
  /// flags upper ones.\n
  /// See OptimalDecomposite2D::getDecompositionFlag() and Object::getDetectionFlag()
  /// for details.
  std::bitset<16> getFlags() const;
  
  friend class SIFFile;

 private:
  NumMatrix<data_t> cartesianCoeffs, errors;
  NumMatrix<complex<data_t> > polarCoeffs;
  PolarTransformation c2p;
  ImageTransformation trafo;
  data_t chisquare, R, noise_mean, noise_rms;
  bool fits;
  History history;
  unsigned int nr, id, flags;
  std::string basefilename;
};

#endif
