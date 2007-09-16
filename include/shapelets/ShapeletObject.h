#ifndef SHAPELETOBJECT_H
#define SHAPELETOBJECT_H

#include <fstream>
#include <string>
#include <NumMatrix.h>
#include <NumVector.h>
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
#include <shapelets/SIFFile.h>


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
  ShapeletObject(std::string sifFile);
  /// Copy constructor.
  ShapeletObject(ShapeletObject& sobj);
  /// Constructor, using cartesian coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on 
  /// given grid.
  ShapeletObject(const NumMatrix<double>& cartesianCoeffs, double beta, const Point2D& xcentroid);
  ///  Constructor, using polar coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on
  /// given grid.
  ShapeletObject(const NumMatrix<complex<double> >& polarCoeffs, double beta, const Point2D& xcentroid);
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
  void setCartesianCoeffs(const NumMatrix<double>& cartesianCoeffs);
  /// Set cartesian coefficient errors.
  void setCartesianCoeffErrors(const NumMatrix<double>& errors);
  /// Set new polar coeficients.
  void setPolarCoeffs(const NumMatrix<complex<double> >& polarCoeffs);
  /// Return active cartesian coefficients.
  const NumMatrix<double>& getCartesianCoeffs();
  /// Return active polar coeficients.
  const NumMatrix<complex<double> >& getPolarCoeffs();

  // methods depending on the decomposition
  /// Return best fit \f$\chi^2\f$ from decomposition.
  /// It will return 0, if ShapeletObject is not constructed from a Fits file.
  double getDecompositionChiSquare();
  /// Return error matrix of the cartesian shapelet coefficients.
  /// It will return empty matrix if ShapeletObject is not constructed from a Fits file.
  const NumMatrix<double>& getDecompositionErrors();

  // methods depending on ImageTransformations
  /// Rotate image counterclockwise by angle \f$\rho\f$.
  void rotate(double rho);
  /// Apply convergence of \f$\kappa\f$ to the image.
  void converge(double kappa);
  /// Shear the image by \f$\gamma_0, \gamma_1\f$.
  /// The dimension of the coefficient matrix will be increased by 2.
  void shear(complex<double> gamma);
  /// Apply flexion to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void flex(const NumMatrix<double>& Dgamma);
  /// Apply lensing operations converge, shear and flex to the image.
  /// The dimension of the coefficient matrix will be increased by 3.
  void lens(double kappa, complex<double> gamma, complex<double> F, complex<double> G);
  /// Translate the image by \f$dx0, dx1\f$.
  /// This is only valid for small translations (less than 1 pixel). 
  /// For larger translations, adjust the centroid by calling setCentroid().\n
  /// \f$dx0, dx1\f$ are assumed to be in pixel scale.
  void translate(double dx0, double dx1);
  /// Circularize the image.
  void circularize();
  /// Flip the image along the X axis
  void flipX();
  /// Brighten the image by the given factor.
  void brighten(double factor);
  /// Convolve the image with another image.
  /// The convolution kernel is given by it's cartesian coefficients and it's beta.
  void convolve(const NumMatrix<double>& KernelCoeffs, double beta_kernel);
  /// Deconvolve the image from another image.
  /// The convolution kernel is given by it's cartesian coefficients and it's beta.
  void deconvolve(const NumMatrix<double>& KernelCoeffs, double beta_kernel);
  /// Rescale the image.
  /// This changes the coefficients such, that they show the same object with
  /// the new scale size.\n
  /// This is only correct for small changes of \f$\beta\f$.
  /// Larger changes of \f$\beta\f$ can be achieved by setBeta().
  void rescale(double newBeta);

  // method for making a profile plot
  /// Return Profile derived from given values from the staring point 
  /// through the ShapeletObject centroid.
  /// The ending point will be symmetric to the starting point.
  void getProfile(const Point2D& start, NumVector<double>& values, int axsize);
  /// Return a Profile from the starting point to the ending point.
  void getProfile(const Point2D& start, const Point2D& end, NumVector<double>& values, int axsize);

  // methods for reading/writing sif files that contain 
  // all necessary information of a shapelet image
  /// Save active set of image parameters to given file.
  void save(std::string filename);
  /// Load active set of image parameters from file.
  void load(std::string filename);

  /// Return history string of image.
  /// This contains all procedure parameters of decomposition, transformations etc.
  std::string getHistory();
  /// Set the history string of the image to an arbitrary string.
  /// This can be used for erasing the history by using a empty string.
  void setHistory(std::string comment);

 private:
  NumMatrix<double> cartesianCoeffs, errors;
  NumMatrix<complex<double> > polarCoeffs;
  PolarTransformation c2p;
  ImageTransformation trafo;
  double chisquare, R;
  bool fits, regularized;
  History history;
  std::ostringstream text;
  char fitsFlag,decompFlag;
};

#endif
