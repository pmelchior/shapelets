#ifndef SHAPELETOBJECT_H
#define SHAPELETOBJECT_H

#include <fstream>
#include <string>
#include <Point2D.h>
#include <Grid.h>
#include <FitsImage.h>
#include <NumMatrix.h>
#include <NumVector.h>
#include <Composite2D.h>
#include <OptimalDecomposite2D.h>
#include <PolarTransformation.h>
#include <ImageTransformation.h>
#include <IO.h>
#include <MatrixManipulations.h>
#include <Profile.h>
#include <SIFFile.h>
#include <Object.h>
#include <History.h>


/** \mainpage ShapeLens++ Documentation
\section Introduction
The Shapelets formalism is a recent approach to image analysis, especially suited for
observational astronomy.\n
In principle one decomposes an arbitrary 2D function into a series
of localized orthogonal basis functions, which are called 'shapelets'. The specialty of this
basis set is that it consists of weighted Hermite polynomials, which correspond to
perturbations arround a Gaussian. The 'shapelets' are also the basis function of the 2D
quantum harmonic oscillator, and thus allow the use of the operator formalism developed
for this problem. E.g. transformations as rotations, translations etc. can be described
by combinations of lowering and raising operators, acting on the 'shapelets'.\n
Since galaxy objects or their constituents have a localized appearance, they are conveniently
described by using a basis set of localized functions. Thus the number of coefficients needed
to describe a galaxy object well will generally be low. This makes the Shapelets formalism
interesting for data reduction and storage in the first place. Further, as it is possible to
to associate physical quantities of the galaxies to functions of their shapelet coefficients,
one can do the physical analysis in the much smaller shapelet space instead of the real 
space, thus saving memory and computation time.
\section References
- Refregier A., 2003, MNRAS, 338, 35 (later called: Paper I)
- Refregier A., Bacon D., 2003, MNRAS, 338, 48 (Paper II)
- Massey R., Refregier A., 2005, MNRAS, 363, 197 (Paper III)
- Melchior, P. et al., 2006 [astro-ph/0608369] (Paper IV)

\section Additional Libraries
- CCFits: Wrapper for cfitsio library to work with FITS files, in C++ (http://heasarc.gsfc.nasa.gov/docs/software/fitsio/ccfits/).
- GNU Scientific Library (GSL): Free and well tested implementation of common 
(and also not so common) mathematical functions and procedures in C 
(http://www.gnu.org/software/gsl/).
- uBLAS: Numeric vector and matrix implementation in C++ with bindings to Linear Algebra 
 libraries (http://www.boost.org/libs/numeric/ublas/doc/).
- ATLAS: Optimized combination of Linear Algebra routines from BLAS and LAPACK, in C 
 (http://math-atlas.sourceforge.net/).
- LAPACK: Linear Algebra solver routines, in Fortran (http://www.netlib.org/lapack/).\n\n

\section Note
In order to use both ATLAS and LAPACK (since ATLAS unfortunately doesn't support 
all LAPACK routines) one has to rename ATLAS' liblapack.a to liblapack-atlas.a.

\author Peter Melchior (pmelchior at ita dot uni-heidelberg dot de)
*/

/// Central class for 2D shapelet objects.
/// Provides all functionalities related to the work with 2D shapelet objects.\n
/// The shapelet objects can be defined by giving a binary SIF (Shapelet Image Format) file, 
/// cartesian or polar coefficients or a FITS image, that is then decomposed into shapelets.\n\n
/// With the coefficients one can then compose images on arbitrary grids for evaluation
/// and transform the image in shapelet space.\n
/// For efficient work with shapelet images the active set of image parameters can be saved
/// in a SIF file and later be loaded from that file again.

class ShapeletObject : public Composite2D {
  typedef complex<double> Complex;
 public:
  /// Default constructor.
  ShapeletObject();
  /// Constructor for reading a SIF file.
  ShapeletObject(std::string sifFile);
  /// Constructor, using cartesian coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on 
  /// given grid.
  ShapeletObject(const NumMatrix<double>& cartesianCoeffs, double beta, const Point2D& xcentroid);
  ///  Constructor, using polar coefficients.
  /// Define image with given \f$\beta\f$, centroid position \f$x_c\f$ on
  /// given grid.
  ShapeletObject(const NumMatrix<Complex>& polarCoeffs, double beta, const Point2D& xcentroid);
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
  /// Set new polar coeficients.
  void setPolarCoeffs(const NumMatrix<Complex>& polarCoeffs);
  /// Return active cartesian coefficients.
  const NumMatrix<double>& getCartesianCoeffs();
  /// Return active polar coeficients.
  const NumMatrix<Complex>& getPolarCoeffs();

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
  /// The default values for the shapelet decomposition.
  /// The values control the behaviour of the constructor ShapeletObject. As these
  /// values affect the behaviour of the constructor only, they have to be modified
  /// before calling the constructor.\n
  /// Example:
  /// \code
  /// Object obj;
  /// ShapeletObject::DEFAULTS::NMAX_HIGH = 24;
  /// ShapeletObject(obj);
  /// \endcode
  struct DEFAULTS {
    /// Lower bound for \f$n_{max}\f$, default = 0.
    static unsigned int NMAX_LOW;
    /// Upper bound for \f$n_{max}\f$, default = 100.
    static unsigned int NMAX_HIGH;
    /// Lower bound for \f$\beta\f$, default = 0.
    static double BETA_LOW;
    /// Upper bound for \f$\beta\f$, default = INFINITY.
    static double BETA_HIGH;
    /// Whether a regularization (see OptimalDecomposite2D::regularize()) should 
    /// be employed, default = 0.
    static bool REGULARIZE;
    /// The upper limit for \f$R\f$ during regularizaton, default = 1e-5.
    static double REG_LIMIT;
  };

 private:
  NumMatrix<double> cartesianCoeffs, errors;
  NumMatrix<Complex> polarCoeffs;
  PolarTransformation c2p;
  ImageTransformation trafo;
  double chisquare, R;
  bool fits, regularized;
  History history;
  std::ostringstream text;
  char fitsFlag,decompFlag;
};

#endif
