#ifndef SHAPELETCONFIG_H
#define SHAPELETCONFIG_H

#include <string>

/// Class for shapelet configuration parameters.
/// This class stores configuration parameters determining the behaviour of the shapelet
/// decomposition process.\n
/// The parameters can be set by a configuration file at construction.\n\n
/// In order to be effective, the parameters have to be set/changed before the shapelet
/// decomposition is called, e.g.
/// \code
/// Object obj;
/// ShapeletConfig::NMAX_HIGH = 24;
/// ShapeletObject(obj);
/// \endcode
///

class ShapeletConfig {
 public:
  /// Default constructor.
  ShapeletConfig();
  /// Argumented constructor.
  /// <tt>filename</tt> is expected to contain configureation parameters
  /// in ASCII format, one keyword/value pair per line, separated by a single or multiple
  /// tab character(s):
  /// \code
  /// NMAX_LOW    4
  /// NMAX_HIGH   24
  /// \endcode
  ShapeletConfig(std::string filename);
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
  /// The SIF filename for the unregularized model (empty string if not employed),
  /// default = "".
  static std::string UNREG_SIFFILE;
  /// Whether flattening of \f$\chi^2\f$ is allowed as termination criterium
  /// for the optimization process, default = 0.
  static bool ALLOW_FLATTENING;
};
#endif
