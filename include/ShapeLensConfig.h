#ifndef SHAPELETCONFIG_H
#define SHAPELETCONFIG_H

#include <string>
#include <Typedef.h>

/// Class for ShapeLens++ configuration parameters.
/// This class stores configuration parameters determining the behaviour of the shapelet
/// decomposition process and the frameing process.\n
/// In order to be effective, the parameters have to be set/changed before the shapelet
/// decomposition is called, e.g.
/// \code
/// Object obj;
/// ShapeLensConfig::NMAX_HIGH = 24;
/// ShapeletObject(obj);
/// \endcode
///
/// Alternatively, the parameters can be set by a configuration file at construction:
/// \code
/// ShapeLensConfig("shapelens.conf");
/// Frame f(somefile);
/// Object obj;
/// ShapeletObject(obj);
/// \endcode;

class ShapeLensConfig {
 public:
  /// Default constructor.
  ShapeLensConfig();
  /// Argumented constructor.
  /// <tt>filename</tt> is expected to contain configureation parameters
  /// in ASCII format, one keyword/value pair per line, separated by a single or multiple
  /// tab character(s):
  /// \code
  /// NMAX_LOW    4
  /// NMAX_HIGH   24
  /// \endcode
  ShapeLensConfig(std::string filename);
  /// Lower bound for \f$n_{max}\f$, default = 0.
  static unsigned int NMAX_LOW;
  /// Upper bound for \f$n_{max}\f$, default = 100.
  static unsigned int NMAX_HIGH;
  /// Lower bound for \f$\beta\f$, default = 0.
  static data_t BETA_LOW;
  /// Upper bound for \f$\beta\f$, default = INFINITY.
  static data_t BETA_HIGH;
  /// Whether a regularization (see OptimalDecomposite2D::regularize()) should 
  /// be employed, default = 0.
  static bool REGULARIZE;
  /// The upper limit for \f$R\f$ during regularizaton, default = 1e-5.
  static data_t REG_LIMIT;
  /// The SIF filename for the unregularized model (empty string if not employed),
  /// default = "".
  static std::string UNREG_SIFFILE;
  /// Whether flattening of \f$\chi^2\f$ is allowed as termination criterium
  /// for the optimization process, default = 0.
  static bool ALLOW_FLATTENING;
  /// Whether objects with <tt>segMap(i) < 0</tt> are replaced by artificial noise,
  /// default = 0.
  static bool FILTER_SPURIOUS;
  /// The amount by which the objects area is enlarged on each side w.r.t. the
  /// objects area in the segmentation map, default = 0.5
  /// (50% increase in size an all sides).
  static data_t ADD_BORDER;
};
#endif
