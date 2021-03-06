#ifndef SHAPELENS_CONFIG_H
#define SHAPELENS_CONFIG_H

#include <string>
#include "Typedef.h"

namespace shapelens {
/// Class for ShapeLens++ configuration parameters.
/// This class stores configuration parameters determining the behaviour of the shapelet
/// decomposition process and the framing process.\n
/// In order to be effective, the parameters have to be set/changed before the method
/// is called, e.g.
/// \code
/// Object obj;
/// ShapeLensConfig::NMAX_HIGH = 24;
/// ShapeletObject(obj);
/// \endcode
///
/// Alternatively, the parameters can be set by a configuration file at construction:
/// \code
/// ShapeLensConfig sconf("shapelens.conf");
/// Frame f(somefile);
/// Object obj;
/// ShapeletObject(obj);
/// \endcode

class ShapeLensConfig {
 public:
  /// Default constructor.
  ShapeLensConfig();
  /// Argumented constructor.
  /// <tt>filename</tt> is expected to contain configuration parameters
  /// in ASCII format, one keyword/value pair per line, separated by a single or multiple
  /// tab character(s):
  /// \code
  /// NMAX_LOW    4
  /// NMAX_HIGH   24
  /// \endcode
  ShapeLensConfig(std::string filename);
  /// Whether the History should be printed to stdout during data processing,
  /// default = 0.
  static bool VERBOSITY;
  /// The minimum number of pixels above <tt>MIN_THRESHOLD</tt> to be regarded as object,
  /// default = 20.
  static unsigned int MIN_PIXELS;
  /// The threshold (in units of the noise variance \f$\sigma_n\f$), which is supposed to 
  /// be the isocontour of minimal significance, default = 1.
  static data_t MIN_THRESHOLD;
  /// The detection threshold (in units of the noise variance \f$\sigma_n\f$), which is the minimum
  /// brightness an oject must have in at least one pixel in order to be detected, default = 2.
  static data_t DETECT_THRESHOLD;
  /// Whether objects with <tt>segMap(i) < 0</tt> are replaced by artificial noise,
  /// default = 0.
  static bool FILTER_SPURIOUS;
  /// Whether the segmented ojects should be checked for overlapping objects
  /// or boundary truncation, default = true.
  static bool CHECK_OBJECT;
  /// The amount by which the objects area is enlarged on each side relative to the
  /// object's area in the segmentation map, default = 0.5
  /// (50% increase in size an all sides).
  static data_t ADD_BORDER;
  /// Lower bound for \f$n_{max}\f$, default = 0.
  static unsigned int NMAX_LOW;
  /// Upper bound for \f$n_{max}\f$, default = 100.
  static unsigned int NMAX_HIGH;
  /// Lower bound for \f$\beta\f$, default = 0.
  static data_t BETA_LOW;
  /// Upper bound for \f$\beta\f$, default = INFINITY.
  static data_t BETA_HIGH;
  /// Relative uncertainty on \f$\beta\f$during optimization, default 0.02.
  static data_t DELTA_BETA;
  /// Whether flattening of \f$\chi^2\f$ is allowed as termination criterium
  /// for the optimization process, default = false.
  static bool ALLOW_FLATTENING;
  /// Whether the value if the shapelet basis functions should be integrated
  /// within pixels when computing the shapelet model, default = false;
  static bool PIXEL_INTEGRATION;
  /// The noise model employed during the decomposition (see Decomposite),
  /// default = "GAUSSIAN"
  static std::string NOISEMODEL;
  /// Whether coordinates are treated as sky coordinates, default = false;
  /// WCS information is read from FITS file header.
  static bool USE_WCS;
};
} // end namespace
#endif
