#ifndef FRAMECONFIG_H
#define FRAMECONFIG_H

#include <string>

/// Class for frame configuration parameters.
/// This class stores configuration parameters determining the behaviour of the framing
/// process.\n
/// The parameters can be set by a configuration file at construction.\n\n
/// In order to be effective, the parameters have to be set/changed before the frame
/// entity is called, e.g.
/// \code
/// FrameConfig("frame.conf");
/// Frame f(somefitsfile);
/// \endcode
///

class FrameConfig {
 public:
  /// Default constructor.
  FrameConfig();
  /// Argumented constructor.
  /// <tt>filename</tt> is expected to contain configureation parameters
  /// in ASCII format, one keyword/value pair per line, separated by a single or multiple
  /// tab character(s):
  /// \code
  /// FILTER_SPURIOUS    0
  /// ADD_BORDER         0.5
  /// \endcode
  FrameConfig(std::string filename);
  /// Whether objects with <tt>segMap(i) < 0</tt> are replaced by artificial noise,
  /// default = 0.
  static bool FILTER_SPURIOUS;
  /// The amount by which the objects area is enlarged on each side w.r.t. the 
  /// objects area in the segmentation map, default = 0.5 
  /// (50% increase in size an all sides).
  static double ADD_BORDER;
 private:
  char checkKeyword(std::string keyword);
};

#endif
