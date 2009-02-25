#ifndef PROPERTY_H
#define PROPERTY_H

#include <Typedef.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>

/// Allowed types for usage in Property.
#include <boost/variant.hpp>
typedef boost::variant<std::string,data_t,int,std::vector<int>,std::vector<data_t> > variant_t;

/// Flexible container class.
/// The class is derived from \p std::map<std::string, variant_t>,
/// where \p variant_t can be either
/// - \p int
/// - \p data_t
/// - \p std::string
/// - \p std::vector<int>
/// - \p std::vector<data_t>
/// 
/// This allows flexible storage containers, e.g.
/// \code
///  Property p;
///  p["e1"] = 0.1;
///  p["index"] = 1;
///  p["info"] = "Some information...";
///  p["vi"] = std::vector<int>(5,42);
/// \endcode
/// The implementation relies on the concept of \p boost::variant, details can be
/// found at http://www.boost.org/doc/libs/1_38_0/doc/html/variant.html.\n\n
/// In comparisons with or assignment from a \p variant_t, the type has to be known
/// \code
/// std::cout << boost::get<int>(p1["index"]) << std::endl;
/// std::cout << boost::get<data_t>(p1["e1"]) << std::endl;
/// \endcode
/// or checked
/// \code
/// for (Property::iterator iter = p1.begin(); iter != p1.end(); iter++) {
///    int i;
///    if (iter->second.type() == typeid(i))
///      std::cout << iter->first << "\t" << boost::get<int>(iter->second) << std::endl;
///    // ... similar checks for other data types ...
///  }
/// \endcode
/// Alternatively, use the concept of \p boost::static_visitor (described as part
/// of the \p boost::variant documentation). 
class Property : public std::map<std::string, variant_t> {
 public:
  /// Constructor.
  Property();
  /// Write contents of property to \p out.
  /// If the entries should be written with a different separator than \p std::endl,
  /// specify \p linesep. This should only be used in cases where the line separator
  /// needs to be masked in order to be written to a certain format.\n
  /// \b CAUTION: \p linesep must not be either \p "\t" or \p "," as they are
  /// needed for the serialization.
  void write(std::ostream& out, std::string linesep = std::string()) const;
  /// Read contents of property from \p out.
  void read(std::istream& in);
};

#endif
