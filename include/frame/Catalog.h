#ifndef CATALOG_H
#define CATALOG_H

#include <map>
#include <string>
#include <bitset>
#include <Typedef.h>

/// Structure to store information of Catalog objects.
/// The names of the parameters are indentical or similar to those used by SExtractor.
///
/// Allowed values for <tt>PARAMETER_NAME</tt> in the file read by Catalog are listed below.
/// In addition to these, the number of the object must be given; 
/// its <tt>PARAMETER_NAME</tt> is one these: <tt>ID,NR,NUMBER</tt>.


struct CatObject {
  /// Minimum <tt>X</tt> coordinate of cutout (allowed: <tt>XMIN*</tt>).
  int XMIN;
  /// Maximum <tt>X</tt> coordinate of cutout (allowed: <tt>XMAX*</tt>).
  int XMAX;
  /// Minimum <tt>Y</tt> coordinate of cutout (allowed: <tt>YMIN*</tt>).
  int YMIN;
  /// Maximum <tt>Y</tt> coordinate of cutout (allowed: <tt>YMAX*</tt>).
  int YMAX;
  /// <tt>X</tt> coordinate of the centroid (alternatives: <tt>X,X_*,XWIN*,XPEAK*</tt>).
  data_t XCENTROID;
  /// <tt>Y</tt> coordinate of the centroid (alternatives: <tt>Y,Y_*,YWIN*,YPEAK*</tt>).
  data_t YCENTROID;
  /// Flux of the Object (allowed: <tt>FLUX*</tt>).
  data_t FLUX;
  /// Flags set by a previous segmentation software (alternatives: none).
  unsigned char FLAGS;
  /// Optional classifier of the Object (alternatives: <tt>CLASS_STAR</tt>).
  data_t CLASSIFIER;
};

/// Class for reading and storing catalogs of objects.
/// A Catalog is a <tt>std::map<unsigned long,CatObject></tt>. That means, 
/// it takes a number and returns a CatObject entity, which lists several
/// properties of this object:
/// \code
/// Catalog cat(catalogfile);
/// CatObject ca = cat[1];
/// std::cout << "FLUX = " << ca.FLUX << std::endl;
/// \endcode
/// 
/// The catalog file will be parsed with the assumption that it starts with 
/// a SExtractor-like header of the form
/// \code
/// # COL_NR PARAMETER_NAME
/// \endcode
/// where <tt>COL_NR</tt> is an integer > 1 and <tt>PARAMETER_NAME</tt> is a string.
/// After that follows the data section. All columns are separated by one or multiple
/// blank characters.\n\n
/// See CatObject for details on mandatory and optional values of <tt>PARAMETER_NAME</tt>.


class Catalog : public std::map<unsigned long, CatObject> {
 public:
  /// Default constructor.
  Catalog();
  /// Argumented contructor to read Catalog from <tt>catfile</tt>.
  /// The file must be given in ASCII and conform with the SExtractor format.
  Catalog(std::string catfile);
  /// Read <tt>catfile</tt>.
  /// The file must be given in ASCII and conform with the SExtractor format.
  void read(std::string catfile);
  /// Save Catalog to file.
  /// The file will be stored in ASCII and conform with the SExtractor format.
  void save(std::string catfile) const;
 private:
  struct CatFormat {
    unsigned short ID;
    unsigned short XMIN;
    unsigned short XMAX;
    unsigned short YMIN;
    unsigned short YMAX;
    unsigned short XCENTROID;
    unsigned short YCENTROID;
    unsigned short FLUX;
    unsigned short FLAGS;
    unsigned short CLASSIFIER;
  };
  CatFormat format;
  std::bitset<10> present;
  bool formatChecked;
  bool checkFormat();
  void setFormatField(std::string name, unsigned short colnr);
};

#endif
