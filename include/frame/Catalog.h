#ifndef SHAPELENS_CATALOG_H
#define SHAPELENS_CATALOG_H

#include <map>
#include <string>
#include <bitset>
#include <list>
#include "../Typedef.h"
#include "CoordinateTransformation.h"
#include "../utils/IO.h"
#include "../utils/Property.h"

namespace shapelens {

/// Structure to store information of Catalog objects.
/// The names of the parameters are indentical or similar to those used by 
/// SExtractor. Coordinates are supposed to be given in image units, 
/// starting at (1/1) in the left-lower corner.
///
/// Allowed values for <tt>PARAMETER_NAME</tt> in the file read by Catalog 
/// are listed below.
/// In addition to these, the number of the object must be given; 
/// its <tt>PARAMETER_NAME</tt> is one these: <tt>ID,NR,NUMBER</tt>.

struct CatObject {
  /// Minimum <tt>X</tt> coordinate of cutout (allowed: <tt>XMIN,XMIN_IMAGE</tt>).
  int XMIN;
  /// Maximum <tt>X</tt> coordinate of cutout (allowed: <tt>XMAX,XMAX_IMAGE</tt>).
  int XMAX;
  /// Minimum <tt>Y</tt> coordinate of cutout (allowed: <tt>YMIN,YMIN_IMAGE</tt>).
  int YMIN;
  /// Maximum <tt>Y</tt> coordinate of cutout (allowed: <tt>YMAX,YMAX_IMAGE</tt>).
  int YMAX;
  /// <tt>X</tt> coordinate of the centroid (alternatives: <tt>X,X_IMAGE,XWIN,XWIN_IMAGE,XPEAK,XCENTROID</tt>).
  data_t XCENTROID;
  /// <tt>Y</tt> coordinate of the centroid (alternatives: <tt>Y,Y_IMAGE,YWIN,YWIN-IMAGE,YPEAK,YCENTROID</tt>).
  data_t YCENTROID;
  /// Flags set by a previous segmentation software (alternatives: none).
  unsigned char FLAGS;
  /// Storage for optional values.
  Property OPT;
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
/// where \p COL_NR is an integer > 1 and \p PARAMETER_NAME is a
/// string. After that follows the data section. All columns are separated by
/// one or multiple blank or tabulator characters.\n\n
/// See CatObject for details on mandatory and optional values of 
/// \p PARAMETER_NAME.
class Catalog : public std::map<unsigned long, CatObject> {
 public:
  /// Default constructor.
  Catalog();
  /// Argumented contructor to read Catalog from <tt>catfile</tt>.
  /// The file can either be given in ASCII and conform with the SExtractor 
  /// format, or be a FITS table.\n
  /// \p optional denotes a list of optional columns in \p catfile, whose
  /// values are inserted into CatObject::OPT.
  Catalog(std::string catfile, const std::list<std::string>& optional = std::list<std::string>());
  /// Read <tt>catfile</tt>.
  /// The file must be given in ASCII and conform with the SExtractor format.
  void read(std::string catfile, const std::list<std::string>& optional = std::list<std::string>());
  /// Save Catalog to file.
  /// The file will be stored in ASCII and conform with the SExtractor format.
  void save(std::string catfile) const;
  /// Apply coordinate transformation to all entries.
  void apply(const CoordinateTransformation& C);
  /// Add a catalog to another.
  /// Entries are identified by map index. Thus, if an entry from \p *this has the same index as
  /// one from \p c, it will be overwritten.
  Catalog operator+ (const Catalog& c);
  /// Add catalog to \p *this.
  /// Behaves identical to operator+.
  void operator+= (const Catalog& c);
  /// Substracts a catalog from another.
  /// Entries are identified by map index.
  Catalog operator- (const Catalog& c);
  /// Substracts a catalog from \p *this.
  /// Entries are identified by map index.
  void operator-= (const Catalog& c);
  /// Match a catalog with another.
  /// The resulting Catalog will only have entries which are in both catalogs; its
  /// indices are the ones from \p *this.
  /// Matches are identified by the spatial proximity of their centroids:
  /// Using a R-tree, objects from \p c are found whose bounding boxes 
  ////(<tt>XMIN..XMAX,YMIN..YMAX</tt>) overlap with the one from objects of \p *this.\n
  /// In case of multiple matches (which can happend when objects are blended in one catalog
  /// but split in the other), a combined centroid 
  /// \f[\vec{x}_c = \frac{\sum_\text{matches} F_i\ \vec{x}_i}{F_c}\f] is computed, where 
  /// \f$F_i\f$ and \f$\vec{x}_i\f$ is the flux and centroid of object \f$i\f$.\n
  /// If \f$\vec{x}_c\f$ matches better than any single object's centroid, all matches are mapped
  /// to the correspondent object in the other catalog; if not, the nearest neighbor is mapped only\n
  //Catalog operator*(const Catalog& c);
  /// Match a catalog with \p *this.
  /// Behaves identical to operator*.
  //void operator*= (const Catalog& c);
  /// Remove a catalog from another.
  /// The resulting Catalog will only have entries which are in \p *this, but not in \p c;
  /// its indices are the ones from \p *this.\n
  /// Matching is done as in operator*.
  //Catalog operator/(const Catalog& c);
  /// Remove a catalog from \p *this.
  /// Behaves identical to operator/.
  //void operator/=(const Catalog& c);
  
 private:
  struct OptFormat {
    unsigned short COLNR;
    int TYPE;
  };
  struct CatFormat {
    unsigned short ID;
    unsigned short XMIN;
    unsigned short XMAX;
    unsigned short YMIN;
    unsigned short YMAX;
    unsigned short XCENTROID;
    unsigned short YCENTROID;
    unsigned short FLAGS;
    std::map<std::string, OptFormat> OPT;
  };
  CatFormat format;
  std::bitset<8> present;
  bool formatChecked;
  bool checkFormat();
  void setFormatField(std::string name, unsigned short colnr, const std::list<std::string>& optional);
  void setFormatFromFITSTable(fitsfile* fptr, const std::list<std::string>& optional);
  void setOptionalField(fitsfile* fptr, int row, const OptFormat& of, const std::string& name, CatObject& co);
  int round(data_t x);
};
} // end namespace
#endif
