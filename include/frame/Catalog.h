#ifndef CATALOG_H
#define CATALOG_H

#include <map>
#include <string>
#include <bitset>
#include <Typedef.h>
//#include <rtree.hh>

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
  /// Optional parent id (alternative: <tt>VECTOR_ASSOC</tt>).
  unsigned long PARENT;
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
    unsigned short PARENT;
  };
  CatFormat format;
  std::bitset<11> present;
  bool formatChecked;
  bool checkFormat();
  void setFormatField(std::string name, unsigned short colnr);
  template <class T>
  class Rectangle {
  public: 
    Rectangle() {
      min[0] = max[0] = min[0] = max[0] = T(0);
    }
    Rectangle (const T& xmin, const T& xmax, const T& ymin, const T& ymax) {
      setCoords(xmin,xmax,ymin,ymax);
    }
    void setCoords(const T& xmin, const T& xmax, const T& ymin, const T& ymax) {
      min[0] = xmin;
      min[1] = ymin;
      max[0] = xmax;
      max[1] = ymax;
    }
    void setCoords(Catalog::const_iterator& iter) {
      min[0] = iter->second.XMIN;
      min[1] = iter->second.YMIN;
      max[0] = iter->second.XMAX;
      max[1] = iter->second.YMAX;
    }
    T* getMin() {
      return min;
    }
    T* getMax() {
      return max;
    }
  private:
    T min[2];
    T max[2];
  };
  //void buildRTree(RTree<unsigned long, unsigned long, 2, data_t>& rtree, const Catalog& c) {
  //  Catalog::Rectangle<unsigned long> rect;
  //  for (Catalog::const_iterator iter = c.begin(); iter != c.end(); iter++) {
  //    rect.setCoords(iter->second.XMIN,iter->second.XMAX,iter->second.YMIN,iter->second.YMAX);
  //    rtree.Insert(rect.getMin(),rect.getMax(),iter->first);
  //  }
  //}
};

#endif
