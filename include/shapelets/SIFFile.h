#ifndef SIFFILE_H
#define SIFFILE_H

#include <fstream>
#include <string>
#include <iostream>
#include <NumMatrix.h>
#include <Typedef.h>
#include <frame/Grid.h>
#include <frame/Point2D.h>
#include <shapelets/ShapeletObject.h>
#include <shapelets/ShapeletObjectList.h>

namespace shapelens {

/// Read/Write methods for SIF files.
/// The SIF (Shapelet Image Format) is used to store and retrieve all information 
/// that specifies a ShapeletObject.

class SIFFile {
 public:
  /// Argumented constructor.
  SIFFile(std::string filename);
  /// Save the given ShapeletObject to SIFFile
  void save(const ShapeletObject& sobj);
  /// Save the ShapeletObject and its unregularized version to SIFFile.
  void save(const ShapeletObject& sobj, const ShapeletObject& unreg);
  /// Save all ShapeletObject within given ShapeletObjectList to SIFFile.
  /// All ShapeletObject entities will be placed in successive extensions
  /// of the file.
  void save(ShapeletObjectList& sl);
  /// Load the shapelet information from filename
  void load(ShapeletObject& sobj, bool preserve_config);
 private:
  std::string filename;
  void saveSObj(fitsfile* fptr, const ShapeletObject& sobj);
};
} // end namespace
#endif
