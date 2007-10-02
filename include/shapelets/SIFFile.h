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
#include <bitset>

/// Definition and usage the ShapeletImageFormat.
/// The SIF (Shapelet Image Format) is used to store and retreive all information 
/// that specifies a ShapeletObject.

class SIFFile {
 public:
  /// Argumented constructor.
  SIFFile(std::string filename);

  /// Save the given information to filename
  void save(const ShapeletObject& sobj);
  /// Load the shapelet information from filename
  void load(ShapeletObject& sobj, bool preserve_config);
 private:
  std::string filename;
};

#endif
