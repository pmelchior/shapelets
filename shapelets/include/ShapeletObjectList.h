#ifndef SHAPELETOBJECTLIST_H
#define SHAPELETOBJECTLIST_H

#include <ShapeletObject.h>
#include <list>
#include <string>

/// Class for managing an ensemble of ShapeletObjects.
/// This class is meant for loading a set of ShapeletObject entities by giving a list
/// of SIFFile names as an ASCII file.\n
/// Apart from that, one can select an arbitrary criterium for any ShapeletObject 
/// to select a subset of these objects by giving an appropriate
/// <tt>selectionFunction</tt>:
/// \code
/// bool selectionFunction(ShapeletObject& sobj) {
///   if (sobj.getDecompositionChiSquare() < 1) 
///     return 1;
///   else
///     return 0;
/// }
/// ShapeletObjectList sl("file.lst",&selectionFunction);
/// \endcode
/// In this example, all SIFFile names stored in <tt>file.lst</tt> are opened;
/// if their decomposition \f$\chi^2\f$ is less than 1, they will be included in
/// the list.


class ShapeletObjectList : public std::list<ShapeletObject*> {
 public:
  /// Constructor for loading all files in <tt>listfile</tt>.
  ShapeletObjectList(std::string listfile);
  /// Constructor for loading a subset of files from <tt>listfile</tt>. 
  /// The ShapeletObject <tt>sobj</tt> entities have to fulfill the criterium 
  /// <tt>selectionFunction (sobj) == 1</tt> to be included in the list.
  ShapeletObjectList(std::string listfile, bool (* selectionFunction) (ShapeletObject&));
  /// Overloaded destructor.
  ~ShapeletObjectList();
 private:
  void readListFile(std::string listfile, bool (* selectionFunction) (ShapeletObject&));
};
#endif
