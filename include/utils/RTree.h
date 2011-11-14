#ifndef SHAPELENS_RTREE_H
#define SHAPELENS_RTREE_H

#ifdef HAS_SPATIALINDEX
#include <spatialindex/SpatialIndex.h>
#include "../Typedef.h"
#include "../frame/Point.h"
#include "../frame/Shapes.h"
#include <vector>
#include <list>

namespace shapelens {
  /// Class for efficient 2D lookups.
  /// Uses a R* tree index from http://libspatialindex.github.com/.
  class RTree {
  public:
    /// Default constructor.
    /// \p indexCapacity and \p leafCapacity denotes the number of nodes
    /// the index or the leaf nodes can hold before it needs to be split.
    RTree(unsigned int indexCapacity = 20, unsigned int leafCapacity = 20);
    /// Insert nodes.
    /// Nodes are Rectangle entities for which the RTree stores their index.\n
    /// \b CAUTION: Don't mix with Point nodes.
    void insertNodes(const std::vector<shapelens::Rectangle<data_t> >& vr);
    /// Insert nodes.
    /// Nodes are Point entities for which the RTree stores their index.
    /// \b CAUTION: Don't mix with Rectangle nodes.
    void insertNodes(const std::vector<shapelens::Point<data_t> >& vr);
    /// Destructor.
    ~RTree();
    /// Find object whose support Rectangle overlaps with \p P.
    /// The list contains the vector indices of those objects whose rectangles
    /// were given at construction time.\n
    /// \b CAUTION: Use this for Rectangle nodes only
    std::list<unsigned long> getMatches(const shapelens::Point<data_t>& P) const;
    /// Find \p k nearest neighbors (assuming Euclidean metric).
    /// The list contains the vector indices of those objects whose points
    /// were given at construction time.\n
    /// \b CAUTION: Use this for Point nodes only.
    std::list<unsigned long> getNearestNeighbors(unsigned int k, const shapelens::Point<data_t>& P) const;

  private:
    SpatialIndex::ISpatialIndex* tree;
    SpatialIndex::IStorageManager* mem;
  };
} // end namespace

#endif // HAS_SPATIALINDEX
#endif
