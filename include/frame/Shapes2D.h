#ifndef SHAPELENS_SHAPES2D_H
#define SHAPELENS_SHAPES2D_H

#include <Typedef.h>
#include <frame/Point2D.h>
#include <list>
#include <stdexcept>

namespace shapelens {
  /// Rectangular patch.
  template <class T>
    class Rectangle {
  public:
    /// Lower-left boundary point.
    Point2D<T> ll;
    /// Top-right boundary point.
    Point2D<T> tr;
  };
  
  /// Edge of Polygon.
  template <class T>
    class Edge {
  public:
    /// First point of edge.
    Point2D<T> p1;
    /// Second point of edge.
    Point2D<T> p2;
    /// Move edge to next point \p p.
    /// Replace p1 by p2 and sets p2 to \p p.
    void moveTo(const Point2D<T>& p) {
      p1 = p2;
      p2 = p;
    }
  };

  /// Class for simple polygons.
  /// The class does not check whether the polygon is in fact simple or complex,
  /// (a complex polygon is one with crossing edges);
  /// this has to be done by the user when constructing it and can be checked
  /// with checkEdges().
  ///
  /// \b CAUTION: For integer-typed polygons, edge-artefacts can occur
  /// when testing isInside().
  template <class T>
  class Polygon {
  public:
    /// Constructor.
    Polygon() {};
    /// Add edge to polygon.
    /// The starting point of \p e must be identical to the end-point of the
    /// prior edge.
    void addEdge(const Edge<T>& e) {
      if (edges.size() == 0)
	edges.push_back(e);
      else if (e.p1 == edges.back().p2)
	edges.push_back(e);
      else
	throw std::invalid_argument("Polygon: Edges do not match");
    }
    data_t getArea() const {
      data_t a = 0;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	a += (iter->p1(0)*iter->p2(1)) - (iter->p2(0)*iter->p1(1));
      }
      return 0.5*a;
    }
    /// Checks whether \p p is inside the polygon.
    /// Uses crossing test: If a ray originating from \p into positive y-direction
    /// crosses an odd number of edges, its inside the polygon.
    bool isInside(const Point2D<T>& p) const {
      unsigned int crossings = 0;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	// x-cordinate of edge-points above and below p
	if ((iter->p1(0) - p(0))*(iter->p2(0) - p(0)) <= 0) {
	  // both y-coordinates above: must cross
	  if ((iter->p1(1) > p(1)) && (iter->p2(1) > p(1)))
	    crossings++;
	  // y-coords differ: crosses only if p is below edge line
	  else if ((iter->p1(1) - p(1))*(iter->p2(1) - p(1)) <= 0) {
	    data_t m = (iter->p2(1) - iter->p1(1))/(iter->p2(0) - iter->p1(0));
	    data_t dy = m*(p(0) - iter->p1(0));
	    if (iter->p1(1) + dy > p(1))
	      crossings++;
	  }
	}
      }
      return bool(crossings%2);
    }
    /// Return Rectangle which bounds the polygon.
    Rectangle<T> getBoundingRectangle() const {
      shapelens::Point2D<double> min,max;
      T mx,my,MX,MY;
      for (typename std::list<Edge<T> >::const_iterator iter = edges.begin(); iter != edges.end(); iter++) {
	mx = std::min(iter->p1(0),iter->p2(0));
	MX = std::max(iter->p1(0),iter->p2(0));
	my = std::min(iter->p1(1),iter->p2(1));
	MY = std::max(iter->p1(1),iter->p2(1));
	if (iter == edges.begin()) {
	  min(0) = mx;
	  min(1) = my;
	  max(0) = MX;
	  max(1) = MY;
	} else {
	  if (mx < min(0))
	    min(0) = mx;
	  if (my < min(1))
	    min(1) = my;
	  if (MX > max(0))
	    max(0) = MX;
	  if (MY > max(1))
	    max(1) = MY;
	}
      }
      shapelens::Rectangle<T> r;
      r.ll = min;
      r.tr = max;
      return r;
    }
    /// Checks whether polygonial chain is closed.
    /// Should be called after completion of the chain to ensure that first
    /// and last edge-points are identical.
    bool checkEdges() {
      if (edges.front().p1 == edges.back().p2)
	return true;
      else
	return false;
    }
    /// Clear all edges of this polygon.
    void clear() {
      edges.clear();
    }
  private:
    std::list<Edge<T> > edges;
  };

  
} // end namespace

#endif
