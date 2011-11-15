#ifndef SHAPELENS_INTERPOLATOR
#define SHAPELENE_INTERPOLATOR

#include "../Typedef.h"
#include "../frame/Point.h"
#include "../frame/WeightFunction.h"
#include "../utils/RTree.h"
#include <stdexcept>

namespace shapelens {
  /// Base class for 2D interpolation of arbitrarily valued data.
  template <class T>
    class Interpolator {
  public:
    /// Constructor to specify \p values at \p points
  Interpolator(const std::vector<Point<data_t> >& points, const std::vector<T>& values) : points(points), values(values) {}
    virtual ~Interpolator();
    /// Interpolates \p values at position \p P.
    virtual T operator()(const Point<data_t>& P) const = 0;
  protected:
    const std::vector<Point<data_t> >& points;
    const std::vector<T>& values;
  };
  
  template <class T>
    inline Interpolator<T>::~Interpolator() {}

  /// Nearest Neighbor Interpolator.
  /// Uses the (weighted) average of the \p k nearest neighbors as interpolator.
  template <class T>
  class NearestNeighborInterpolator : public Interpolator<T> {
  public:
    /// Constructor.
    /// Specifies \p  \p values at \p points. \p k denotes the number of nearest
    /// neighbors to consider, \p w their relative weight.\n
    /// \b NOTE: Since \p w is given at construction time, it needs to
    /// allow modification (e.g. recentering to the point of interest).
  NearestNeighborInterpolator(const std::vector<Point<data_t> >& points, const std::vector<T>& values, unsigned int k, LocalWeightFunction& w) : 
    Interpolator<T>(points, values), w(w), k(k) {
      tree.insertNodes(points);
      if (k < 1)
	throw std::invalid_argument("NearestNeighborInterpolator: k must be >= 1!");
	}
    
    /// Interpolates \p values at position \p P.
    virtual T operator()(const Point<data_t>& P) const {
      w.setCentroid(P);
      std::list<unsigned long> neighbors = tree.getNearestNeighbors(k, P);
      // initialize mean and sum of weights
      std::list<unsigned long>::iterator iter = neighbors.begin();
      data_t sum_w = w(Interpolator<T>::points[*iter]);
      T mean = Interpolator<T>::values[*iter];
      mean *= sum_w;
      iter++;
      for (iter; iter != neighbors.end(); iter++) {
	data_t weight = w(Interpolator<T>::points[*iter]);
	T tmp = Interpolator<T>::values[*iter];
	tmp *= weight;
	mean += tmp;
	sum_w += weight;
      }
      mean /= sum_w;
      return mean;
    }
  private:
    LocalWeightFunction& w;
    RTree tree;
    unsigned int k;
  };
} // end namespace

#endif
