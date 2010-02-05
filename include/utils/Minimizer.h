#ifndef SHAPELENS_MINIMIZER_H
#define SHAPELENS_MINIMIZER_H

#include "../Typedef.h"
#include <numla/NumVector.h>

namespace shapelens {
  /// Class for multi-dimensional minimizers.
  class Minimizer {
  public:
    /// Abstract base class for functions to be minimized.
    /// The example below shows the computation of \f$\chi^2\f$ for 
    /// a Gaussian distribution.
    /// \code
    /// class chi2_func : public Minimizer::Functor {
    /// public:
    ///  chi2_func(const NumVector<data_t>& p_, const NumMatrix<data_t>& sigma) : 
    ///    p(p_), sigma_1(sigma.invert()) {
    ///  }
    ///  data_t operator()(const NumVector<data_t>& x) {
    ///    NumVector<data_t> d = x;
    ///    d -= p;
    ///    return d*(sigma_1*d);
    ///  }
    /// private:
    ///  NumVector<data_t> p;
    ///  NumMatrix<data_t> sigma_1;
    /// };
    /// \endcode
    /// Note that pointers or references to \p p and \p sigma would avoid 
    /// copying at construction time, but require these structures to exist
    /// during the minimization.
  class Functor {
  public:
    /// Virtual function evaluation.
    virtual data_t operator () (const NumVector<data_t>&) = 0;
  };
  
  /// Minimize func with Powell's algorithm.
  /// \p p is the initial guess of the minimum of \p func and updated by
  /// the method. The return value is the found minimum of \p func within the
  /// uncertainty of \p ftol.\p itmax is the maximum number of iterations
  /// (for both the number of line searches and the calls during each
  /// line search).\n\n
  /// Employs the modified algorithm from Numerical Recipes (Press et al., 1992)
  static data_t Powell(Functor& func, NumVector<data_t>& p, data_t ftol, unsigned int itmax = 100);
  static data_t Simplex(Functor& func, NumVector<data_t>& p, data_t ftol, unsigned int itmax = 100);
  };
} // end namespace


#endif
