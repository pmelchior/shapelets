#ifndef SHAPELENS_KSB_H
#define SHAPELENS_KSB_H

#include "../frame/Object.h"
#include "../frame/Moments.h"

namespace shapelens {
  class KSB {
  public :
    /// Constructor.
    KSB(const Object& obj);
    complex<data_t> chi;
    complex<data_t> gamma() const;
    complex<data_t> gamma(const KSB& psf) const;
    
  private:
    data_t trQ, trQ_, M,
      a_, a__, b_, b__, c_, c__, d_, e_, f_, f__, g_, g__;
    NumMatrix<data_t> P_sh, P_sm;
    
    complex<data_t> __chi(const Moment2& Q) const;
    data_t __trQ(const Moment2& Q) const;
    data_t __a(const Moment2& Q2, const Moment4& Q4) const;
    data_t __b(const Moment1& Q1, const Moment3& Q3) const;
    data_t __c(const Moment2& Q) const;
    data_t __d(const Moment2& Q) const;
    data_t __e(const Moment1& Q) const;
    data_t __f(const Moment1& Q1, const Moment3& Q3) const;
    data_t __g(const Moment4& Q) const;

    complex<data_t> __p(const KSB& star) const;
  };

} // end namespace shapelens

#endif
