#ifndef SHAPELENS_KSB_H
#define SHAPELENS_KSB_H

#include "../frame/Object.h"
#include "../frame/Moments.h"

namespace shapelens {
  class KSB {
  public :
    KSB(const Object& obj);
    complex<data_t> chi;
    complex<data_t> gamma() const;
    complex<data_t> gamma(const KSB& psf) const;
    
    data_t trQ, trQ_, M,
      mu_, mu__, psi_, psi__, pi_, pi__, nu_, nu__;
    NumMatrix<data_t> P_sh, P_sm, e_sh, e_sm;
    
  private:
    complex<data_t> __chi(const Moment2& Q) const;
    data_t __trQ(const Moment2& Q) const;
    data_t __psi(const Moment4& Q) const;
    data_t __mu(const Moment4& Q) const;
    data_t __nu(const Moment4& Q) const;
    data_t __pi(const Moment4& Q) const;

    complex<data_t> __p(const KSB& star) const;
  };

} // end namespace shapelens

#endif
