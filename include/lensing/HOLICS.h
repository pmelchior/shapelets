#ifndef SHAPELENS_HOLICS_H
#define SHAPELENS_HOLICS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"

namespace shapelens {
  class HOLICS {
  public :
    HOLICS(const Object& obj);
    complex<data_t> zeta;
    complex<data_t> delta;
    complex<data_t> F() const;
    complex<data_t> G() const;
    complex<data_t> F(const HOLICS& psf) const;
    complex<data_t> G(const HOLICS& psf) const;
    
    data_t sigma, 
      trQ, trQ_, trQ__, 
      xi, xi_, xi__, xi___, 
      M, M_,
      v0_, v0___,
      C_0_zeta, C_Delta_zeta, C_0_delta,
      Delta_0_L,
      P_0_zeta, P_0_D, P_Delta_D, P_Delta_zeta, P_0_delta;
    
  private:
    ...
    data_t __trQ(const Moment2& Q) const;
    data_t __xi(const Moment4& Q) const;
    data_t __v0(const Moment6& Q) const;
    complex<data_t> __zeta(const Moment3& Q) const;
    complex<data_t> __delta(const Moment3& Q) const;
    complex<data_t> __zeta_q(const HOLICS& star) const;
    complex<data_t> __delta_q(const HOLICS& star) const;
  };

} // end namespace shapelens

#endif
