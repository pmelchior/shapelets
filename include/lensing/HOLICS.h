#ifndef SHAPELENS_HOLICS_H
#define SHAPELENS_HOLICS_H

#include "../frame/Object.h"
#include "../frame/Moments.h"

namespace shapelens {
  class HOLICS {
  public :
    /// Constructor.
    HOLICS(const Object& obj, data_t sigma);
    complex<data_t> zeta;
    complex<data_t> delta;
    complex<data_t> F() const;
    complex<data_t> G() const;
    complex<data_t> F(const HOLICS& star) const;
    complex<data_t> G(const HOLICS& star) const;
    
  private:
    data_t sigma, 
      trQ, trQ_, trQ__, 
      xi, xi_, xi__, xi___, 
      M, M_,
      v0_, v0___;
    data_t __trQ(const Moment2& Q) const;
    data_t __xi(const Moment4& Q) const;
    data_t __v0(const Moment6& Q) const;
    data_t __C_0_zeta() const;
    data_t __C_Delta_zeta() const;
    data_t __Delta_0_L() const;
    data_t __P_0_zeta() const;
    data_t __P_0_D() const;
    data_t __P_Delta_D() const;
    data_t __P_Delta_zeta() const;
    data_t __C_0_delta() const;
    data_t __P_0_delta() const;
    complex<data_t> __zeta(const Moment3& Q) const;
    complex<data_t> __delta(const Moment3& Q) const;
    complex<data_t> __zeta_q(const HOLICS& star) const;
    complex<data_t> __delta_q(const HOLICS& star) const;
    
    // helper class for Gaussian weight function
    class W {
    public:
      W(data_t sigma, const Point<data_t>& centroid);
      data_t f(const Point<data_t>& P) const;
      data_t f_(const Point<data_t>& P) const;
      data_t f__(const Point<data_t>& P) const;
      data_t f___(const Point<data_t>& P) const;
    private:
      data_t sigma2;
      data_t f(data_t r) const;
      Point<data_t> C;
    };
  };

} // end namespace shapelens

#endif
