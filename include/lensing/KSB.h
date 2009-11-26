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
    complex<data_t> gamma_nl(data_t accuracy = 0.001) const;
    complex<data_t> gamma_nl(const KSB& psf, data_t accuracy = 0.001) const;

  private:
    complex<data_t> __chi(const Moment2& Q) const;
    data_t __trQ(const Moment2& Q) const;
    data_t __psi(const Moment4& Q) const;
    data_t __mu(const Moment4& Q) const;
    data_t __nu(const Moment4& Q) const;
    data_t __pi(const Moment4& Q) const;
    data_t __lambda(const Moment6& Q) const;
    data_t __omega(const Moment6& Q) const;
    data_t __sigma(const Moment6& Q) const;
    data_t __rho(const Moment6& Q) const;
    data_t __ix(const Moment6& Q) const;
    data_t __delta(const Moment6& Q) const;
    
    complex<data_t> __p(const KSB& star) const;

  public:
    data_t trQ, trQ_, M,
      mu_, mu__, psi_, psi__, pi_, pi__, nu_, nu__, lambda__, omega__, sigma__,rho__,ix__,delta__;
    NumMatrix<data_t> P_sh, P_sm, e_sh, e_sm;

    // helper class for tensor of rank 3:
    // glue two NumMatrices togehter...
    class NumTensor {
    public:
      NumTensor();
      data_t& operator()(bool i, bool j, bool k);
      const data_t& operator()(bool i, bool j, bool k) const;
    private:
      NumMatrix<data_t> M0, M1;
    };
    NumTensor R;
    
  };

} // end namespace shapelens

#endif
