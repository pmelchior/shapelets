#include "../../include/lensing/KSB.h"
#include "../../include/utils/Minimizer.h"

namespace shapelens {
  // helper function
  complex<data_t> operator*(const NumMatrix<data_t>& M, complex<data_t> d) {
    return complex<data_t> (M(0,0)*real(d)+M(0,1)*imag(d),M(1,0)*real(d)+M(1,1)*imag(d));
  }
  // helper class
  KSB::NumTensor::NumTensor() : M0(2,2), M1(2,2) {
  }
  data_t& KSB::NumTensor::operator()(bool i, bool j, bool k) {
    if (i)
      return M1(j,k);
    else
      return M0(j,k);
  }
  const data_t& KSB::NumTensor::operator()(bool i, bool j, bool k) const {
    if (i)
      return M1(j,k);
    else
      return M0(j,k);
  }

  KSB::KSB(const Object& obj_) {
    // for changing obj.w, we need to have write permissions...
    Object& obj = const_cast<Object&>(obj_);
    Moment0 M0(obj);
    M = M0(0);
    if (obj.w.getType() == 0)
      throw std::invalid_argument("KSB: usage of flat weights is not permitted");
    Moment2 Q2(obj);
    trQ = __trQ(Q2);
    chi = __chi(Q2);

    // measure moments with W' (w.r.t r^2)
    obj.w.setDerivative(-1);
    Moment2 Q2_(obj);
    Moment4 Q4_(obj);
    trQ_ = __trQ(Q2_);
    psi_ = __psi(Q4_);
    mu_ = __mu(Q4_);
    pi_ = __pi(Q4_);
    nu_ = __nu(Q4_);

    // measure moments with W'' (w.r.t r^2)
    obj.w.setDerivative(-2);
    Moment4 Q4__(obj);
    mu__ = __mu(Q4__);
    psi__ = __psi(Q4__);
    pi__ = __pi(Q4__);
    nu__ = __nu(Q4__);

    Moment6 Q6__(obj);
    lambda__ = __lambda(Q6__);
    omega__  = __omega(Q6__);
    sigma__  = __sigma(Q6__);
    rho__    = __rho(Q6__);
    ix__     = __ix(Q6__);
    delta__  = __delta(Q6__);
    
    // measure moments with W''' (w.r.t r^2)
    obj.w.setDerivative(-3);
    Moment8 Q8___(obj);
    a1___=___a1(Q8___);
    a2___=___a2(Q8___);
    a3___=___a3(Q8___);
    a4___=___a4(Q8___);
    a5___=___a5(Q8___);

    // reset to zero-th derivative since obj is const
    obj.w.setDerivative(0);


    R(0,0,0) = (real(chi)/trQ)*(2*rho__ + 4*psi_) - 2*lambda__/trQ - 4*pi_/trQ;
    R(0,0,1) = (real(chi)/trQ)*(4*ix__  + 8*mu_ ) - 4*omega__/trQ;
    R(0,1,0) = R(0,0,1) - 8*nu_/trQ;
    R(0,1,1) = (real(chi)/trQ)*(8*delta__ + 16*Q4_(0,0,1,1)) - 8*sigma__/trQ;
    
    R(1,0,0) = (imag(chi)/trQ)*(2*rho__ + 4*psi_) - 4*omega__/trQ;
    R(1,0,1) = (imag(chi)/trQ)*(4*ix__ + 8*mu_) - 8*sigma__/trQ - 4*pi_/trQ;
    R(1,1,0) = (imag(chi)/trQ)*(4*ix__ + 8*mu_) - 8*sigma__/trQ; 
    R(1,1,1) = (imag(chi)/trQ)*(8*delta__ + 16*Q4_(0,0,1,1)) - 16*Q6__(0,0,0,1,1,1)/trQ - 8*nu_/trQ;
    

    R1(0,0,0) =  - 2*lambda__/trQ - 4*pi_/trQ;
    R1(0,0,1) =  - 4*omega__/trQ;
    R1(0,1,0) = R(0,0,1) - 8*nu_/trQ;
    R1(0,1,1) =  - 8*sigma__/trQ;
    
    R1(1,0,0) =  - 4*omega__/trQ;
    R1(1,0,1) =  - 8*sigma__/trQ - 4*pi_/trQ;
    R1(1,1,0) =  - 8*sigma__/trQ; 
    R1(1,1,1) =  - 16*Q6__(0,0,0,1,1,1)/trQ - 8*nu_/trQ;

    // Compute U (involved 8-th moments and III derivative of
    // W. Responsible for terms g*g*g)
    U1(0,0,0)=a1___/trQ;
    U2(0,0,0)=a5___/trQ;
    U1(1,0,0)=U1(0,1,0)=U1(0,0,1)=U2(0,0,0);
    U2(1,0,0)=a4___/trQ;
    U2(0,0,1)=U2(0,1,0)=U1(0,1,1)=U1(1,1,0)=U1(1,0,1)=U2(1,0,0);
    U2(1,1,0)=a3___/trQ;
    U2(0,1,1)=U2(1,0,1)=U1(1,1,1)=U2(1,1,0);
    U2(1,1,1)=a2___/trQ;
   
    // Compute P_sh 
    Pa.resize(2,2);
    Pa(0,0) = -2*real(chi)*real(chi);//+2*psi_/trQ + 2.0 ;
    Pa(0,1) = -2*imag(chi)*real(chi);//+4*mu_/trQ;
    Pa(1,0) = -2*imag(chi)*real(chi);//+4*mu_/trQ;
    Pa(1,1) = -2*imag(chi)*imag(chi);//+8*Q4_(0,0,1,1)/trQ+2.0;

    Pb.resize(2,2);
    Pb(0,0) = 2*psi_/trQ;;
    Pb(0,1) = +4*mu_/trQ;
    Pb(1,0) = 4*mu_/trQ;
    Pb(1,1) = 8*Q4_(0,0,1,1)/trQ;

    Pc.resize(2,2);
    Pc(0,0) =  -2*real(chi)*pi_/trQ;
    Pc(0,1) = -4*real(chi)*nu_/trQ;
    Pc(1,0) = -2*imag(chi)*pi_/trQ;
    Pc(1,1) =  -4*imag(chi)*nu_/trQ;
    
    P_sh.resize(2,2);
    P_sh(0,0) = -2*real(chi)*pi_/trQ-2*real(chi)*real(chi)+2*psi_/trQ + 2.0 ;
    P_sh(0,1) = -4*real(chi)*nu_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,0) = -2*imag(chi)*pi_/trQ-2*imag(chi)*real(chi)+4*mu_/trQ;
    P_sh(1,1) = -4*imag(chi)*nu_/trQ-2*imag(chi)*imag(chi)+8*Q4_(0,0,1,1)/trQ+2.0;
    
    //Compute K
    K.resize(2,2);
    K(0,0)=rho__;
    K(0,1)=2*ix__;
    K(1,0)=K(0,1);
    K(1,1)=4*delta__;

    //Compute B
    B.resize(2,2);
    B(0,0)=psi_;
    B(0,1)=2*mu_;
    B(1,0)=B(0,1);
    B(1,1)=4*Q4_(0,0,1,1);


    //Compute P_sh (stricktly I order)
    P1.resize(2,2);
    P1(0,0) = 2*B(0,0)/trQ + 2. ;
    P1(0,1) = 2*B(0,1)/trQ;
    P1(1,0) = P1(0,1);
    P1(1,1) = 2*B(1,1)/trQ + 2. ;
        
    //Compute P_sh (stricktly II order)
    P2.resize(2,2);
    P2(0,0) = -2*real(chi)*pi_/trQ;
    P2(0,1) = -4*real(chi)*nu_/trQ;
    P2(1,0) = -2*imag(chi)*pi_/trQ;
    P2(1,1) = -4*imag(chi)*nu_/trQ;
    
    
    // compute P_sm
    P_sm.resize(2,2);
    P_sm(0,0) = M + 2*trQ_ + psi__;
    P_sm(0,1) = 2*mu__;
    P_sm(1,0) = P_sm(0,1);
    P_sm(1,1) = M + 2*trQ_ + 4*Q4__(0,0,1,1);
    P_sm /= trQ;
    e_sm.resize(1,2);
    e_sm(0,0) = (2*(Q2_(0,0)-Q2_(1,1)) + pi__)/trQ;
    e_sm(0,1) = (4*Q2_(0,1) + 2*nu__)/trQ;
    P_sm(0,0) -= real(chi)*e_sm(0,0);
    P_sm(0,1) -= real(chi)*e_sm(0,1);
    P_sm(1,0) -= imag(chi)*e_sm(0,0);
    P_sm(1,1) -= imag(chi)*e_sm(0,1);
  }

  complex<data_t> KSB::__chi(const Moment2& Q) const {
    return complex<data_t>(Q(0,0)-Q(1,1),2*Q(0,1))/trQ;
  }
  data_t KSB::__trQ(const Moment2& Q) const {
    return Q(0,0) + Q(1,1);
  }
  data_t KSB::__psi(const Moment4& Q) const {
    return Q(0,0,0,0) - 2*Q(0,0,1,1) + Q(1,1,1,1);
  }
  data_t KSB::__mu(const Moment4& Q) const {
    return Q(0,0,0,1) - Q(0,1,1,1);
  }
  data_t KSB::__nu(const Moment4& Q) const {
    return Q(0,0,0,1) + Q(0,1,1,1);
  }
  data_t KSB::__pi(const Moment4& Q) const {
    return Q(0,0,0,0) - Q(1,1,1,1);
  }
  data_t KSB::__lambda(const Moment6& Q) const {
    return Q(0,0,0,0,0,0) - 3*Q(0,0,0,0,1,1) + 3*Q(0,0,1,1,1,1) - Q(1,1,1,1,1,1);
  }
  data_t KSB::__omega(const Moment6& Q) const {
    return Q(0,0,0,0,0,1) + Q(0,1,1,1,1,1) - 2*Q(0,0,0,1,1,1);
  }
  data_t KSB::__sigma(const Moment6& Q) const {
    return Q(0,0,0,0,1,1) - Q(0,0,1,1,1,1);
  }
  data_t KSB::__rho(const Moment6& Q) const {
    return Q(0,0,0,0,0,0) - Q(0,0,0,0,1,1) - Q(0,0,1,1,1,1) + Q(1,1,1,1,1,1);
  }
  data_t KSB::__ix(const Moment6& Q) const {
    return Q(0,0,0,0,0,1) - Q(0,1,1,1,1,1);
  }
  data_t KSB::__delta(const Moment6& Q) const {
    return Q(0,0,0,0,1,1) + Q(0,0,1,1,1,1);
  }

  data_t KSB::___a1(const Moment8& Q) const {
    return Q(0,0,0,0,0,0,0,0)-4.*Q(0,0,0,0,0,0,1,1)+6.*Q(0,0,0,0,1,1,1,1)-4.*Q(0,0,1,1,1,1,1,1)+Q(1,1,1,1,1,1,1,1);
  }
  
  data_t KSB::___a2(const Moment8& Q) const {
    return 16.*Q(0,0,0,0,1,1,1,1);
  }
  data_t KSB::___a3(const Moment8& Q) const {
   return 8.*(Q(0,0,0,0,0,1,1,1)-Q(0,0,0,1,1,1,1,1));
  }
  
  data_t KSB::___a4(const Moment8& Q) const {
    return 4*(Q(0,0,0,0,0,0,1,1)-2*Q(0,0,0,0,1,1,1,1)+Q(0,0,1,1,1,1,1,1));
  }

  data_t KSB::___a5(const Moment8& Q) const {
    return 2*(Q(0,0,0,0,0,0,0,1)-3*Q(0,0,0,0,0,1,1,1)+3.*Q(0,0,0,1,1,1,1,1)-Q(0,1,1,1,1,1,1,1));
  }

  
  // eq. (22)
  complex<data_t> KSB::__p(const KSB& star) const {
    return star.P_sm.invert() * star.chi;
  }
  // eq. (12)
  complex<data_t> KSB::gamma() const {
    return P_sh.invert()*chi;
  }
  // P_sh stricktly to I order
  complex<data_t> KSB::gamma_first() const {
    return P1.invert()*chi;
  }
  //using trace P^g
  complex<data_t> KSB::gammaTr()const{
    return (2./(P_sh(0,0)+P_sh(1,1)))*chi;
  }
  //using trace first order
  complex<data_t> KSB::gammaTr1()const{
    return (2./(P1(0,0)+P1(1,1)))*chi;
  }
  // eqs. (23) and (27, no delta)
  complex<data_t> KSB::gamma(const KSB& psf) const {
    NumMatrix<data_t> P_gamma = P_sh - psf.P_sh*psf.P_sm.invert()*P_sm;
    complex<data_t> p = __p(psf);
    return P_gamma.invert()*(chi - P_sm*p);
  }

  complex<data_t> KSB::gamma_first(const KSB& psf) const {
    NumMatrix<data_t> P_gamma = P1 - psf.P1*psf.P_sm.invert()*P_sm;
    complex<data_t> p = __p(psf);
    return P_gamma.invert()*(chi - P_sm*p);
  }

  complex<data_t> KSB::gammaTr(const KSB& psf) const {
    NumMatrix<data_t> P_gamma = P_sh - psf.P_sh*psf.P_sm.invert()*P_sm;
    complex<data_t> p = __p(psf);
    return (2./(P_gamma(0,0)+P_gamma(1,1)))*(chi - P_sm*p);
  }
  
  complex<data_t> KSB::gammaTr1(const KSB& psf) const {
    NumMatrix<data_t> P_gamma = P1 - psf.P1*psf.P_sm.invert()*P_sm;
    complex<data_t> p = __p(psf);
    return (2./(P_gamma(0,0)+P_gamma(1,1)))*(chi - P_sm*p);
  }


  // *** KSB goes non-linear here ! ***

  // helper struct
  struct KSBTuple {
    const KSB& ksb1;
    const KSB& ksb2;
  };

  
  // Functor for nonlinear P_sh treatment
  class __Psh_nl : public Minimizer::Functor {
  public :
    __Psh_nl(const KSBTuple& kt_, bool conv_,bool flag_) :
      kt(kt_), convolved(conv_), flag(flag_) {}

    data_t operator()(const NumVector<data_t>& gamma) {
      if (convolved)
	{
	  if(flag)
	    { return deltaChi(gamma);}
	  else
	    {return deltaChi(gamma);}
	}
      else
	{
	  if(flag)
	    {return deltaChi_sh(gamma);}
	  else
	    
	    {return deltaChi_exact(gamma);}
	}
    }
    
 
      
  private:
    bool convolved;
    bool flag;
    const KSBTuple kt;
    
    
    complex<data_t> __chi_sh_exact(const KSB& ksb, const complex<data_t>& gamma) {
      
      //double chi2=real(ksb.chi)*real(ksb.chi)+imag(ksb.chi)*imag(ksb.chi);
      
      //double chi1=real(ksb.chi);
      //double chi2=imag(ksb.chi);
      //double g1=real(gamma);
      //double g2=imag(gamma);
      double g2=abs(gamma*gamma);//real(gamma)*real(gamma)+imag(gamma)*imag(gamma);
      complex<data_t> chi_con=conj(ksb.chi);
	    
      double den=1.+g2-real(gamma*chi_con)-real(gamma*chi_con);
      
      return (ksb.chi*g2-ksb.chi*real(gamma*chi_con)-ksb.chi*real(gamma*chi_con)+gamma+gamma-gamma*gamma*chi_con)/den;
      //return ksb.chi*g2-ksb.chi*real(gamma*chi_con)-ksb.chi*real(gamma*chi_con)+gamma+gamma-gamma*g2-gamma*g2+4.*gamma*real(gamma*chi_con)-gamma*gamma*chi_con;  
    }
    
   
    
    complex<data_t> __chi_sh_nl(const KSB& ksb, const complex<data_t>& gamma) {
      
      double g1=real(gamma);
      double g2=imag(gamma);
      double chi1=real(ksb.chi);
      double chi2=imag(ksb.chi);
      double p11=ksb.P_sh(0,0);
      double p12=ksb.P_sh(0,1);
      double p21=p12;
      double p22=ksb.P_sh(1,1);
      double p11_0=ksb.P1(0,0);
      double p12_0=ksb.P1(0,1);
      double p21_0=p12_0;
      double p22_0=ksb.P1(1,1);
      double p11_1=ksb.P2(0,0);
      double p12_1=ksb.P2(0,1);
      double p21_1=p12_1;
      double p22_1=ksb.P2(1,1);
      double B11=ksb.B(0,0)/ksb.trQ;
      double B12=ksb.B(0,1)/ksb.trQ;
      double B21=B12;
      double B22=ksb.B(1,1)/ksb.trQ;
      double r111=ksb.R(0,0,0);
      double r112=ksb.R(0,0,1);
      double r121=ksb.R(0,1,0);
      double r122=ksb.R(0,1,1);
      double r211=ksb.R(1,0,0);
      double r212=ksb.R(1,0,1);
      double r221=ksb.R(1,1,0);
      double r222=ksb.R(1,1,1);
      double s11=2*ksb.K(0,0)/ksb.trQ;
      double s12=2*ksb.K(0,1)/ksb.trQ;
      double s21=s12;
      double s22=2*ksb.K(1,1)/ksb.trQ;
      double r111_0=ksb.R1(0,0,0);
      double r112_0=ksb.R1(0,0,1);
      double r121_0=ksb.R1(0,1,0);
      double r122_0=ksb.R1(0,1,1);  
      double r211_0=ksb.R1(1,0,0);
      double r212_0=ksb.R1(1,0,1);
      double r221_0=ksb.R1(1,1,0);
      double r222_0=ksb.R1(1,1,1);

      double u1111=(4./3.)*ksb.U1(0,0,0);
      double u1112=(4./3.)*ksb.U1(0,0,1);
      double u1121=(4./3.)*ksb.U1(0,1,0);
      double u1122=(4./3.)*ksb.U1(0,1,1);
      double u1211=(4./3.)*ksb.U1(1,0,0);
      double u1212=(4./3.)*ksb.U1(1,0,1);
      double u1221=(4./3.)*ksb.U1(1,1,0);
      double u1222=(4./3.)*ksb.U1(1,1,1);
      double u2111=(4./3.)*ksb.U2(0,0,0);
      double u2112=(4./3.)*ksb.U2(0,0,1);
      double u2121=(4./3.)*ksb.U2(0,1,0);
      double u2122=(4./3.)*ksb.U2(0,1,1);
      double u2211=(4./3.)*ksb.U2(1,0,0);
      double u2212=(4./3.)*ksb.U2(1,0,1);
      double u2221=(4./3.)*ksb.U2(1,1,0);
      double u2222=(4./3.)*ksb.U2(1,1,1);

      
      double gg=real(gamma)*real(gamma)+imag(gamma)*imag(gamma);
      complex<data_t> chi_con=conj(ksb.chi);
      double e1=(2./ksb.trQ)*(real(gamma)*ksb.pi_+2.*imag(gamma)*ksb.nu_);
      double e2=-gg-(2.0/ksb.trQ)*(real(gamma)*real(gamma)*ksb.K(0,0)+2*real(gamma)*imag(gamma)*ksb.K(0,1)+imag(gamma)*imag(gamma)*ksb.K(1,1))+2.0*(real(gamma)*real(ksb.chi)+imag(gamma)*imag(ksb.chi))-(4.0/ksb.trQ)*(real(gamma)*real(gamma)*ksb.B(0,0)+2*real(gamma)*imag(gamma)*ksb.B(0,1)+imag(gamma)*imag(gamma)*ksb.B(1,1));
      
      double a=1;
      double main1 = g1*(p11-2*B11*gg+g1*(r111+g1*(s11+u1111)+g2*u1112)+g2*(r112+g1*(s12+u1121)+g2*u1122)) + g2*(p12-2*B12*gg+g1*(r121+g1*(s21+u1211)+g2*u1212)+g2*(r122+g1*(s22+u1221)+g2*u1222));
      double cor1  = g1*(p11_0*(e1+e2)+p11_1*e1+g1*r111_0*e1+g2*r112_0*e1)+g2*(p12_0*(e1+e2)+p12_1*e1+g1*r121_0*e1+g2*r122_0*e1)+chi1*gg*e1;
      double other1= chi1*gg+(g2*g2-g1*g1)*(chi1-a*2*B11*g1-a*2*B12*g2)-2*g1*g2*(chi2-2*a*B12*g1-2*a*B22*g2);

      double main2 = g1*(p21-2*B21*gg+g1*(r211+g1*u2111+g2*(s11+u2112))+g2*(r212+g1*u2121+g2*(s12+u2122))) + g2*(p22-2*B22*gg+g1*(r221+g1*u2211+g2*(s21+u2212))+g2*(r222+g1*u2221+g2*(s22+u2222)));
      double cor2  = g1*(p21_0*(e1+e2)+p21_1*e1+g1*r211_0*e1+g2*r212_0*e1)+g2*(p22_0*(e1+e2)+p22_1*e1+g1*r221_0*e1+g2*r222_0*e1)+chi2*gg*e1;
      double other2= chi2*gg-(g2*g2-g1*g1)*(chi2-a*2*B21*g1-a*2*B22*g2)-2*g1*g2*(chi1-2*a*B11*g1-2*a*B12*g2);
    
      return complex<data_t>(main1+cor1+other1, main2+cor2+other2);
    }
    
   
    
    complex<data_t> __chi_g(const KSB& ksb, const KSB& star, const complex<data_t>& gamma) {
      complex<data_t> chi_sh_star = star.P1*gamma;//__chi_sh(star,gamma);
      //complex<data_t> chi_sh_star = __chi_sh_nl(star,gamma);
      return ksb.P_sm*star.P_sm.invert()*chi_sh_star;
    }


    ////////////////////////////////////////////////////////////////////////////

    complex<data_t> __chi_sm(const KSB& ksb, const KSB& star) {
      NumMatrix<data_t> P_sm_star_1 = star.P_sm.invert();
      return complex<data_t>(// 1st component
			     ksb.P_sm(0,0)*(P_sm_star_1(0,0)*real(star.chi) + P_sm_star_1(0,1)*imag(star.chi)) +
			     ksb.P_sm(0,1)*(P_sm_star_1(1,0)*real(star.chi) + P_sm_star_1(1,1)*imag(star.chi)),
			     // 2nd component
			     ksb.P_sm(1,0)*(P_sm_star_1(0,0)*real(star.chi) + P_sm_star_1(0,1)*imag(star.chi)) +
			     ksb.P_sm(1,1)*(P_sm_star_1(1,0)*real(star.chi) + P_sm_star_1(1,1)*imag(star.chi)));
    }
    ///////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////
    double deltaChi_sh(const NumVector<data_t>& g) {
      const complex<data_t>* gamma = reinterpret_cast<const complex<data_t>*>(g.c_array());
      complex<data_t> chi_sh = __chi_sh_nl(kt.ksb1,*gamma);
      return abs(kt.ksb1.chi - chi_sh);
    }

     double deltaChi_exact(const NumVector<data_t>& g) {
      const complex<data_t>* gamma = reinterpret_cast<const complex<data_t>*>(g.c_array());
      complex<data_t> chi_sh = __chi_sh_exact(kt.ksb1,*gamma);
      return abs(kt.ksb1.chi - chi_sh);
    }

    double deltaChi(const NumVector<data_t>& g) {
      const complex<data_t>* gamma = reinterpret_cast<const complex<data_t>*>(g.c_array());
      complex<data_t> chi_sh = __chi_sh_nl(kt.ksb1,*gamma);
      complex<data_t> chi_sm = __chi_sm(kt.ksb1,kt.ksb2);
      complex<data_t> chi_g = __chi_g(kt.ksb1,kt.ksb2,*gamma);
      return abs(kt.ksb1.chi - (chi_sh + chi_sm - chi_g));
    }
  };


  complex<data_t> KSB::gamma_nl(data_t accuracy) const {
    const KSBTuple kt = {*this,*this};
    // initial value from linear treatment
    complex<data_t> gamma0 = kt.ksb1.gamma();
    // accuracy relative w.r.t. gamma0
    accuracy*=abs(gamma0);
      
    // initialize the search parameters and their stepsizes
    NumVector<data_t> gamma(2), stepsize(2);
    gamma(0) = real(gamma0);
    gamma(1) = imag(gamma0);
    // quadratic dependence on gamma0 and in opposite direction
    stepsize(0) = -0.5*real(gamma0)*fabs(real(gamma0));
    stepsize(1) = -0.5*imag(gamma0)*fabs(imag(gamma0));
    
    // set up functor
    __Psh_nl p(kt,false,true);
    // call the minimizert
    Minimizer::Simplex(p, gamma, stepsize, accuracy);
    real(gamma0) = gamma(0);
    imag(gamma0) = gamma(1);
    return gamma0;
  }

  complex<data_t> KSB::gamma_exact(data_t accuracy) const {
    const KSBTuple kt = {*this,*this};
    // initial value from linear treatment
    complex<data_t> gamma0 = kt.ksb1.gamma();
    // accuracy relative w.r.t. gamma0
    accuracy*=abs(gamma0);
      
    // initialize the search parameters and their stepsizes
    NumVector<data_t> gamma(2), stepsize(2);
    gamma(0) = real(gamma0);
    gamma(1) = imag(gamma0);
    // quadratic dependence on gamma0 and in opposite direction
    stepsize(0) = -0.5*real(gamma0)*fabs(real(gamma0));
    stepsize(1) = -0.5*imag(gamma0)*fabs(imag(gamma0));
    
    // set up functor
    __Psh_nl p(kt,false,false);
    // call the minimizert
    Minimizer::Simplex(p, gamma, stepsize, accuracy);
    real(gamma0) = gamma(0);
    imag(gamma0) = gamma(1);
    return gamma0;
  }

  
  complex<data_t> KSB::gamma_second(data_t accuracy) const {
    const KSBTuple kt = {*this,*this};
    // initial value from linear treatment
    complex<data_t> gamma0 = kt.ksb1.gamma();
    // accuracy relative w.r.t. gamma0
    accuracy*=abs(gamma0);
      
    // initialize the search parameters and their stepsizes
    NumVector<data_t> gamma(2), stepsize(2);
    gamma(0) = real(gamma0);
    gamma(1) = imag(gamma0);
    // quadratic dependence on gamma0 and in opposite direction
    stepsize(0) = -0.5*real(gamma0)*fabs(real(gamma0));
    stepsize(1) = -0.5*imag(gamma0)*fabs(imag(gamma0));
    
    // set up functor
    __Psh_nl p(kt,false,false);
    // call the minimizert
    Minimizer::Simplex(p, gamma, stepsize, accuracy);
    real(gamma0) = gamma(0);
    imag(gamma0) = gamma(1);
    return gamma0;
  }

  complex<data_t> KSB::gamma_nl(const KSB& psf, data_t accuracy) const {
    const KSBTuple kt = {*this,psf};
    // initial value from linear treatment
    complex<data_t> gamma0 = kt.ksb1.gamma(kt.ksb2);

    // accuracy relative w.r.t. gamma0
    accuracy*=abs(gamma0);
      
    // initialize the search parameters and their stepsizes
    NumVector<data_t> gamma(2), stepsize(2);
    gamma(0) = real(gamma0);
    gamma(1) = imag(gamma0);
    // quadratic dependence on gamma0 and in opposite direction
    stepsize(0) = -0.5*real(gamma0)*fabs(real(gamma0));
    stepsize(1) = -0.5*imag(gamma0)*fabs(imag(gamma0));
      
    // set up functor
    __Psh_nl p(kt,true,true);
    // call the minimizert
    Minimizer::Simplex(p, gamma, stepsize, accuracy);
    real(gamma0) = gamma(0);
    imag(gamma0) = gamma(1);
    return gamma0;
  }


//   data_t devChi(const KSB& ksb, const KSB& star, complex<data_t>& gamma) {
//     complex<data_t> chi_sh = __chi_sh_nl(ksb,gamma);
//     complex<data_t> chi_sm = __chi_sm(ksb,star);
//     complex<data_t> chi_g = __chi_g(ksb,star,gamma);
//     return abs(ksb.chi - (chi_sh + chi_sm - chi_g));
//   }
//     complex<data_t> gamma;
//     for (data_t g1 = -1; g1 <= 1; g1+=0.005) {
//       real(gamma) = g1;
//       for (data_t g2 = -1; g2 <= 1; g2+=0.005) {
// 	imag(gamma) = g2;
// 	if (g2 >= -sqrt(1-g1*g1) && g2 <= sqrt(1-g1*g1))
// 	  std::cout << g1 << "\t" << g2 << "\t" << devChi(*this,psf,gamma) << std::endl;
// 	else
// 	  std::cout << g1 << "\t" << g2 << "\t" << 0 << std::endl;
//       }
//     }

} // end namespace
