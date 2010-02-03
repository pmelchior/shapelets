#include "../../include/utils/Minimizer.h"
#include <gsl/gsl_multimin.h>
#include <numla/NumMatrix.h>

namespace shapelens {

  // helper structure for linmin
  struct linmin_params {
    NumVector<data_t>& p;
    NumVector<data_t>& p_;
    NumVector<data_t>& xi;
    Minimizer::mini_functor& f;
  };
  
  // evaluates functor at p + x*xi
  double f1dim(const gsl_vector* v, void* params) {
    double x = gsl_vector_get(v,0);
    linmin_params* lp = reinterpret_cast<linmin_params*>(params);
  for (int i=0; i < lp->p.size(); i++)
    lp->p_(i) = lp->p(i) + x*lp->xi(i);
  return (lp->f)(lp->p_);
}

// minimize function along given direction
// this is only a externalized simplex minimizer
// set up by Powell
data_t linmin(gsl_multimin_fminimizer * s, gsl_multimin_function * f, gsl_vector * x, gsl_vector * ss, data_t ftol, unsigned int itmax) {
     
  size_t iter = 0;
  int status;
  double size;
     
  /* Starting point */
  gsl_vector_set (x, 0, 0.5);
  gsl_vector_set_all (ss, 0.2);
  
  gsl_multimin_fminimizer_set (s, f, x, ss);

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status) 
      break;
     
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, ftol);
     
  } while (status == GSL_CONTINUE && iter < itmax);

  // set minimal position
  double xmin = gsl_vector_get (s->x, 0);
  linmin_params* lp = reinterpret_cast<linmin_params*>(f->params);
  for (int i=0; i < lp->p.size(); i++) {
    lp->p(i) = lp->p(i) + xmin*lp->xi(i);
    lp->xi(i) *= xmin;
  }
  return s->fval;
}

// the real minimizer
// taken from Numerical Recipes (sect. 10.5)
data_t Minimizer::Powell(mini_functor& func, NumVector<data_t>& p, data_t ftol, unsigned int itmax) {
    
  // initialization
  int n = p.size(), i, iter = 0,ibig,j;
  data_t del,fp,fptt,t,fret, tiny = 1e-15;
  NumVector<data_t> pt = p, ptt(n), xit(n);
  NumMatrix<data_t> xi(n,n);
  for (int i=0; i<n; i++)
    xi(i,i) = 1;
  
  fret=func(p);

  // initialization for linmin
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  s = gsl_multimin_fminimizer_alloc (T, 1);
  gsl_vector *ss, *x;
  x = gsl_vector_alloc (1);
  ss = gsl_vector_alloc (1);
  gsl_multimin_function f;
  f.n = 1;
  f.f = f1dim;
  NumVector<data_t> p_(p.size());
  linmin_params lp = {p,p_,xit,func};
  f.params = &lp;
     
  while (iter < itmax) {
    fp = fret;
    ibig = 0;
    del = 0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) 
	xit(j)=xi(j,i);
      fptt=fret;
      fret = linmin(s,&f,x,ss,ftol,itmax);
      if (fabs(fptt-fret) > del) {
	del = fptt-fret;
	ibig = i;
      }
    }
    if (2.0*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+tiny) {
      break;
    }

    for (j=0; j<n; j++) {
      ptt(j) = 2.0*p(j)-pt(j);
      xit(j) = p(j)-pt(j);
      pt(j) = p(j);
    }
    fptt=func(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(fret)+fptt)*sqrt(fp-fret-del)-del*sqrt(fp-fptt);
      if (t < 0.0) {
	fret = linmin(s,&f,x,ss,ftol,itmax);
	for (j=0; j<n; j++) {
	  xi(j,ibig) = xi(j,n);
	  xi(j,n) = xit(j);
	}
      }
    }
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  return fret;
}


} // end namespace
