#include "../../include/utils/Minimizer.h"
#include <gsl/gsl_multimin.h>
#include <numla/NumMatrix.h>

namespace shapelens {

  // helper structure for linmin
  struct linmin_params {
    NumVector<data_t>& p;
    NumVector<data_t>& p_;
    NumVector<data_t>& xi;
    Minimizer::Functor& f;
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
// this is an externalized simplex minimizer
data_t linmin(gsl_multimin_function& f, linmin_params& lp, data_t ftol, unsigned int itmax) {
  // set up simplex
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  s = gsl_multimin_fminimizer_alloc (T, 1);
  gsl_vector *ss, *x;
  x = gsl_vector_alloc (1);
  ss = gsl_vector_alloc (1);
  gsl_vector_set (x, 0, 0.5);
  gsl_vector_set (ss, 0, 0.2);
  gsl_multimin_fminimizer_set (s, &f, x, ss);

  size_t iter = 0;
   int status;
  do {
    iter++;
     status = gsl_multimin_fminimizer_iterate(s);
    
    if (status) 
      break;
     
     status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s), ftol);
     
  } while (status == GSL_CONTINUE && iter < itmax);

  // set minimal position
  double xmin = gsl_vector_get (s->x, 0);
  for (int i=0; i < lp.p.size(); i++) {
      lp.p(i) +=  xmin*(lp.xi(i));
      lp.xi(i) *= xmin;
  }
  
  data_t fval = s->fval;
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free(x);
  gsl_vector_free(ss);

  return fval;
}

// the real minimizer
// taken from Numerical Recipes (sect. 10.5)
data_t Minimizer::Powell(Functor& func, NumVector<data_t>& p, data_t ftol, unsigned int itmax) {
    
  // initialization
  int n = p.size(), iter = 0,ibig;
  data_t del,fp,fptt,t,fret, tiny = 1e-15;
  NumVector<data_t> pt = p, ptt(n), xit(n);
  NumMatrix<data_t> xi(n,n);
  for (int i=0; i<n; i++)
    xi(i,i) = 1;
  fret=func(p);

  // initialization for linmin
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
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) 
	xit(j)=xi(j,i);
      fptt=fret;
      fret = linmin(f,lp,ftol,itmax);
      if (fabs(fptt-fret) > del) {
	del = fptt-fret;
	ibig = i;
      }
    }
    if (2.0*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+tiny) {
      break;
    }

    for (int j=0; j<n; j++) {
      ptt(j) = 2.0*p(j)-pt(j);
      xit(j) = p(j)-pt(j);
      pt(j) = p(j);
    }
    fptt=func(ptt);
    if (fptt < fp) {
      t=2.0*(fp-2.0*(fret)+fptt)*sqrt(fp-fret-del)-del*sqrt(fp-fptt);
      if (t < 0.0) {
	fret = linmin(f,lp,ftol,itmax);
	for (int j=0; j<n; j++) {
	  xi(j,ibig) = xi(j,n-1);
	  xi(j,n-1) = xit(j);
	}
      }
    }
    iter++;
  }

  return fret;
}



/// ******* simplex stuff here ********

struct simplex_params {
  NumVector<data_t>& p;
  Minimizer::Functor& func;
};

double functor2gsl(const gsl_vector* v, void* params) {
  simplex_params* sp = reinterpret_cast<simplex_params*>(params);
  for (int i =0; i < sp->p.size(); i++)
    sp->p(i) = gsl_vector_get(v,i);
  return sp->func(sp->p);
}

data_t Minimizer::Simplex(Functor& func, NumVector<data_t>& p, const NumVector<data_t>& steps, data_t ftol, unsigned int itmax) {
  int N = p.size();
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  s = gsl_multimin_fminimizer_alloc (T, N);
  gsl_vector *ss, *x;
 
  // set up x to point to data section of p
  x = gsl_vector_alloc (N);
  double * xp = x->data;
  x->data = p.c_array();

  // the stepsizes: 10% of p(i)
  ss = gsl_vector_alloc (N);
  for (int i=0; i < N; i++)
    gsl_vector_set(ss,i,steps(i));

  // minimization function
  gsl_multimin_function f;
  f.n = N;
  f.f = functor2gsl;
  simplex_params sp = {p,func};
  f.params = &sp;

  // the minimization
  gsl_multimin_fminimizer_set (s, &f, x, ss);
  int iter = 0, status;
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    
    if (status) 
      break;
     
    status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s), ftol);
     
  } while (status == GSL_CONTINUE && iter < itmax);

  // set minimal position
  for (int i=0; i < p.size(); i++) {
    p(i) = gsl_vector_get(s->x,i);
  }
  data_t min = s->fval;
  x->data = xp;
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  return min;
}

} // end namespace
