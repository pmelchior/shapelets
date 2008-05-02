#include <ShapeLensConfig.h>
#include <shapelets/OptimalDecomposite2D.h>
#include <shapelets/CoefficientVector.h>
#include <shapelets/ImageTransformation.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <map>
// for chi^2 minimization
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>

using namespace std;

OptimalDecomposite2D::OptimalDecomposite2D(const Object& O, int innmaxLow, int innmaxHigh, data_t inbetaLow, data_t inbetaHigh) : 
Decomposite2D(2,(O.getSize(0) + O.getSize(1))/(2*8),O), obj(O) {
  // set limits for nmax and beta 
  nmaxLow = GSL_MIN_INT(innmaxLow,innmaxHigh);
  nmaxHigh = GSL_MAX_INT(innmaxLow,innmaxHigh);
  betaLow = GSL_MIN_DBL(inbetaLow, inbetaHigh);
  betaHigh = GSL_MAX_DBL(inbetaLow,inbetaHigh);
  
  // estimators for beta
  beta = (obj.getSize(0) + obj.getSize(1))/(2*8);

  // take the minimum of the axis sizes to get a limit for theta_max
  // and garantee orthogonality during minimization
  image_dimension = GSL_MIN(obj.getSize(0),obj.getSize(1));

  npixels = obj.size();
  optimized = nmaxTrouble = 0;
  flags.reset();

  // whether correlation function has to be considered as termination criterium
  if (ShapeLensConfig::NOISEMODEL == "COVARIANCE")
    noise_correlated = 1;
  else
    noise_correlated = 0; 
  
  // set precision in history for history to 4
  history << std::setprecision(4);
}
  
void OptimalDecomposite2D::optimize() {
  int status;
  time_t t0,t1;
  t0 = time(NULL);
  
  // set nmax at start to 2
  // unless nmaxHigh is smaller than 2 or nmaxLow is larger than 2
  optimalNMax = 2;
  if (nmaxHigh < 2) 
    optimalNMax = nmaxHigh;
  else if (nmaxLow > 2)
    optimalNMax = nmaxLow;
  Decomposite2D::setNMax(optimalNMax);

  bestChiSquare = 0;
  // step 1) find optimal beta for the (2/2) shapelet basis
  status = findOptimalBeta(1);
  // increase orders if minimization doesn't converge
 if (status != 0) {
    Decomposite2D::setNMax(Decomposite2D::getNMax()+2);
    status = findOptimalBeta(1);
    if (status != 0) {
      history << "Minimization doesn't converge (well), probably due to image distortions close to the object." <<  endl;
      flags[7] = 1;
    }
  }
 if (!flags.test(7)) {
   // step 2) increase nmax until chi^2 = 1 or is flat
   // and store optimalNMax for comparison in step 4
   findOptimalNMax(2);
   //Decomposite2D::setNMax(optimalNMax);
   // step 3) find optimal beta for new shapelet basis
   status = findOptimalBeta(3);
   // step 4)
   // we could have run into problems during step 2, therefore chi^2
   // might be larger than 1
   if (!noise_correlated && optimalChiSquare >  1 + Decomposite2D::getChiSquareVariance()) {
     history << "#" << endl << "# Decomposition stopped before chi^2 = 1!" << endl;
     flags[0] = 1;
   }
   if (nmaxTrouble)
     flags[6] = 1;
   // step 5) opposite: if chisquare is too good, reduce nmax
   // and increase in steps of 1 to find fit with chisquare close to 1
   if (!noise_correlated) {
     while (optimalChiSquare < 1) {
       int oldoptimalNMax = Decomposite2D::getNMax();
       int newNMAX = (int) floor(0.75*oldoptimalNMax) - 1;
       if (newNMAX < nmaxLow) newNMAX = nmaxLow;
       Decomposite2D::setNMax(newNMAX);
       history << "#" << endl << "# Checking for lower n_max: chi^2 < 1" << endl;
       history << "# Restarting with n_max = " << newNMAX << "."<< endl;
       findOptimalNMax(5);
       Decomposite2D::setNMax(optimalNMax);
       // step 6) if optimal nmax has decreased now look again for beta and xc
       if (Decomposite2D::getNMax() < oldoptimalNMax) {
	 history << "# Found lower n_max." << endl;
	 status = findOptimalBeta(6);
       }
       else if (Decomposite2D::getNMax() == oldoptimalNMax) {
	 break;
       }
       else 
	 optimalNMax = oldoptimalNMax;
     }
   }
   // step 7) when we use the correlation function of the residuals as termination criterium
   // during 2) or 3) we have found chi^2>1 but corr_res < corr.
   // search for nmax and beta such that corr_res >= corr
   else {
     if (comp_corr < 0) {
       history << "#" << endl << "# Lowering n_max for residuals compatible with expectation" << endl;
       int iter = 1;
       while (comp_corr < 0) {
	 int newNMAX = Decomposite2D::getNMax() - 1;
	 Decomposite2D::setNMax(newNMAX);
	 data_t chisquare = Decomposite2D::getChiSquare();
	 data_t variance = Decomposite2D::getChiSquareVariance();
	 checkCorrelationFunctionFromResiduals();
	 history << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
	 history << chisquare << "\t" << variance << "\t" + comp_corr_string << endl;
	 if (comp_corr >= 0)
	   status = findOptimalBeta(7);
       }
       optimalNMax = Decomposite2D::getNMax();
     }
   }
 }

  if (status == 0) optimized = 1;
  t1 = time(NULL);
  history << "#" << endl << "# Computation time for decomposition: " << t1 - t0 << " seconds" << endl;
}


// see Paper III, eq. 39
void OptimalDecomposite2D::findOptimalNMax(unsigned char step) {
  int iter = 1;
  history << "#" << endl << "# Finding optimal decomposition order n_max";
  history << ", beta = " << beta << std::endl;
  history << "# iter.\tn_max\tchi^2\tsigma(chi^2)";
  if (noise_correlated)
    history<<"\txi_res - xi";
  history<< endl;

  data_t chisquare,newChisquare,derivative_chi2,variance;
  chisquare = Decomposite2D::getChiSquare();
  variance = Decomposite2D::getChiSquareVariance();

  history << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
  history << chisquare << "\t" << variance;
  if (noise_correlated) {
    checkCorrelationFunctionFromResiduals();
    history<<"\t\t" + comp_corr_string;
  }
  history << endl;

// during the first steps of the search only use even orders for a fast decomposition
  // only when decomposition is too good ( chisquare < 1 - sigma), go back and
  // use also odd nmax.
  int increment;
  if (step == 5) increment = 1;
  else increment = 2;
  
  if (chisquare <= 1) {
    optimalNMax = Decomposite2D::getNMax();
    history << "# Optimal decomposition order n_max = " << optimalNMax;
    history << " already reached" << endl;
  }
  while (chisquare > 1) {
    // break during step 2 if order gets too high
    // for a relatively fast findOptimalBeta() here
    if (step == 2 && (Decomposite2D::getNMax() == 6 || Decomposite2D::getNMax()%12 == 0)) {
      history << "# Interrupting here for better estimation of beta" << endl;
      if (Decomposite2D::getNMax() == 6)
	findOptimalBeta(2);
      else
	findOptimalBeta(3);
      if (optimalChiSquare <= 1) {
	optimalNMax = Decomposite2D::getNMax();
	history << "# Optimal decomposition order n_max = " << optimalNMax << endl;
	break;
      }
      else {
	history << "#" << endl << "# Continuing search for optimal n_max" << endl;
	history << "# iter.\tn_max\tchi^2\tsigma(chi^2)";
	if (noise_correlated)
	  history<<"\txi_res - xi";
	history<< endl;
      }
    }
    
    // reached the limit of nmax
    if (Decomposite2D::getNMax()+increment > nmaxHigh) {
      history << "# Stopping at nmax limit = " << nmaxHigh << endl;
      flags[2] = 1;
      if (!nmaxTrouble) {
	optimalNMax = Decomposite2D::getNMax();
      } else {
	Decomposite2D::setNMax(optimalNMax);
	history << "# Returning to previous best fit n_max = " << optimalNMax << endl;
      }
      break;
    }

    // increase shapelet order 
    Decomposite2D::setNMax(Decomposite2D::getNMax()+increment);
    newChisquare = Decomposite2D::getChiSquare();
    variance = Decomposite2D::getChiSquareVariance();

    history << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
    history << newChisquare << "\t" << variance;
    if (noise_correlated) {
      checkCorrelationFunctionFromResiduals();
      history<<"\t\t" + comp_corr_string;
    }
    history << endl;

    // depending on result of chi^2:
    // chisquare is smaller than 1: we've reached the goal
    if (newChisquare <= 1) {
      optimalNMax = Decomposite2D::getNMax();
      history << "# Optimal decomposition order n_max = " << optimalNMax << endl;
      break;
    }

    // correlation function of residuals goes below correlation function of noise
    if (noise_correlated && comp_corr < 0) {
      history << "# Stopping here: pixel correlation in residuals becomes less than expected" << endl;
      optimalNMax = Decomposite2D::getNMax();
      flags[1] = 1;
      break;
    }

    // don't do this during the refinement procedure when chi^2 was already low

    if (step != 5) {
      // flattening: chi^2 does improves less than sigma(chi^2)
      if (ShapeLensConfig::ALLOW_FLATTENING && !nmaxTrouble && fabs(newChisquare - chisquare)/increment < variance) {
	  bestChiSquare = chisquare = newChisquare;
	  optimalNMax = Decomposite2D::getNMax();
	  history << "# chi^2 becomes flat. Stopping search at n_max = " << optimalNMax << "." << endl;
	  flags[1] = 1;
	  break;
      }

      // now decomposition get worse, not a good sign
      // save best nmax and chi2
      if (!nmaxTrouble && newChisquare > chisquare) {
	nmaxTrouble = 1;
	bestChiSquare = chisquare;
	optimalNMax = Decomposite2D::getNMax() - increment;
	derivative_chi2 = (newChisquare - chisquare)/increment;
	history << "# chi^2 becomes worse! Saving best fit n_max = " << optimalNMax << " now." << endl;
      }
      // if chi2 becomes better again remember new best values
      if (nmaxTrouble && newChisquare >= 0 && newChisquare < bestChiSquare) {
	nmaxTrouble = 0;
	bestChiSquare = newChisquare;
	optimalNMax = Decomposite2D::getNMax();
	history << "# Better n_max = " <<  optimalNMax << " found now. Continuing search for optimal n_max." << endl;
      }
      // if increase in chi^2 becomes even bigger: break
      if (nmaxTrouble && (newChisquare - chisquare)/increment > derivative_chi2) {
	history << "# chi^2 becomes increasingly worse. Object underconstrained." << endl;
	history << "# Returning to best fit n_max = " << optimalNMax << endl;
	Decomposite2D::setNMax(optimalNMax);
	flags[4] = 1;
	break;
      }
      // if chi^2 gets negative (nCoeffs >= nPixels): go back to last nmax
      if (newChisquare < 0 || newChisquare == INFINITY) {
	if (!nmaxTrouble) {
	  optimalNMax = Decomposite2D::getNMax() - increment;
	  nmaxTrouble = 1;
	}
	Decomposite2D::setNMax(optimalNMax);
	history << "# nPixels <= nCoeffs! Object underconstrained." << endl;
	history << "# Returning to best fit n_max = " << optimalNMax << endl;
	flags[4] = 1;
	break;
      }
      // if theta_min becomes too small -> undersampling
      if (2*optimalBeta/sqrt(Decomposite2D::getNMax()+1.) < 1) {
	if (!nmaxTrouble) {
	  optimalNMax = Decomposite2D::getNMax() - increment;
	  nmaxTrouble = 1;
	}
	Decomposite2D::setNMax(optimalNMax);
	history << "# 2 Theta_min < 1! Image becomes undersampled." << endl;
	history << "# Returning to best fit n_max = " << optimalNMax << endl;
	flags[3] = 1;
	break;
      }
    }
    // we have to continue...
    chisquare = newChisquare;
    iter++;
  }
}

// used as parameters for getChiSquare_Beta().
// this essentially links back to *this
struct parameters {
  Decomposite2D& d;
};

// Since the minimizer is GSL function written in C, it is not able
// to call a member function by reference.
// So this one is global.
double getChiSquare_Beta (double beta, void *p) {
  parameters * decomp = (parameters *)p;
  decomp->d.setBeta(beta);
  return decomp->d.getChiSquare();
}

bool OptimalDecomposite2D::testBetaLowerLimit(data_t& beta) {
  // undersampling is minimum beta to avoid 2*theta_min < 1
  // see Paper IV, eq. (13)
  data_t undersampling = 0.5*sqrt(Decomposite2D::getNMax()+1.);
  if (beta <= GSL_MAX(undersampling,betaLow)) {
    beta = GSL_MAX(undersampling,betaLow);
    return 0;
  } else
    return 1;
}

bool OptimalDecomposite2D::testBetaUpperLimit(data_t& beta) {
  // geometric is maximum beta to avoid theta_max > image_dimension
  data_t geometric = 0.5*image_dimension/sqrt(Decomposite2D::getNMax() +1.0);
  if (beta >= GSL_MIN(geometric, betaHigh)) {
    beta = GSL_MIN(geometric, betaHigh);
    return 0;
  } else
    return 1;
}

// searches for beta which minimizes chi^2
// beta is contrained by image_dimension/2
// return 0 if minimum has been found, -2 if not, and other numbers for error codes
int OptimalDecomposite2D::findOptimalBeta(unsigned char step) {
  history << "#" << endl << "# Finding optimal beta";
  history << ", n_max = " << Decomposite2D::getNMax() << endl;

  // first check if we already did a minimization at this nmax
  map<int,data_t>::iterator iter = bestBeta.find(Decomposite2D::getNMax());
  if (iter != bestBeta.end()) {
    optimalBeta = beta = iter->second;
    Decomposite2D::setBeta(optimalBeta);
    iter = bestChi2.find(Decomposite2D::getNMax());
    optimalChiSquare = iter->second;
    history<< "# Using already found minimum: chi^2 = " << optimalChiSquare << " at beta = ";
    history << beta;
    if (noise_correlated) {
      checkCorrelationFunctionFromResiduals();
      history<<", xi_res - xi = (" + comp_corr_string +")";
    }
    history << endl; 
    return 0;
  }
  
  // not minimized yet
  else {
    size_t iter = 0;
    int status;
    data_t a,b, accuracy;

    // in step 1) correct beta is completely unknown, 
    // therefore define very loose bounds a,b; 
    // a is chosen to be the minimum scale size for nmax = 2 not subject
    // to random sampling (2 theta_min > 1)
    // since this probably not be the last call of the function, 
    // coarse accuracy is enough
    // in step 2) beta is expected to decrease since nmax is 
    // larger than before, but things are still quite uncertain
    // in steps 3-7) accuracy is improving as is the a priori 
    // knowledge of beta
    // in steps 6,7) we try to find lower nmax, which raises beta, 
    // thus b higher
    switch (step) {
    case 1: a = 0.5*sqrt(3.); b = 2*beta; accuracy = 0.02*beta; break;
    case 2: beta *=0.75; a = 0.66*beta; b = 2*beta; accuracy = 0.02*beta; break;
    case 3: a = 0.8*beta; b = 1.1*beta; accuracy = 0.02*beta; break;
    case 6: a = 0.9*beta; b= 1.2*beta; accuracy = 0.02*beta; break;
    case 7: a = 0.9*beta; b= 1.2*beta; accuracy = 0.02*beta; break;
    }
    
    history << "# initial setup: a = " << a << ", beta = " << beta << ", b = " << b << endl;

    // this ensures beta is within reasonable bounds
    // according to undersampling, orthogonality and evtl. specified bounds on beta
    // in addition, we have to ensure that a < beta < b for the minimizer below
    bool testa = testBetaLowerLimit(a), testb = testBetaUpperLimit(b);
    bool limitTrouble = 0;
    data_t mina=0, maxb=image_dimension;
    testBetaLowerLimit(mina); // mina = minimal allowed scale size
    testBetaUpperLimit(maxb); // maxb = maximal allowed scale size

    if (a < beta && beta < b); // as it should be
    else { // testa == 1 or testb == 1
      if (a==b) {
	// if upper and lower bound on beta are identical, we don't have to do anything
	optimalBeta = beta = a;
	Decomposite2D::setBeta(optimalBeta);
	optimalChiSquare = Decomposite2D::getChiSquare();
	bestBeta[Decomposite2D::getNMax()] = optimalBeta;
	bestChi2[Decomposite2D::getNMax()] = optimalChiSquare;
	history << "# Beta = " << beta << " by external constraints: chi^2 = " << optimalChiSquare << endl;
	return 0;
      }
      else if (b < a) { // one of the constraints (on a or b) is severe
	if (testa) { // a can still be lowered
	  if (mina < b)
	    a = 0.5*(mina+b); // a < b 
	  else
	    limitTrouble = 1;
	}
	else if (testb) { // b can be lowered
	  if (maxb > a)
	    b = 0.5*(maxb+a); // a < b
	  else
	    limitTrouble = 1;
	}
	else
	  limitTrouble = 1;
	if (limitTrouble) {
	  history << "# Constraints on beta to severe. Minimization impossible!" << endl;
	  return -2;
	}
      }
      if (beta < a || beta > b)
	beta = 0.5*(a+b); // a < beta < b
    }
    history << "# after checking limits: a = " << a << ", beta = " << beta << ", b = " << b << endl;

    // for the bracketing used in the minimizer we have to also ensure that
    // f(a) > f(beta) and f(b) > f(beta);
    Decomposite2D::setBeta(beta);
    data_t f = Decomposite2D::getChiSquare();
    Decomposite2D::setBeta(a);
    data_t fa = Decomposite2D::getChiSquare();
    Decomposite2D::setBeta(b);
    data_t fb = Decomposite2D::getChiSquare();
    history << "# starting values: f(a) = " << fa << ", f(beta) = " << f << ", f(b) = " << fb << endl;
    
    if (fa < fb) {
      while (fa <= f) {
	// shift everything towards a
	if (fa != f)
	  b = beta; fb = f;
	beta = a; f = fa;
	// check if a gets out of bounds
	a /=1.5;
	if (testBetaLowerLimit(a)) {
	  Decomposite2D::setBeta(a);
	  fa = Decomposite2D::getChiSquare();
	}
	// if yes: shift it just a little further apart, but make
	// its result so worse that will not be considered later
	else {
	  a = mina;
	  a -= accuracy/2;
	  fa = fb; // since fb > fa and f = fa, this terminates the loop
	}
      }
    } else if (fb < fa) {
      while (fb <= f) {
	// or shift everything towards b
	if (fb != f)
	  a = beta; fa = f;
	beta = b; f = fb;
	// check if b gets out of bounds
	b *= 1.5;
	if (testBetaUpperLimit(b)) {
	  Decomposite2D::setBeta(b);
	  fb = Decomposite2D::getChiSquare();
	}
	// as above: shift it slightly and make it worse enough
	else {
	  b = maxb;
	  b += accuracy/2;
	  fb = fa;
	}
      }
    } else {
      // fa == fb here: this object has a flat chi2 function
      std::cerr << "# chi^2 is flat; detection of minimum failed!" << endl;
      std::terminate();
    }

    history << "# after bracketing: f(a) = " << fa << ", f(beta) = " << f << ", f(b) = " << fb << endl;
    history << "# Setting limits: " << a << " <= beta <= " << b << endl;
    history << "# Starting with beta = " << beta << endl;
    history << "# iter.\tchi^2\tbeta\tdelta(beta)" << endl;

    size_t max_iter = 100;
    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);

    // define the function which should be minimized and its parameters
    gsl_function F;
    F.function = &getChiSquare_Beta;
    parameters params = { *this };
    F.params = &params;

    // initialize the minizer with already found values from above
    gsl_min_fminimizer_set_with_values (s, &F, beta, f, a, fa, b, fb);

    do {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
     
      beta = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);
     
      // the accuracy comes in here:
      // when size is smaller than accuracy we have convergence
      status = gsl_min_test_interval (a, b, accuracy, 0.0);
     
      history << "# " << iter << "\t" << gsl_min_fminimizer_f_minimum(s) << "\t" << beta << "\t" << b-a << endl;

      if (status == GSL_SUCCESS) {
	optimalBeta = beta;
	Decomposite2D::setBeta(beta);
	optimalChiSquare = Decomposite2D::getChiSquare();
	history<< "# Converged to minimum: chi^2 = " << optimalChiSquare << " at beta = ";
	history << optimalBeta;
	if (noise_correlated) {
	  checkCorrelationFunctionFromResiduals();
	  history<<", xi_res - xi = (" + comp_corr_string +")";
	}
	history << endl; 
      }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    // if minimizer did not converge within max_iter:
    // save best solution found by now
    if (status == GSL_CONTINUE) {
      bestBeta[Decomposite2D::getNMax()] = beta;
      bestChi2[Decomposite2D::getNMax()] = gsl_min_fminimizer_f_minimum(s);
    }
    else {
      bestBeta[Decomposite2D::getNMax()] = optimalBeta;
      bestChi2[Decomposite2D::getNMax()] = optimalChiSquare;
      // update weight map for upcoming decompositions
      Decomposite2D::updateWeightMap();
    }

    gsl_min_fminimizer_free (s);
    return status;
  }
}

const NumVector<data_t>& OptimalDecomposite2D::getResiduals() {
  if (!optimized) optimize();
  return Decomposite2D::getResiduals();
}

const NumVector<data_t>& OptimalDecomposite2D::getModel() {
  if (!optimized) optimize();
  return Decomposite2D::getModel();
}

// computes correlation function from residuals and compares it point by point
// with the one stored in obj
void OptimalDecomposite2D::checkCorrelationFunctionFromResiduals() {
  const CorrelationFunction& xi = obj.getCorrelationFunction();
  // look for the bandwidth of V to decide on the size of the box within
  // which to compute the correlation function
  unsigned int bandwidth = obj.getPixelCovarianceMatrix().getBandwidth();
  unsigned int size;
  switch(bandwidth) {
  case 1: size=0; break;
  case 5: size=1; break;
  case 9: size=1; break;
  default: size=2; break;
  }
  
  CorrelationFunction xi_res(Decomposite2D::getResiduals(),obj.getGrid(),size);
  NumVector<data_t> corr = xi.getCorrelationFunction(), sigma = xi.getCorrelationError(),
    corr_res = xi_res.getCorrelationFunction(), sigma_res = xi_res.getCorrelationError();

  comp_corr = 0;
  comp_corr_string = "";
  for (uint i=0; i<corr_res.size(); i++) {
    if (corr(i) + sigma(i) < corr_res(i) - sigma_res(i)) {
      comp_corr++;
      comp_corr_string += "+";
    }
    else if (corr(i) - sigma(i) > corr_res(i) + sigma_res(i)) {
      comp_corr--;
      comp_corr_string += "-";
    }
    else 
      comp_corr_string += "=";
  }
}
 
const CoefficientVector<data_t>& OptimalDecomposite2D::getCoeffs() {
  if (!optimized) optimize();
  return Decomposite2D::getCoeffs();
}

NumMatrix<data_t> OptimalDecomposite2D::getCovarianceMatrix() {
  if (!optimized) optimize();
  return Decomposite2D::getCovarianceMatrix();
}

data_t OptimalDecomposite2D::getOptimalBeta() {
  if (!optimized) optimize();
  return optimalBeta;
}
data_t OptimalDecomposite2D::getOptimalChiSquare() {
  if (!optimized) optimize();
  return optimalChiSquare;
}

const bitset<8>& OptimalDecomposite2D::getDecompositionFlags() {
  return flags;
}


// ### REGULARIZATION ###
// this part implements the regularization procedure which minimizes the
// negative flux in the reconstruction

struct reg_parameters {
  Decomposite2D& d;
  History& history;
  const NumVector<data_t>& reco;
  const NumVector<data_t>& residual;
  data_t& beta;
  int& nmax;
  data_t& f;
  data_t& chi2;
  data_t& lambda;
  data_t& H;
  data_t& R;
  data_t& wantedR;
  bool updateCoeffs;
  data_t betaLow;
  data_t betaHigh;
};

void computeFluxes(const NumVector<data_t>& model, data_t& negFlux, data_t& posFlux, data_t& R, data_t& H) {
  negFlux = 0;
  posFlux = 0;
  for(int i=0; i< model.size(); i++) {
    if (model(i) < 0) {
      negFlux += model(i);
    }
    if (model(i) > 0) {
      posFlux += model(i);
    }
  }
  if (negFlux<0)
    R = -negFlux/posFlux;
  else 
    R = 0;
  // H is the penalty function;
  // it should be non-linear and react fast on small values of R
  // to speed up the optimization
  H = acosh(1+R);
}

double f_beta (const gsl_vector *v, void *params) {
  reg_parameters * p = (reg_parameters *)params;

  data_t beta = gsl_vector_get(v, 0);
  const NumVector<data_t>& reco = p->reco;
  const NumVector<data_t>& residual = p->residual;

  // constraints on beta: betaLow <= beta <= betaHigh
  if (beta >= p->betaLow && beta <= p->betaHigh) {
    p->d.setBeta(beta);
    p->d.setNMax(p->nmax);
    // at this position: update coeffs for the chi^2 problem
    // assumption: they are very similar to best fit for f
    // as long as lamda is small
    p->d.updateCoeffs(p->updateCoeffs);

    // now compute chi^2, model values and residuals
    p->chi2 = p->d.getChiSquare();
  
    data_t posFlux, negFlux;
    computeFluxes(reco,negFlux,posFlux, p->R, p->H);
    p->f = (p->chi2 + (p->lambda)*p->H);
    p->beta = beta;
  } else 
    p->f = INFINITY ; // make f worse if beta is out of bounds
  
  return (p->f);
}

double f_coeffs(const gsl_vector *v, void *params) {
  reg_parameters * p = (reg_parameters *)params;
  NumVector<data_t>& coeffs = p->d.accessCoeffs();

  // copy coeffs from vector into Decomposite2D
  for (int i=0; i< v->size; i++)
    coeffs(i) = gsl_vector_get(v, i);
  
  // force update of model and residuals
  p->d.updateModelResiduals();
  // now compute chi^2
  data_t chi2 = p->d.getChiSquare();
  
  data_t posFlux, negFlux, R, f, H;
  computeFluxes(p->reco,negFlux,posFlux,R, H);
  f = chi2 + (p->lambda)*H;
  // since we want to check for the values at best f
  // we store these only if new f is better than actual f
  if (f < p->f) {
    p->f = f;
    p->chi2 = chi2;
    p->R = R;
    p->H = H;
  }
  return f;
}

void minimize_beta(reg_parameters& p) {
  p.history << "#" << endl << "# Minimize w.r.t beta" << endl;
  gsl_vector *x, *ss;
  gsl_multimin_function my_func;
  my_func.f = &f_beta;
  my_func.n = 1;
  my_func.params = &p;
  
  data_t beta = p.beta;
  x = gsl_vector_alloc (1);
  ss = gsl_vector_alloc (1);
  gsl_vector_set (x, 0, beta);
  gsl_vector_set (ss, 0, 0.05*beta);

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, x->size);
  gsl_multimin_fminimizer_set (s, &my_func, x, ss);

  int iter =0, status;
  do {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      if (status) {
	std::cout << gsl_strerror (status) << std::endl;
	break;
      }
      // required accuracy: 1% beta
      status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s),0.01*beta);
      
      if (status == GSL_SUCCESS) {
	p.history << "# new minimum found:" << endl;
      }
  } while (status == GSL_CONTINUE); 
  // update all quantities with new beta
  // this is neccessary, since p contains values from last call, not from best f
  f_beta(s->x,&p);
  p.history << "# " << iter << ":\t beta = " << gsl_vector_get (s->x, 0)<< "\t f = " << s->fval;
  p.history << "\t chi^2 = " << p.chi2 <<  "\t R = "<< p.R << std::endl;

  // clean up
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);
}

void minimize_coeffs(reg_parameters& p) {
  p.history << "#" << endl << "# Minimizing w.r.t. shapelet coefficients" << endl;
  const NumVector<data_t>& coeffVector = p.d.getCoeffs();
  int nCoeffs = coeffVector.size();
  gsl_vector *x, *ss;

  gsl_multimin_function my_func;
  my_func.f = &f_coeffs;
  my_func.n = nCoeffs;
  my_func.params = &p;
  
  x = gsl_vector_alloc (nCoeffs);
  ss = gsl_vector_alloc (nCoeffs);
  for (int i=0; i<nCoeffs; i++) {
    gsl_vector_set (x, i, coeffVector(i));
    gsl_vector_set (ss, i, 0.05*coeffVector(i));
  }
  
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, nCoeffs);
  gsl_multimin_fminimizer_set (s, &my_func, x, ss);

  int iter =0, status;
  do {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);

      if (status) {
	std::cout << gsl_strerror (status) << std::endl;
	break;
      }
      
      if (p.R <= p.wantedR && p.chi2 <= 1)
	status = GSL_SUCCESS;
      else
 	status = GSL_CONTINUE;


      if (status == GSL_SUCCESS) {
      	p.history << "# new minimum found:" << endl;
      }
  } while (status == GSL_CONTINUE && iter < 1000);
 
  // update coeffs and all other quantities in p
  f_coeffs(s->x,&p);

  p.history << "# " << iter << ":\t f = " << s->fval << "\t chi^2 = " << p.chi2 << "\t R = " << p.R << std::endl;
  
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);
}

// append actual regularization results to vector
void OptimalDecomposite2D::appendRegResults(std::vector<regResults>& results, int nmax, data_t f, data_t lambda, data_t beta, data_t chi2, data_t R, const NumVector<data_t>& coeffs) {
  regResults reg = { nmax, f, lambda, beta, chi2, R , coeffs};
  results.push_back(reg);
}

// find nmax where f is minimal in results
int OptimalDecomposite2D::findNMaxofBestF(std::vector<regResults>& results) {
  int maxnum = results.size()-1;
  data_t bestf = results[maxnum].f;
  int bestnmax = results[maxnum].nmax;
  for (int i=maxnum-1; i>=0; i--) {
    if (results[i].f < bestf) {
      bestf = results[i].f;
      bestnmax = results[i].nmax;
    }
  }
  return bestnmax;
}

// find result set where R is minimal and chi2<=1
int OptimalDecomposite2D::findSuboptimalResultIndex(std::vector<regResults>& results) {
  int maxnum = results.size()-1;
  data_t bestR = results[0].R;
  int bestIndex = 0;
  for (int i=1; i<=maxnum; i++) {
    if (results[i].R < bestR && results[i].chi2 <= 1) {
      bestR = results[i].R;
      bestIndex = i;
    }
  }
  return bestIndex;
}
  
data_t OptimalDecomposite2D::regularize(data_t wantedR) {
  if (!optimized) optimize();

  int nmax;
  data_t negFlux, posFlux, R, initialR, H, initialH, chi2, initialChi2, beta;
  const NumVector<data_t>& model = Decomposite2D::getModel();
  const NumVector<data_t>& residuals = Decomposite2D::getResiduals();
  const NumVector<data_t>& coeffs = Decomposite2D::getCoeffs();

  computeFluxes(model, negFlux, posFlux, R, H);
  initialR = R;
  initialH = H;
  chi2 = initialChi2 = optimalChiSquare;
  beta = optimalBeta;
  nmax = optimalNMax;
  history << "#" << std::endl << "# Regularization: minimizing negative Flux " << std::endl;
  history << "# Requirement: R = negative_Flux/positive_Flux < " << wantedR << std::endl;
  history << "# Minimizing f = chi^2 + lambda*H, with H = acosh(1+R)" << std::endl;
  history << "# initial R = " << initialR << ", => initial H = " << H << std::endl;

  // nothing to do anymore?
  if (wantedR >= R ) {
    history << "# Regularization not neccessary: R already lower than wanted R = " << wantedR << std::endl;
  }
  else if (!noise_correlated && optimalChiSquare > 1) {
    history << "# Regularization not appropriate: chi^2 > 1"<< std::endl;
  }
  else {
    // lambda is adjusted such that chi2 and R are compatible;
    // for good performance lower lambda by factor 3
    // because model is less strong restricted by data/chi2
    data_t lambda = 0.3*chi2/H;
    data_t initialLambda = lambda;
    data_t f = optimalChiSquare + lambda*H;
    history << "# Starting with lambda = " << lambda << std::endl;

    // set up vector for storing the numerous regularization results
    std::vector<regResults> results;
    // store values before regularization now: lambda = 0
    appendRegResults(results,nmax,f,0,beta,chi2,R,coeffs);

    reg_parameters par = { *this, history, model, residuals, beta, nmax, f, chi2, lambda, H, R, wantedR, 1, betaLow, betaHigh};

    time_t t0,t1;
    t0 = time(NULL);
    minimize_beta(par);
    minimize_coeffs(par);
    
    if (noise_correlated) 
      checkCorrelationFunctionFromResiduals();
    appendRegResults(results,nmax,f,lambda,beta,chi2,R,coeffs);
    t1 = time(NULL);
    bool trouble = 0;
    // check if further iterations are neccessary
    while ((!noise_correlated &&(R > wantedR || chi2 > 1)) || (noise_correlated && comp_corr>=0)) {
      // introduce arbitrary time constraint
      if (t1-t0 > 900) {
	history << "# Regularization takes too long, stopping here!" << endl;
	trouble = 1;
	break;
      }	
      // model is to strongly constrained by data
      // or just not capable of fulfilling both conditions:
      // increase nmax to get more flexibility
      if ((!noise_correlated && chi2>1 && R > wantedR) || (noise_correlated && R > wantedR && comp_corr >=0)) {
	// find out if increasing nmax would actually help:
	// only if best fit f was obtained at this nmax,
	// it's likely that increasing nmax helps
	int bestnmax = findNMaxofBestF(results);
	if (bestnmax == nmax) {
	  if (nmax + 1 <= nmaxHigh) {
	    nmax+=1;
	    history << "# Object too strongly constrained by data, setting nmax = " << nmax << " and lambda = " << lambda << std::endl;
	    par.updateCoeffs = 1;
	  } else {
	    history << "# Object too strongly constrained by data, but nmax limit reached." << endl;
            history << "# Stopping here and reverting to previous best fit values:" << endl;
	    trouble = 1;
	    break;
	  }
	} else {
	  history << "# f increased during last increase of nmax: regularization constraints to strict for this object." << endl;
          history << "# Reverting to previous best values:" << endl;
	  trouble = 1;
	  break;
	}
      }
      else {
	// no good fit to data:
	// regularization is too extreme, lower lambda
	// also allow use of best fit (w.r.t chi2) coeffs when searching for beta
	// otherwise chi2 cannot be lowered sufficiently
	if (chi2 > 1) {
	  lambda /= GSL_MAX(2,1 + fabs(initialChi2-chi2)*100);
	  history << "# Regularization too strong, setting lambda = " << lambda << std::endl;
	  par.updateCoeffs = 1;
	}
	// reg. too weak, problems is given by high R:
	// fix coeffs when searching best beta. It uses the coeffs from the last
	// search for optimal coeffs (w.r.t. f)
	else if (R > wantedR && chi2 <= 1) {
	  lambda *= GSL_MAX(1.5,initialH/fabs(initialH-H));
	  history << "# Regularization too weak, setting lambda = " << lambda << std::endl;
	  par.updateCoeffs = 0;
	}
      }
      minimize_beta(par);
      minimize_coeffs(par);
      if (noise_correlated) 
	checkCorrelationFunctionFromResiduals();
      appendRegResults(results,nmax,f,lambda,beta,chi2,R,coeffs);
      t1 = time(NULL);
    }
    if (!trouble) {
	optimalNMax = nmax;
	optimalBeta = beta;
	optimalChiSquare = chi2;
    	history << "# Regularization converged: chi^2 = " << optimalChiSquare << ", R = "<< R << std::endl;
    }
    // we have to revert to those reg. parameters nmax and lambda, where 
    // chi2 <=1 and R minimal
    else {
      int bestIndex = findSuboptimalResultIndex(results);
      regResults bestResult = results[bestIndex];
      // reset to those values
      lambda = bestResult.lambda;
      f = bestResult.f;
      R = bestResult.R;
      optimalChiSquare = chi2 = bestResult.chi2;
      optimalNMax = nmax = bestResult.nmax;
      optimalBeta = beta = bestResult.beta;
      Decomposite2D::setNMax(nmax);
      Decomposite2D::setBeta(beta);
      Decomposite2D::accessCoeffs() = bestResult.coeffs;
      Decomposite2D::updateCoeffs(0);
      history << "# nmax = " << optimalNMax << ", lambda = " << lambda << " -> beta = " << beta << ", f = " << f << ", chi^2 = " << chi2 << ", R = " << R << std::endl;
    }
    history << "#" << endl << "# Computation time for regularization: " << t1 - t0 << " seconds" << std::endl;
  }
  return R;
}

std::string OptimalDecomposite2D::getHistory() {
  return history.str();
}
