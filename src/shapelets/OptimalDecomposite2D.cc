#include "../../include/ShapeLensConfig.h"
#include "../../include/shapelets/OptimalDecomposite2D.h"
#include "../../include/shapelets/CoefficientVector.h"
#include "../../include/shapelets/ImageTransformation.h"
#include <math.h>
#include <time.h>
#include <iomanip>
#include <map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>

using namespace shapelens;
using namespace std;

OptimalDecomposite2D::OptimalDecomposite2D(const Object& obj, Composite2D& c) : 
  Decomposite2D(obj,c)
{
  // take the minimum of the axis sizes to get a limit for theta_max
  // and garantee orthogonality during minimization
  image_dimension = GSL_MIN(obj.getSize(0),obj.getSize(1));
  npixels = obj.size();

  // whether correlation function has to be considered as termination criterium
  if (ShapeLensConfig::NOISEMODEL == "COVARIANCE")
    noise_correlated = 1;
  else
    noise_correlated = 0; 
  
  // set precision in history for history to 4
  history << std::setprecision(4);
  
  // start optimization
  nmaxTrouble = 0;
  flags.reset();
  optimize();
}
  
void OptimalDecomposite2D::optimize() {
  int status;
  time_t t0,t1;
  t0 = time(NULL);
  
  // set nmax at start to 2
  // unless ShapeLensConfig::NMAX_HIGH is smaller than 2 or ShapeLensConfig::NMAX_LOW is larger than 2
  optimalNMax = 2;
  if (ShapeLensConfig::NMAX_HIGH < 2) 
    optimalNMax = ShapeLensConfig::NMAX_HIGH;
  else if (ShapeLensConfig::NMAX_LOW > 2)
    optimalNMax = ShapeLensConfig::NMAX_LOW;
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
       if (newNMAX < 0) break;
       if (newNMAX < ShapeLensConfig::NMAX_LOW) newNMAX = ShapeLensConfig::NMAX_LOW;
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
     // step 5a) if chi^2 < 0 at nmax = 0: try what happens with chi2
     // if we set the only shapelet coefficient to zero: no model = pure data
     if (optimalNMax == 0 && optimalChiSquare < 1) {
       history << "#" << endl << "# Testing model evidence: ";
       CoefficientVector<data_t> coeffs = Decomposite2D::C2D.getCoeffs();
       data_t c0 = coeffs(0);
       coeffs(0) = 0;
       Decomposite2D::setCoeffs(coeffs);
       Decomposite2D::fixCoeffs(true);
       data_t newChi2 = Decomposite2D::getChiSquare();
       if (newChi2 < 1) {
	 history << "model insignificant and set to zero" << endl;
	 flags[7] = 1;
       }
       else { // go back to old coeffs again
	 history << "model significant" << endl;
	 coeffs(0) = c0;
	 Decomposite2D::setCoeffs(coeffs);
       }
       Decomposite2D::fixCoeffs(false);
     }
   }
   // step 6) when we use the correlation function of the residuals as termination criterium
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
 t1 = time(NULL);
 history << "#" << endl << "# Computation time for decomposition: " << t1 - t0 << " seconds" << endl;
}


// see Paper III, eq. 39
void OptimalDecomposite2D::findOptimalNMax(unsigned char step) {
  int iter = 1;
  history << "#" << endl << "# Finding optimal decomposition order n_max";
  history << ", beta = " << Decomposite2D::getBeta() << std::endl;
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
    if (Decomposite2D::getNMax()+increment > ShapeLensConfig::NMAX_HIGH) {
      history << "# Stopping at nmax limit = " << ShapeLensConfig::NMAX_HIGH << endl;
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
    if (newChisquare <= 1 && newChisquare > 0) {
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
      // decomposition get worse, not a good sign
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
      if (2*Decomposite2D::getBeta()/sqrt(Decomposite2D::getNMax()+1.) < 1) {
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
      // flattening: chi^2 does improves less than sigma(chi^2)
      if (ShapeLensConfig::ALLOW_FLATTENING && !nmaxTrouble && fabs(newChisquare - chisquare)/increment < variance) {
	  bestChiSquare = chisquare = newChisquare;
	  optimalNMax = Decomposite2D::getNMax();
	  history << "# chi^2 becomes flat. Stopping search at n_max = " << optimalNMax << "." << endl;
	  flags[1] = 1;
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
  double betaLow;
  double betaHigh;
};

// Since the minimizer is GSL function written in C, it is not able
// to call a member function by reference.
// So this one is global.
double getChiSquare_Beta (const gsl_vector *v, void *p) {
  double beta, result;
  parameters * decomp = (parameters *)p;
  beta = gsl_vector_get(v, 0);
  // constraints on beta: betaLow <= beta <= betaHigh
  if (beta >= decomp->betaLow && beta <= decomp->betaHigh) {
    decomp->d.setBeta(beta);
    result =  decomp->d.getChiSquare();
  } else
    result = INFINITY ; // make chi^2 worse if beta is out of bounds
   return result;
}

data_t OptimalDecomposite2D::getBetaLimit(bool upper) {
  if (!upper) {
    // undersampling is minimum beta to avoid 2*theta_min < 1
    // see Paper IV, eq. (13)
    data_t undersampling = 0.5*sqrt(Decomposite2D::getNMax()+1.);
    return GSL_MAX(undersampling,ShapeLensConfig::BETA_LOW);
  } else {
    // geometric is maximum beta to avoid theta_max > image_dimension
    data_t geometric = 0.5*image_dimension/sqrt(Decomposite2D::getNMax() +1.0);
    return GSL_MIN(geometric,ShapeLensConfig::BETA_HIGH);
  }
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
    optimalBeta = iter->second;
    Decomposite2D::setBeta(optimalBeta);
    Decomposite2D::computeModel(); // we need to update coeffs, model, covariance ...
    iter = bestChi2.find(Decomposite2D::getNMax());
    optimalChiSquare = iter->second;
    history<< "# Using already found minimum: chi^2 = " << optimalChiSquare << " at beta = ";
    history << optimalBeta;
    if (noise_correlated) {
      checkCorrelationFunctionFromResiduals();
      history<<", xi_res - xi = (" + comp_corr_string +")";
    }
    history << endl; 
    return 0;
  }
  else if (ShapeLensConfig::BETA_LOW == ShapeLensConfig::BETA_HIGH) {
    // if upper and lower bound on beta are identical, we don't have to do anything
    optimalBeta = ShapeLensConfig::BETA_LOW;
    Decomposite2D::setBeta(optimalBeta);
    optimalChiSquare = Decomposite2D::getChiSquare();
    bestBeta[Decomposite2D::getNMax()] = optimalBeta;
    bestChi2[Decomposite2D::getNMax()] = optimalChiSquare;
    history << "# Beta = " << optimalBeta << " by external constraints: chi^2 = " << optimalChiSquare << endl;
    flags[5] = 1;
    return 0;
  }
  // not minimized yet
  else {
    // in step 1) correct beta is completely unknown, 
    // therefore we set it the the average of the upper and lower limit on beta
    // in step 2) beta is expected to decrease since nmax is 
    // larger than before, but things are still quite uncertain
    // in steps 3-7) accuracy is improving as is the a priori 
    // knowledge of beta
    // in steps 6,7) we try to find lower nmax
    data_t beta = Decomposite2D::getBeta(), betaMin = getBetaLimit(0), betaMax = getBetaLimit(1);
    data_t stepsize;

    switch (step) {
    case 1: beta = 0.5*(betaMax + betaMin); stepsize = 0.21*(betaMax + betaMin); break;
    case 2: beta *= 0.75; stepsize = 0.2 * beta; break;
    case 3: stepsize = 0.1 * beta; break;
    case 6: stepsize = 0.05 * beta; break;
    case 7: stepsize = 0.05 * beta; break;
    }
    data_t accuracy = ShapeLensConfig::DELTA_BETA * beta;

    history << "# Starting with beta = " << beta << endl;
    history << "# iter.\tchi^2\tbeta\tdelta(beta)" << endl;

    // initialize minimizer
    size_t iter = 0, max_iter = 25, np = 1;
    int status;
    data_t size;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x, *min;
    parameters params = { *this, betaMin, betaMax };
    gsl_multimin_function F;

    // define initial vertex vector
    ss = gsl_vector_alloc (np);
    // define primary stepsize
    gsl_vector_set (ss,0,stepsize);

    // Starting point
    x = gsl_vector_alloc (np);
    gsl_vector_set (x, 0, beta);

    // define the function which should be minimized and its parameters
    F.f = &getChiSquare_Beta;
    F.n = np;
    F.params = &params;
    s = gsl_multimin_fminimizer_alloc (T, np);
    gsl_multimin_fminimizer_set (s, &F, x, ss);

    do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
     
      if (status)
	break;

      // the accuracy comes in here:
      // when size is smaller than accuracy we have convergence
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size,accuracy);

      history << "# " << iter << "\t" << s->fval << "\t" << gsl_vector_get (s->x,0) << "\t" << size << endl;
      if (status == GSL_SUCCESS) {
	min = gsl_multimin_fminimizer_x(s);
	optimalBeta = gsl_vector_get (min,0);
	optimalChiSquare = gsl_multimin_fminimizer_minimum(s);
	Decomposite2D::setBeta(optimalBeta);
	Decomposite2D::computeModel();
	history << "# Converged to minimum: chi^2 = " << optimalChiSquare << " at beta = " << optimalBeta;
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
      min = gsl_multimin_fminimizer_x(s);
      bestBeta[Decomposite2D::getNMax()] = gsl_vector_get (min,0);
      bestChi2[Decomposite2D::getNMax()] = gsl_multimin_fminimizer_minimum(s);
    }
    else {
      bestBeta[Decomposite2D::getNMax()] = optimalBeta;
      bestChi2[Decomposite2D::getNMax()] = optimalChiSquare;
      // update weight map for upcoming decompositions
      Decomposite2D::updateWeightMap();
    }

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return status;
  }
}

// computes correlation function from residuals and compares it point by point
// with the one stored in obj
void OptimalDecomposite2D::checkCorrelationFunctionFromResiduals() {
  const CorrelationFunction& xi = (Decomposite2D::obj).xi;
  
  CorrelationFunction xi_res(Decomposite2D::getResiduals(),0,xi.getMaxLength());
  const std::map<Point2D<int>, data_t>& corr = xi.getCorrelationFunction(), sigma = xi.getCorrelationError(), corr_res = xi_res.getCorrelationFunction(), sigma_res = xi_res.getCorrelationError();

  comp_corr = 0;
  comp_corr_string = "";
  for (std::map<Point2D<int>, data_t>::const_iterator iter = corr.begin(); iter != corr.end(); iter++) {
    if (iter->second + sigma.find(iter->first)->second < corr_res.find(iter->first)->second - sigma_res.find(iter->first)->second) {
      comp_corr++;
      comp_corr_string += "+";
    }
    else if (iter->second - sigma.find(iter->first)->second > corr_res.find(iter->first)->second + sigma_res.find(iter->first)->second) {
      comp_corr--;
      comp_corr_string += "-";
    }
    else 
      comp_corr_string += "=";
  }
}

data_t OptimalDecomposite2D::getOptimalChiSquare() {
  return optimalChiSquare;
}

const bitset<8>& OptimalDecomposite2D::getDecompositionFlags() {
  return flags;
}

const History& OptimalDecomposite2D::getHistory() {
  return history;
}


