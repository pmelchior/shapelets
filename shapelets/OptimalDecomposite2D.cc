#include <OptimalDecomposite2D.h>
#include <MatrixManipulations.h>
#include <ImageTransformation.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <map>
// for chi^2 minimization
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>

using namespace std;

OptimalDecomposite2D::OptimalDecomposite2D(const Object& O, int innmaxLow, int innmaxHigh, double inbetaLow, double inbetaHigh) : 
Decomposite2D(2,(O.getSize(0) + O.getSize(1))/(2*8),O), obj(O) {
  // set limits for nmax and beta 
  nmaxLow = GSL_MIN_INT(innmaxLow,innmaxHigh);
  nmaxHigh = GSL_MAX_INT(innmaxLow,innmaxHigh);
  betaLow = GSL_MIN_DBL(inbetaLow, inbetaHigh);
  betaHigh = GSL_MAX_DBL(inbetaLow,inbetaHigh);

  // estimators for beta from FitsImage
  beta = (obj.getSize(0) + obj.getSize(1))/(2*8);

  // take the minimum of the axis sizes to get a limit for theta_max
  // and garantee orthogonality during minimization
  image_dimension = GSL_MIN(obj.getSize(0),obj.getSize(1));

  npixels = obj.size();
  optimized = nmaxTrouble = 0;
  flag = 0;

  // whether correlation function has to be considered as termination criterium
  if (obj.getNoiseModel().compare("COVARIANCE")==0)
    noise_correlated = 1;
  else
    noise_correlated = 0; 
  
  // set precision in text for history to 4
  text << std::setprecision(4);
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
      history.append("Minimization doesn't converge (well), probably due to image distortions close to the object.\n");
      flag = 5;
    }
  }
 if (flag != 5) {
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
     history.append("#\n# Decomposition stopped before chi^2 = 1!\n");
     if (flag==0) flag = 4;
   }
   // step 5) opposite: if chisquare is too good, reduce nmax
   // and increase in steps of 1 to find fit with chisquare close to 1
   if (!noise_correlated) {
     while (optimalChiSquare < 1) {
       int oldoptimalNMax = Decomposite2D::getNMax();
       int newNMAX = (int) floor(0.75*oldoptimalNMax)-1;
       if (newNMAX < nmaxLow) newNMAX = nmaxLow;
       Decomposite2D::setNMax(newNMAX);
       text << "#" << endl << "# Checking for lower n_max: chi^2 < 1" << endl;
       text << "# Restarting with n_max = " << newNMAX << "."<< endl;
       history.append(text);
       findOptimalNMax(5);
       Decomposite2D::setNMax(optimalNMax);
       // step 6) if optimal nmax has decreased now look again for beta and xc
       if (Decomposite2D::getNMax() < oldoptimalNMax) {
	 history.append("# Found lower n_max.\n");
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
       history.append("#\n# Lowering n_max for residuals compatible with expectation\n");
       int iter = 1;
       while (comp_corr < 0) {
	 int newNMAX = Decomposite2D::getNMax() - 1;
	 Decomposite2D::setNMax(newNMAX);
	 double chisquare = Decomposite2D::getChiSquare();
	 double variance = Decomposite2D::getChiSquareVariance();
	 checkCorrelationFunctionFromResiduals();
	 text << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
	 text << chisquare << "\t" << variance << "\t" + comp_corr_string << endl;
	 history.append(text);
	 if (comp_corr >= 0)
	   status = findOptimalBeta(7);
       }
       optimalNMax = Decomposite2D::getNMax();
     }
   }
 }
  // setting up vector -> matrix trafo stuff (for getting coeffs and errors)
  nCoeffs = getNCoeffs(optimalNMax);
  makeNVector(nVector,nCoeffs,optimalNMax);

  if (status == 0) optimized = 1;
  t1 = time(NULL);
  text << "#" << endl << "# Computation time for decomposition: " << t1 - t0 << " seconds" << endl;
  history.append(text);
}


// see Paper III, eq. 39
void OptimalDecomposite2D::findOptimalNMax(unsigned char step) {
  int iter = 1;
  text << "#" << endl << "# Finding optimal decomposition order n_max";
  text << ", beta = " << beta << std::endl;
  text << "# iter.\tn_max\tchi^2\tsigma(chi^2)";
  if (noise_correlated)
    text<<"\txi_res - xi";
  text<< endl;
  history.append(text);

  double chisquare,newChisquare,derivative_chi2,variance;
  chisquare = Decomposite2D::getChiSquare();
  variance = Decomposite2D::getChiSquareVariance();

  text << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
  text << chisquare << "\t" << variance;
  if (noise_correlated) {
    checkCorrelationFunctionFromResiduals();
    text<<"\t\t" + comp_corr_string;
  }
  text << endl;
  history.append(text);

// during the first steps of the search only use even orders for a fast decomposition
  // only when decomposition is too good ( chisquare < 1 - sigma), go back and
  // use also odd nmax.
  int increment;
  if (step == 5) increment = 1;
  else increment = 2;
  
  if (chisquare <= 1) {
    optimalNMax = Decomposite2D::getNMax();
    text << "# Optimal decomposition order n_max = " << optimalNMax;
    text << " already reached" << endl;
    history.append(text);
  }
  while (chisquare > 1) {
    // break during step 2 if order gets too high
    // for a relatively fast findOptimalBeta() here
    if (step == 2 && Decomposite2D::getNMax() == 6) {
      history.append("# Interrupting here for better estimation of beta\n");
      findOptimalBeta(2);
      if (optimalChiSquare <= 1) {
	optimalNMax = Decomposite2D::getNMax();
	text << "# Optimal decomposition order n_max = " << optimalNMax << endl;
	history.append(text);
	break;
      }
      else {
	history.append("#\n# Continuing search for optimal n_max\n");
	text << "# iter.\tn_max\tchi^2\tsigma(chi^2)";
	if (noise_correlated)
	  text<<"\txi_res - xi";
	text<< endl;
	history.append(text);
      }
    }
    
    else if (step == 2 && Decomposite2D::getNMax()%12 == 0) {
      history.append("# Interrupting here for better estimation of beta\n");
      findOptimalBeta(3);
      if (optimalChiSquare <= 1) {
	optimalNMax = Decomposite2D::getNMax();
	text << "# Optimal decomposition order n_max = " << optimalNMax << endl;
	history.append(text);
	break;
      }
      else {
	history.append("#\n# Continuing search for optimal n_max\n");
	text << "# iter.\tn_max\tchi^2\tsigma(chi^2)";
	if (noise_correlated)
	  text<<"\txi_res - xi";
	text<< endl;
	history.append(text);
      }
    }
    
    // reached the limit of nmax
    if (Decomposite2D::getNMax()+increment > nmaxHigh) {
      if (!nmaxTrouble) {
	optimalNMax = Decomposite2D::getNMax();
	flag = 1;
      }
      text << "# Stopping at nmax limit = " << nmaxHigh << endl;
      history.append(text);
      break;
    }

    // increase shapelet order 
    Decomposite2D::setNMax(Decomposite2D::getNMax()+increment);
    newChisquare = Decomposite2D::getChiSquare();
    variance = Decomposite2D::getChiSquareVariance();

    text << "# " << iter << "\t" << Decomposite2D::getNMax() << "\t";
    text << newChisquare << "\t" << variance;
    if (noise_correlated) {
      checkCorrelationFunctionFromResiduals();
      text<<"\t\t" + comp_corr_string;
    }
    text << endl;
    history.append(text);

    // depending on result of chi^2:
    // chisquare is smaller than 1 + sigma: we've reached the goal
    if (newChisquare <= 1) {
      optimalNMax = Decomposite2D::getNMax();
      text << "# Optimal decomposition order n_max = " << optimalNMax << endl;
      history.append(text);
      break;
    }

    // correlation function of residuals goes below correlation function of noise
    if (noise_correlated && comp_corr < 0) {
      history.append("# Stopping here: pixel correlation in residuals becomes less than expected\n");
      optimalNMax = Decomposite2D::getNMax();
      break;
    }

    // don't do this during the refinement procedure when chi^2 was already low
    if (step != 5) {

      // now decomposition get worse, not a good sign
      // save best nmax and chi2
      if (!nmaxTrouble && newChisquare > chisquare) {
	nmaxTrouble = 1;
	bestChiSquare = chisquare;
	optimalNMax = Decomposite2D::getNMax() - increment;
	derivative_chi2 = (newChisquare - chisquare)/increment;
	text << "# chi^2 becomes worse! Saving best fit n_max = " << optimalNMax << " now." << endl;
	history.append(text);
      }
      // if chi2 becomes better again remember new best values
      if (nmaxTrouble && newChisquare >= 0 && newChisquare < bestChiSquare) {
	nmaxTrouble = 0;
	bestChiSquare = newChisquare;
	optimalNMax = Decomposite2D::getNMax();
	text << "# Better n_max = " <<  optimalNMax << " found now. Continuing search for optimal n_max." << endl;
	history.append(text);
      }
      // if increase in chi^2 becomes even bigger: break
      if (nmaxTrouble && (newChisquare - chisquare)/increment > derivative_chi2) {
	text << "# chi^2 becomes increasingly worse. Object underconstrained." << endl;
	text << "# Returning to best fit n_max = " << optimalNMax << endl;
	history.append(text);
	Decomposite2D::setNMax(optimalNMax);
	break;
      }
      // if chi^2 gets negative (nCoeffs >= nPixels): go back to last nmax
      if (newChisquare < 0 || newChisquare == INFINITY) {
	if (!nmaxTrouble) {
	  optimalNMax = Decomposite2D::getNMax() - increment;
	  nmaxTrouble = 1;
	}
	Decomposite2D::setNMax(optimalNMax);
	text << "# nPixels <= nCoeffs! Object underconstrained." << endl;
	text << "# Returning to best fit n_max = " << optimalNMax << endl;
	history.append(text);
	flag = 3;
	break;
      }
      // if theta_min becomes too small -> undersampling
      if (2*optimalBeta/sqrt(Decomposite2D::getNMax()+1.) < 1) {
	if (!nmaxTrouble) {
	  optimalNMax = Decomposite2D::getNMax() - increment;
	  nmaxTrouble = 1;
	}
	Decomposite2D::setNMax(optimalNMax);
	text << "# 2 Theta_min < 1! Image becomes undersampled." << endl;
	text << "# Returning to best fit n_max = " << optimalNMax << endl;
	history.append(text);
	flag = 2;
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

// searches for beta which minimizes chi^2
// beta is contrained by image_dimension/2
// return 0 if minimum has been found, -2 if not, and other numbers for error codes
int OptimalDecomposite2D::findOptimalBeta(unsigned char step) {
  
  // first check if we already did a minimization at this nmax
  map<int,double>::iterator iter = bestBeta.find(Decomposite2D::getNMax());
  if (iter != bestBeta.end()) {
    optimalBeta = beta = iter->second;
    Decomposite2D::setBeta(optimalBeta);
    iter = bestChi2.find(Decomposite2D::getNMax());
    optimalChiSquare = iter->second;
    text << "#" << endl << "# Finding optimal beta";
    text << ", n_max = " << Decomposite2D::getNMax() << endl;
    text<< "# Using already found minimum: chi^2 = " << optimalChiSquare << " at beta = ";
    text << beta;
    if (noise_correlated) {
      checkCorrelationFunctionFromResiduals();
      text<<", xi_res - xi = (" + comp_corr_string +")";
    }
    text << endl; 
    history.append(text);
    return 0;
  }

  // not minimized yet
  else {
    size_t np = 1;
    size_t iter = 0;
    int status;
    double stepsize, size, accuracy;
    double geometric_constraint = image_dimension/(2*sqrt(Decomposite2D::getNMax() +1.0));

    // in step 1) correct beta is completely unknown, therefore define very loose bounds a,b
    // since this probably not be the last call of the function, coarse accuracy is enough
    // in step 2) beta is expected to decrease since nmax is larger than before, but things
    // are still quite uncertain
    // in steps 3-7) accuracy is improving as is the a priori knowledge of beta
    switch (step) {
    case 1: stepsize = 0.3*beta; accuracy = 0.05*beta; break;
    case 2: beta *=0.75; stepsize = 0.2*beta; accuracy = 0.01*beta; break;
    case 3: stepsize = 0.1*beta; accuracy = 0.01*beta; break;
    case 6: stepsize = 0.05*beta; accuracy = 0.01*beta; break;
    case 7: stepsize = 0.05*beta; accuracy = 0.01*beta; break;
    }

    // this ensures beta is within reasonable bounds
    // according to undersampling and orthogonality
    if (betaLow < 0.2) betaLow = 0.2;
    if (betaHigh > geometric_constraint) betaHigh = geometric_constraint;
    if (beta < betaLow) beta = betaLow;
    if (beta > betaHigh) beta = betaHigh;

    text << "#" << endl << "# Finding optimal beta";
    text << ", n_max = " << Decomposite2D::getNMax() << endl;
    text << "# Setting limits: " << betaLow << " <= beta <= " << betaHigh << endl;;
    text << "# Starting with beta = " << beta << endl;
    text << "# iter.\tchi^2\tbeta\tdelta(beta)" << endl;
    history.append(text);

    size_t max_iter = 100;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    parameters params = { *this, betaLow, betaHigh };
    gsl_multimin_function F;

    // define initial vertex vector
    ss = gsl_vector_alloc (np);
    // define primary stepsize
    gsl_vector_set_all (ss,stepsize);

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

      text << "# " << iter << "\t" << s->fval << "\t" << gsl_vector_get (s->x,0) << endl;
      history.append(text);

      if (status == GSL_SUCCESS) {
	optimalBeta = gsl_vector_get (s->x, 0);
	Decomposite2D::setBeta(beta);
	optimalChiSquare = Decomposite2D::getChiSquare();
	text<< "# Converged to minimum: chi^2 = " << optimalChiSquare << " at beta = ";
	text << optimalBeta;
	if (noise_correlated) {
	  checkCorrelationFunctionFromResiduals();
	  text<<", xi_res - xi = (" + comp_corr_string +")";
	}
	text << endl; 
	history.append(text);
      }
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    bestBeta[Decomposite2D::getNMax()] = optimalBeta;
    bestChi2[Decomposite2D::getNMax()] = optimalChiSquare;

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return status;
  }
}

const NumVector<double>& OptimalDecomposite2D::getResiduals() {
  if (!optimized) optimize();
  return Decomposite2D::getResiduals();
}

const NumVector<double>& OptimalDecomposite2D::getModel() {
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
  NumVector<double> corr = xi.getCorrelationFunction(), sigma = xi.getCorrelationError(),
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
 
void OptimalDecomposite2D::getShapeletCoeffs(NumMatrix<double>& coeffs) {
  if (!optimized) optimize();
  if (coeffs.getRows() != optimalNMax+1 || coeffs.getColumns() != optimalNMax+1)
    coeffs = NumMatrix<double>(optimalNMax+1,optimalNMax+1);

  // getting the active coeff vector
  const NumVector<double>& coeffVector = Decomposite2D::getCoeffs();
  // putting it into matrix form
  vectorMapping(coeffVector,coeffs,nVector,nCoeffs);
}

void OptimalDecomposite2D::getShapeletErrors(NumMatrix<double>& errors) {
  if (!optimized) optimize();
  if (errors.getRows() != optimalNMax+1 || errors.getColumns() != optimalNMax+1)
    errors = NumMatrix<double>(optimalNMax+1,optimalNMax+1);

  // getting error vector from noise only
  const NumVector<double>& noiseError = Decomposite2D::getErrors();
  // now get error from uncertainty in beta
  const NumVector<double>& coeffVector = Decomposite2D::getCoeffs();
  NumVector<double> errorVector(nCoeffs);
  getCoeffErrorFromBeta(coeffVector,errorVector);

  // combining the two by Gauss error propagation
  for (int i = 0; i < nCoeffs; i++)
    errorVector(i) = sqrt(errorVector(i)*errorVector(i) + noiseError(i)*noiseError(i));
  // putting it into matrix form
  vectorMapping(errorVector,errors,nVector,nCoeffs);
}

double OptimalDecomposite2D::getOptimalBeta() {
  if (!optimized) optimize();
  return optimalBeta;
}
double OptimalDecomposite2D::getOptimalChiSquare() {
  if (!optimized) optimize();
  return optimalChiSquare;
}

// compute the error on the coeffs coming from uncertainty in beta
// assume delta(beta) = 10% beta at maximum
// to derive 3 sigma errors on coeffs
// FIXME: find good assumptions on beta error?
void OptimalDecomposite2D::getCoeffErrorFromBeta(const NumVector<double>& coeffVector, NumVector<double>& errorVector) {
  ImageTransformation *trafo =  new ImageTransformation();
  NumMatrix<double> betaTrafo(nCoeffs,nCoeffs);
  trafo->makeRescalingMatrix(betaTrafo,beta*1.1,beta,nCoeffs,nVector);
  //getBetaTrafoMatrix(betaTrafo,beta*1.1,beta);
  errorVector = betaTrafo*coeffVector;
  trafo->makeRescalingMatrix(betaTrafo,beta*0.9,beta,nCoeffs,nVector);
  //getBetaTrafoMatrix(betaTrafo,beta*0.90,beta);
  NumVector<double> lowerCoeffs;
  lowerCoeffs = betaTrafo*coeffVector;
  // subtract them to get difference
  errorVector -= lowerCoeffs;
  for (int i=0; i < nCoeffs; i++)
    errorVector(i) = fabs(errorVector(i))/3; // assume the error to be 3 sigma
  delete trafo;
}

char OptimalDecomposite2D::getDecompositionFlag() {
  return flag;
}


// ### REGULARIZATION ###
// this part implements the regularization procedure which minimizes the
// negative flux in the reconstruction

struct reg_parameters {
  Decomposite2D& d;
  History& history;
  const NumVector<double>& reco;
  const NumVector<double>& residual;
  double& beta;
  int& nmax;
  double& f;
  double& chi2;
  double& lambda;
  double& H;
  double& R;
  double& wantedR;
  bool updateCoeffs;
  double betaLow;
  double betaHigh;
};

void computeFluxes(const NumVector<double>& model, double& negFlux, double& posFlux, double& R, double& H) {
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

  double beta = gsl_vector_get(v, 0);
  const NumVector<double>& reco = p->reco;
  const NumVector<double>& residual = p->residual;

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
  
    double posFlux, negFlux;
    computeFluxes(reco,negFlux,posFlux, p->R, p->H);
    p->f = (p->chi2 + (p->lambda)*p->H);
    p->beta = beta;
  } else 
    p->f = INFINITY ; // make f worse if beta is out of bounds
  
  return (p->f);
}

double f_coeffs(const gsl_vector *v, void *params) {
  reg_parameters * p = (reg_parameters *)params;
  NumVector<double>& coeffs = p->d.accessCoeffs();

  // copy coeffs from vector into Decomposite2D
  for (int i=0; i< v->size; i++)
    coeffs(i) = gsl_vector_get(v, i);
  
  // force update of model and residuals
  p->d.updateModelResiduals();
  // now compute chi^2
  double chi2 = p->d.getChiSquare();
  
  double posFlux, negFlux, R, f, H;
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
  p.history.append("#\n# Minimize w.r.t beta\n");
  gsl_vector *x, *ss;
  gsl_multimin_function my_func;
  my_func.f = &f_beta;
  my_func.n = 1;
  my_func.params = &p;
  
  double beta = p.beta;
  x = gsl_vector_alloc (1);
  ss = gsl_vector_alloc (1);
  gsl_vector_set (x, 0, beta);
  gsl_vector_set (ss, 0, 0.05*beta);

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, x->size);
  gsl_multimin_fminimizer_set (s, &my_func, x, ss);

  int iter =0, status;
  std::ostringstream text;
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
	p.history.append("# new minimum found:\n");
      }
  } while (status == GSL_CONTINUE); 
  // update all quantities with new beta
  // this is neccessary, since p contains values from last call, not from best f
  f_beta(s->x,&p);
  text << "# " << iter << ":\t beta = " << gsl_vector_get (s->x, 0)<< "\t f = " << s->fval;
  text << "\t chi^2 = " << p.chi2 <<  "\t R = "<< p.R << std::endl;
  p.history.append(text);

  // clean up
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);
}

void minimize_coeffs(reg_parameters& p) {
  p.history.append("#\n# Minimizing w.r.t. shapelet coefficients\n");
  const NumVector<double>& coeffVector = p.d.getCoeffs();
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
  std::ostringstream text;
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
      	p.history.append("# new minimum found:\n");
      }
  } while (status == GSL_CONTINUE && iter < 1000);
 
  // update coeffs and all other quantities in p
  f_coeffs(s->x,&p);

  text << "# " << iter << ":\t f = " << s->fval << "\t chi^2 = " << p.chi2 << "\t R = " << p.R << std::endl;
  p.history.append(text);
  
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (ss);
}

// append actual regularization results to vector
void OptimalDecomposite2D::appendRegResults(std::vector<regResults>& results, int nmax, double f, double lambda, double beta, double chi2, double R, const NumVector<double>& coeffs) {
  regResults reg = { nmax, f, lambda, beta, chi2, R , coeffs};
  results.push_back(reg);
}

// find nmax where f is minimal in results
int OptimalDecomposite2D::findNMaxofBestF(std::vector<regResults>& results) {
  int maxnum = results.size()-1;
  double bestf = results[maxnum].f;
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
  double bestR = results[0].R;
  int bestIndex = 0;
  for (int i=1; i<=maxnum; i++) {
    if (results[i].R < bestR && results[i].chi2 <= 1) {
      bestR = results[i].R;
      bestIndex = i;
    }
  }
  return bestIndex;
}
  
double OptimalDecomposite2D::regularize(double wantedR) {
  if (!optimized) optimize();

  int nmax;
  double negFlux, posFlux, R, initialR, H, initialH, chi2, initialChi2, beta;
  const NumVector<double>& model = Decomposite2D::getModel();
  const NumVector<double>& residuals = Decomposite2D::getResiduals();
  const NumVector<double>& coeffs = Decomposite2D::getCoeffs();

  computeFluxes(model, negFlux, posFlux, R, H);
  initialR = R;
  initialH = H;
  chi2 = initialChi2 = optimalChiSquare;
  beta = optimalBeta;
  nmax = optimalNMax;
  text << "#" << std::endl << "# Regularization: minimizing negative Flux " << std::endl;
  text << "# Requirement: R = negative_Flux/positive_Flux < " << wantedR << std::endl;
  text << "# Minimizing f = chi^2 + lambda*H, with H = acosh(1+R)" << std::endl;
  text << "# initial R = " << initialR << ", => initial H = " << H << std::endl;
  history.append(text);

  // nothing to do anymore?
  if (wantedR >= R ) {
    text << "# Regularization not neccessary: R already lower than wanted R = " << wantedR << std::endl;
    history.append(text);
  }
  else if (!noise_correlated && optimalChiSquare > 1) {
    text << "# Regularization not appropriate: chi^2 > 1"<< std::endl;
    history.append(text);
  }
  else {
    // lambda is adjusted such that chi2 and R are compatible;
    // for good performance lower lambda by factor 3
    // because model is less strong restricted by data/chi2
    double lambda = 0.3*chi2/H;
    double initialLambda = lambda;
    double f = optimalChiSquare + lambda*H;
    text << "# Starting with lambda = " << lambda << std::endl;
    history.append(text);

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
	history.append("# Regularization takes too long, stopping here!\n");
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
	    text << "# Object too strongly constrained by data, setting nmax = " << nmax << " and lambda = " << lambda << std::endl;
	    history.append(text);
	    par.updateCoeffs = 1;
	  } else {
	    history.append("# Object too strongly constrained by data, but nmax limit reached.\n# Stopping here and reverting to previous best fit values:\n");
	    trouble = 1;
	    break;
	  }
	} else {
	  history.append("# f increased during last increase of nmax: regularization constraints to strict for this object.\n# Reverting to previous best values:\n");
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
	  text << "# Regularization too strong, setting lambda = " << lambda << std::endl;
	  history.append(text);
	  par.updateCoeffs = 1;
	}
	// reg. too weak, problems is given by high R:
	// fix coeffs when searching best beta. It uses the coeffs from the last
	// search for optimal coeffs (w.r.t. f)
	else if (R > wantedR && chi2 <= 1) {
	  lambda *= GSL_MAX(1.5,initialH/fabs(initialH-H));
	  text << "# Regularization too weak, setting lambda = " << lambda << std::endl;
	  history.append(text);
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
    	text << "# Regularization converged: chi^2 = " << optimalChiSquare << ", R = "<< R << std::endl;
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
      text << "# nmax = " << optimalNMax << ", lambda = " << lambda << " -> beta = " << beta << ", f = " << f << ", chi^2 = " << chi2 << ", R = " << R << std::endl;
    }
    text << "#" << endl << "# Computation time for regularization: " << t1 - t0 << " seconds" << std::endl;
    history.append(text);
  }
  return R;
}

const History& OptimalDecomposite2D::getHistory() {
  return history;
}
