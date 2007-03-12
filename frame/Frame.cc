#include <Frame.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>

using namespace std;

Frame::Frame(string filename) : FitsImage<double>(filename) {
  text << "# Reading FITS file " << filename << endl;
  text << "# Image properties: size = "<< FitsImage<double>::getSize(0) << "/" << FitsImage<double>::getSize(1) << endl; 
  history.append(text);
  subtractedBG = estimatedBG = foundObjects= 0;
  noise_rms = noise_mean = 0;
  flag = 0;
  numberofObjects = 0;
}

// estimate noise by iterative sigma clipping
void Frame::estimateNoise() {
  flag=0;
  // for GSL sorting functions we have to copy NumVector to double*
  int npixels = FitsImage<double>::getNumberOfPixels();
  double* D = (double*) malloc(npixels*sizeof(double));
  for (int i=0; i < npixels; i++)
    D[i] = (FitsImage<double>::getData())(i);

  // first check left border of the image for pixel value variations
  // if its 0, its a (simulated) image with 0 noise
  // FIXME: background_variance is set to the minmal value above 0
  // Is this corect? 
  if (gsl_stats_variance(D,1,FitsImage<double>::getSize(0)) == 0) {
    noise_mean = 0;
    double min = gsl_stats_max(D,1,npixels);
    for (int i=0; i < npixels; i++)
      if (D[i] < min && D[i] > 0) min = D[i];
    noise_rms = sqrt(min*min);
    flag = 1;
  } 
  else {
    gsl_sort(D, 1, npixels);
    double sigma = gsl_stats_sd(D,1,npixels);
    double median = gsl_stats_median_from_sorted_data (D,1,npixels);

    // sigma clipping here
    int j, jmax = npixels;
    while (1) {
      j=0;
       for (int i = 0; i < jmax; i++ ) {
	 // if data is 3 sigma arround iterative median, keep it
	 // maybe additional constrain: D[i] > 0?
	 if (D[i] > median-3*sigma && D[i] < median+3*sigma)  {
	   D[j] = D[i];
	   j++;
	 }
       }
      // next time only work on the first jmax = j, all others are not sorted
      if (j >= 1) {
	gsl_sort(D, 1, j-1);
 	median = gsl_stats_median_from_sorted_data (D,1,j-1);
 	sigma = gsl_stats_sd(D,1,j-1);
	// no change in the selected pixels: converged
	// this is not as fast as 
	// if (fabs(newmedian-median) < median/npixels)
	if (jmax == j)
	  break;
	else {
	  jmax = j;
	}
      }
      else {
	history.append("Sky background estimation did not converge!\n");
	flag = 3;
      }
    }
    if (flag != 3) {
      // if distribution becomes definitely skewed to the brighter values
      // create histogram of pixel values between median+-3*sigma
      // if the pixels to the right become too bright (+1 error = sqrt(count))
      // set it to the values of the appropriate left bin
      // to symmetrize the histogram around the mean.
      if(gsl_stats_skew(D,1,j-1) > 0.1) {
	int nbins = j/100;
	if (nbins%2 ==1) nbins++;
	gsl_histogram * h = gsl_histogram_alloc (nbins);
	gsl_histogram_set_ranges_uniform (h,median-3*sigma,median+3*sigma);
	for(int i =0; i<j; i++)
	  gsl_histogram_increment (h, D[i]);
	for(int binnr=nbins/2 +1; binnr<nbins; binnr++) {
	  if (gsl_histogram_get(h,binnr) > gsl_histogram_get(h,nbins-binnr)+sqrt(gsl_histogram_get(h,nbins-binnr)))
	    h->bin[binnr] = gsl_histogram_get(h,nbins-binnr);
	}
	noise_mean = gsl_histogram_mean(h);
	noise_rms = gsl_histogram_sigma(h);
	gsl_histogram_free (h);
	flag = GSL_MAX_INT(flag,2);
      }
      // no skewness, image contains enough noise for a solid noise estimation
      // using sigma-clipping only
      else {
	noise_rms = sigma;
	noise_mean = median;
      }
    }
  }
  free(D);
  text << "# Background estimation: mean = " << noise_mean;
  text << ", sigma = " << noise_rms << std::endl;
  history.append(text);
  estimatedBG = 1;
}
unsigned short Frame::getProcessingFlag() {
  return flag;
}

double Frame::getNoiseMean() {
  if (!estimatedBG) estimateNoise();
  return noise_mean;
}

double Frame::getNoiseRMS() {
  if (!estimatedBG) estimateNoise();
  return noise_rms;
}

void Frame::setNoiseMeanRMS(double mean, double rms) {
  estimatedBG = 1;
  noise_mean = mean;
  noise_rms = rms;
  text << "# Noise estimates explicitly set: mean = " << mean << ", rms = " << rms << std::endl;
  history.append(text);
}
  
void Frame::subtractBackground() {
  if (!estimatedBG) estimateNoise();
  if (!subtractedBG) {
    for (int i=0; i < FitsImage<double>::getNumberOfPixels(); i++) {
      (FitsImage<double>::accessData())(i) -= noise_mean;
    }
    subtractedBG = 1;
    text << "# Background subtraction: noise level = " << noise_mean << std::endl;
    history.append(text);
    noise_mean = 0;
  }
}

void Frame::findObjects(unsigned int minPixels, double significanceThresholdSigma, double detectionThresholdSigma) {
  if (!foundObjects && flag != 3) {
    const NumVector<double>& data = FitsImage<double>::getData();
    int counter = 1;
    unsigned int npixels = FitsImage<double>::getNumberOfPixels();
    pixelmap = NumVector<int>(npixels);
    list<int>::iterator iter;

    // check the bounds for the thresholds
    if (significanceThresholdSigma < 1)
      significanceThresholdSigma = 1.5;
    if (detectionThresholdSigma < 1)
      detectionThresholdSigma = 3;
    if (minPixels == 0)
      minPixels = 50;

    // define threshold with global noise values and supplied threshold levels
    if (!estimatedBG) estimateNoise();
    double lowThreshold = noise_mean + significanceThresholdSigma*noise_rms;
    double highThreshold = noise_mean + detectionThresholdSigma*noise_rms;

    // find the maximum pixel in the field
    double maxvalue=0;
    int maxindex = 0;
    for (int i=0; i< npixels; i++) {
      if (data(i) > maxvalue) {
	maxvalue = data(i);
	maxindex = i;
      }
    }
    // there is a detected maximum (not a flat (noisy) image)
    if (data(maxindex) > highThreshold) {
      // search for all pixels connected to the maximum one
      definePixelMap(maxindex,counter,lowThreshold);
      // if maximum has enough pixels define it as object 1
      if (pixellist.size() >= minPixels) {
	text << "# Maximum Object 1 detected with " << pixellist.size() << " significant pixels at (" << maxindex%(FitsImage<double>::getSize(0)) << "/" << maxindex/(FitsImage<double>::getSize(0)) << ")" << std::endl;
	history.append(text);
      // if not, mark the pixels int the active list with -1 (seen, but not big enough)
      } else {
	for(iter = pixellist.begin(); iter != pixellist.end(); iter++ )
	  pixelmap(*iter) = -1;
	counter--;
      }

      // now look for all other objects apart from the brightest one
      for (int i =0; i < npixels; i++) {
	if (data(i) > highThreshold && pixelmap(i) == 0) {
	  counter++;
	  definePixelMap(i,counter,lowThreshold);
	  if (pixellist.size() >= minPixels) {
	    text << "# Object " << counter << " detected with " << pixellist.size() << " significant pixels at (" << i%(FitsImage<double>::getSize(0)) << "/" << i/(FitsImage<double>::getSize(0)) << ")"  << std::endl;
	    history.append(text);
	  }
	  else {
	    for(iter = pixellist.begin(); iter != pixellist.end(); iter++ )
	      pixelmap(*iter) = -1;
	    counter--;
	  }
	}
      }
    } else {
      counter = 0;
    }

    numberofObjects = counter;
    if (numberofObjects >= 1)
      foundObjects = 1;
  }
}

unsigned int Frame::getNumberOfObjects() {
  return numberofObjects;
}

// cut the image to a small region around the object
// and set all pixels to zero than are not related to the image
void Frame::fillObject(Object& O) {
  if (O.getID() > 0 && O.getID() <= numberofObjects) {

    // for artificial noise for area contaminated by different object
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    // define the points (xmin/ymin) and (xmax/ymax) 
    // that enclose the object
    int axsize0, axsize1,xmin, xmax, ymin, ymax;
    xmin = axsize0 = FitsImage<double>::getSize(0);
    xmax = 0;
    ymin = axsize1 = FitsImage<double>::getSize(1);
    ymax = 0;
    // this could be faster if we had a list of pixellists from findObjects() and
    // definePixelMap()
    for(int i=0; i < axsize0*axsize1; i++) {
      if (pixelmap(i) == O.getID()) {
	int x = i%axsize0;
	int y = i/axsize0;
	if (x < xmin) xmin = x;
	if (y < ymin) ymin = y;
	if (x > xmax) xmax = x;
	if (y > ymax) ymax = y;
      }
    }

    O.history = history;
    text << "# Extracting Object " << O.getID() << ", ";
    text << "found in the area (" << xmin << "/" << ymin << ") to (";
    text << xmax << "/" << ymax << ")" << endl;
    O.history.append(text);
    
    // check if outer sizes of the object are identical to the image
    // boundary, since then the objects is cutted 
    if (xmin == 0 || ymin == 0 || xmax == axsize0 -1 || ymax == axsize1 - 1) {
      O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),4));
      O.history.append("# Object cut off at the image boundary!\n");
    }

    // set estimator of scale size to (average object size)/8
    // since theta_max = object_size/2 = sqrt(nmax+1)*beta
    // we assume here nmax = 15
    //scaleSize = (double)(xmax - xmin + ymax - ymin)/(2*8);
    
    // add border around object for including the edges
    addFrameBorder(xmin,xmax,ymin,ymax);
    text << "# Extending the area around object to (" << xmin << "/" << ymin << ") to (";
    text << xmax << "/" << ymax << ")" << endl;
    O.history.append(text);

    // check if object was close to the image boundary so that noise has to be added
    if (xmin < 0 || ymin < 0 || xmax >= axsize0 || ymax >= axsize1) {
      O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),2));
      O.history.append("# Object close to image boundary: Possible cut-off. Extending frame with noise.\n");
    }

    // define new object data set, find 1-sigma noise oscillations with more 
    // than 4 pixels and set their pixelmap flag to -2
    // int the end only object data into new vector of smaller size, the rest will
    // filled up with artificial noise
    NumVector<double>& objdata = O.accessData();

    objdata = NumVector<double>((xmax-xmin+1)*(ymax-ymin+1));
    const NumVector<double>& data = FitsImage<double>::getData();
    list<int>::iterator iter;
    for (int i =0; i < objdata.size(); i++) {
      // old coordinates derived from new pixel index i
      int axis0 = xmax-xmin+1;
      int x = i%axis0 + xmin;
      int y = i/axis0 + ymin;
      int j = x+y*axsize0;

      // if pixel is out of image region, fill noise
      if (x < 0 || y < 0 || x >= axsize0 || y >= axsize1)
	objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
      else {
	// if pixel is noise: check if its starting point of a larger noise oscillation 
	if (pixelmap(j) == 0) {
	  // either positive ...
	  if (data(j) > noise_mean + noise_rms) {
	    refinePixelMap(j,1,xmin,xmax,ymin,ymax);
	    // found no larger oscillation: restore pixelmap
	    if (pixellist.size() <= 8)
	      for(iter = pixellist.begin(); iter != pixellist.end(); iter++)
		pixelmap(*iter) = 0;
	  }
	  else if (data(j) < noise_mean - noise_rms) {
	    refinePixelMap(j,0,xmin,xmax,ymin,ymax);
	    if (pixellist.size() <= 8)
	      for(iter = pixellist.begin(); iter != pixellist.end(); iter++)
		pixelmap(*iter) = 0;
	  }
	}
	// select only object and noise ...
	if (pixelmap(j) == O.getID() || pixelmap(j) == 0)
	  objdata(i) = data(j);
	// but mask different objects and large noise oscillations with noise
	else {
	  objdata(i) = noise_mean + gsl_ran_gaussian (r, noise_rms);
	  if (pixelmap(j) > 0) { // other object in the frame
	    O.setDetectionFlag(GSL_MAX_INT(O.getDetectionFlag(),1));
	    //O.history << "# Another object nearby, but not overlapping." << endl;
	  }
	}
      }
    }
    gsl_rng_free (r);
    
    // now paint objects border into objectMap
    //paintSegmentBorder(xmin,xmax,ymin,ymax);

    // Grid will be changed but not shifted (all pixels stay at their position)
    O.accessGrid() = Grid(xmin,xmax,1,ymin,ymax,1);

    // Fill other quantities into Object
    O.setNoiseMeanRMS(noise_mean,noise_rms);
    O.setNoiseModel("POISSONIAN");
    O.setBaseFilename(FitsImage<double>::getFilename());
    // this calculates flux and centroid;
    O.getFlux();
  } else {
    text.str("# Frame: This Object does not exist!\n");
    O.history.append(text);
    terminate();
  }
}

NumVector<int>& Frame::getObjectMap() {
  return pixelmap;
}

void Frame::definePixelMap(int startpixel, int counter, double threshold) {
  const NumVector<double>& data = FitsImage<double>::getData();
  // get iteratively all connected pixels that are above threshold
  // starting from maxindex and put it in list
  
  pixellist.clear();
  pixellist.push_back(startpixel);
  pixelmap(startpixel) = counter;
  //int iter = 0;
  list<int>::iterator theIterator;
  theIterator = pixellist.begin();
  int& theEnd = pixellist.back();
  int pixelnumber = 0;
  while (pixelnumber < pixellist.size()) {
    int pixel = *theIterator;
    // loop over all direct neighbors and check if they are above threshold
    for (unsigned int dir = 0; dir < 8 ; dir++) {
      int neighbor = neighborpixel(pixel,dir);
      if (neighbor != -1) {
	if (data(neighbor) > threshold) {
	  if (pixelmap(neighbor) == 0) {
	    pixellist.push_back(neighbor);
	    pixelmap(neighbor) = counter;
	  }
	}
      }
    }
    theIterator++;
    pixelnumber++;
  }
}

// this refines the pixelmap within the bounds of the smaller object frame
// it searches for noise oscillations above the 1 sigma level
void Frame::refinePixelMap(int startpixel, bool positive, int xmin, int xmax, int ymin, int ymax) {
  int axsize0 = FitsImage<double>::getSize(0);
  const NumVector<double>& data = FitsImage<double>::getData();
  pixellist.clear();
  pixellist.push_back(startpixel);
  list<int>::iterator theIterator;
  theIterator = pixellist.begin();
  int pixelnumber = 0;
  while (pixelnumber < pixellist.size()) {
    int pixel = *theIterator;
    int x = pixel%axsize0;
    int y = pixel/axsize0;
    if (x > xmin && x < xmax && y > ymin && y < ymax) {
      for (unsigned int dir = 0; dir < 8 ; dir++) {
	int neighbor = neighborpixel(pixel,dir);
	int xneighbor = neighbor%axsize0;
	int yneighbor = neighbor/axsize0;
	if (xneighbor >= xmin && xneighbor <= xmax && yneighbor >= ymin && yneighbor <= ymax) {
	  if ((pixelmap(neighbor) == 0 || pixelmap(neighbor) == -4) && positive && 
	      data(neighbor) > noise_mean+noise_rms) {
	    pixellist.push_back(neighbor);
	    pixelmap(neighbor) = -2;
	  }
	  else if ((pixelmap(neighbor) == 0 || pixelmap(neighbor) == -4) && !positive && 
		   data(neighbor) < noise_mean-noise_rms) {
	    pixellist.push_back(neighbor);
	    pixelmap(neighbor) = -3;
	  }
	}
      }
    }
    theIterator++;
    pixelnumber++;
  }
}


int Frame::neighborpixel(int pixel,unsigned int direction) {
  int axsize0 = FitsImage<double>::getSize(0);
  int axsize1 = FitsImage<double>::getSize(1);
  int x = pixel%axsize0;
  int y = pixel/axsize0;
  int index;
  switch(direction) {
  case 0: 
    if (x<axsize0-1) index = y*axsize0 + x + 1;  // right neighbour
    else index = -1;
    break;
  case 1:
    if (y<axsize1-1 && x<axsize0-1) index = (y+1)*axsize0 + x + 1;  // right bottom
    else index = -1;
    break;
  case 2: 
    if (y<axsize1-1) index = (y+1)*axsize0 + x ;  // bottom
    else index = -1;
    break;
  case 3: 
    if (y<axsize1-1 && x>0) index = (y+1)*axsize0 + x - 1;  // left bottom
    else index = -1;
    break;  
  case 4: 
    if (x>0) index = y*axsize0 + x - 1; // left
    else index = -1;
    break;
  case 5: 
    if (y>0 && x>0) index = (y-1)*axsize0 + x - 1;  // left top
    else index = -1;
    break;   
  case 6: 
    if (y>0) index = (y-1)*axsize0 + x;  // top
    else index = -1;
    break;
  case 7: 
    if (y>0 && x<axsize0-1) index = (y-1)*axsize0 + x + 1;  // right top
    else index = -1;
    break;  
  }
  return index;
}

// now extend to region around the object by
// typically objectsize/4, minimum 8 pixels
void Frame::addFrameBorder(int& xmin, int& xmax, int& ymin, int& ymax) { 
  int xrange, yrange, xborder, yborder;
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  switch(xrange%4) {
  case 1: xmax--; break;
  case 2: xmin++; xmax++; break;
  case 3: xmax++; break;
  }
  switch(yrange%4) {
  case 1: ymax--; break;
  case 2: ymin++; ymax++; break;
  case 3: ymax++; break;
  }
  xrange = xmax - xmin;
  yrange = ymax - ymin;
  // make the object frame square, because of const beta in both directions
  if (xrange < yrange) {
    yborder = GSL_MAX_INT(yrange/4, 10);
    xborder = yborder + (yrange - xrange)/2;
  } else {
    xborder = GSL_MAX_INT(xrange/4, 10);
    yborder = xborder + (xrange - yrange)/2;
  }
  xmin -= xborder;
  xmax += xborder-1;
  ymin -= yborder;
  ymax += yborder-1;
}

void Frame::paintSegmentBorder(int xmin, int xmax, int ymin, int ymax) {
  int axsize0 = FitsImage<double>::getSize(0);
  int axsize1 = FitsImage<double>::getSize(1);
  // check if corners are within frame
  if (xmin<0) xmin=0;
  if (xmax>=axsize0) xmax=axsize0-1;
  if (ymin<0) ymin=0;
  if (ymax>=axsize1) xmax=axsize1-1;
  // low border
  for (int i=(xmin+ymin*axsize0); i<=(xmax+ymin*axsize0); i++)
    if (pixelmap(i) == 0)
      pixelmap(i) = -4;
  // high border
  for (int i=(xmin+ymax*axsize0); i<=(xmax+ymax*axsize0); i++)
    if (pixelmap(i) == 0)
      pixelmap(i) = -4;
  // left border
  for (int i=(xmin+(ymin+1)*axsize0);i<=(xmin+(ymax-1)*axsize0); i+=axsize0)
    if (pixelmap(i) == 0)
      pixelmap(i) = -4;
  // right border
  for (int i=(xmax+(ymin+1)*axsize0);i<=(xmax+(ymax-1)*axsize0); i+=axsize0)
    if (pixelmap(i) == 0)
      pixelmap(i) = -4;
}

const History& Frame::getHistory () {
  return history;
}
