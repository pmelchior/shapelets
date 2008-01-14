#include <IO.h>
#include <frame/Point3D.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/timeb.h>
#include <boost/tokenizer.hpp>
#include <fstream>

namespace ublas = boost::numeric::ublas;

fitsfile* openFITSFile(std::string filename, bool write) {
  // check if file alredy exists
  // if not: call createFITSFile()
  bool exists = 0;
  std::ifstream fin;
  fin.open(filename.c_str());
  if(fin.is_open())
    exists = 1;
  fin.close();
  
  int status = 0;
  fitsfile* outfptr;
  if (exists) {
    fits_open_file(&outfptr, filename.c_str(), (int) write, &status);
  } 
  else if (write)
    outfptr = createFITSFile(filename);
  return outfptr;
}

fitsfile* createFITSFile(std::string filename) {
  int status = 0;
  fitsfile *outfptr;
  filename = "!"+filename; // overwrite existing file if necessary
  // create fits file
  fits_create_file(&outfptr,filename.c_str(), &status); 
  return outfptr;
}

int closeFITSFile(fitsfile* fptr) {
  int status = 0;
  fits_close_file(fptr, &status);
  return status;
}

int updateFITSKeywordString(fitsfile *outfptr, std::string keyword, std::string value, std::string comment) {
  int status = 0;
  fits_write_key (outfptr, getFITSDataType(value), const_cast<char *>(keyword.c_str()), const_cast<char *>(value.c_str()), const_cast<char *>(comment.c_str()), &status);
  return status;
}

int appendFITSHistory(fitsfile *outfptr, std::string history) {
  // since it is too long the be saved in one shot
  // split it line by line
  // and for each line insert a HISTORY line in the header
  boost::char_separator<char> sep("\n");
  boost::tokenizer<boost::char_separator<char> > tok(history, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator tok_iter;
  int status = 0, spaces = 8;
  for(tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter) {
    // replace tab by spaces because tab cannot be saved in FITS header card
    // in order to mimick the tab style, an adaptive amount of spaces is inserted
    std::string card = (*tok_iter);
    char tabReplacement(' ');
    while ((card).find("\t") != std::string::npos) {
      int index = (card).find("\t");
      card.replace(index,1,spaces - (index%spaces), tabReplacement);
    }
    fits_write_history (outfptr,const_cast<char*>(card.c_str()), &status);
  }
  return status;
}

int readFITSKeywordString(fitsfile *fptr, std::string key, std::string& val) {
  int status = 0;
  char* comment = NULL;
  char value[FLEN_CARD];
  fits_read_key (fptr,getFITSDataType(val), const_cast<char *>(key.c_str()),&value, comment, &status);
  val = std::string(value);
  return status;
}

int readFITSKeyCards(fitsfile *fptr, std::string key, std::string& value) {
  int status = 0, nkeys, keylen=0;
  char* comment = NULL;
  char card[FLEN_CARD], keyword[FLEN_KEYWORD];
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  for (int i = 1; i <= nkeys; i++) {
    fits_read_record(fptr, i, card, &status);
    fits_get_keyname(card,keyword,&keylen,&status);
    if (std::string(keyword) == key) {
      // append all but the first 8 characters
      value += card + 8;
      value += "\n";
    }
  }
  return status;
}


void addUniformNoise(NumVector<data_t>& data, data_t noisemean, data_t noiselimit) {
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  // seed the random number generator with the actual time in microsecs
  // otherwise possible noise correlations
  unsigned int seeder = time(NULL)*1000; 
  struct timeb tp;
  ftime(&tp);
  seeder += tp.millitm;
  gsl_rng_set(r,seeder);
  for (int i=0; i < data.size(); i++)
    data(i) += noisemean + noiselimit*gsl_rng_uniform (r);
  gsl_rng_free (r);
}

void addGaussianNoise(NumVector<data_t>& data, data_t noisemean, data_t noisesigma) {
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  // seed the random number generator with the actual time in microsecs
  // otherwise possible noise correlations
  unsigned int seeder = time(NULL)*1000; 
  struct timeb tp;
  ftime(&tp);
  seeder += tp.millitm;
  gsl_rng_set(r,seeder);
  for (int i=0; i < data.size(); i++) {
    data(i) += noisemean + gsl_ran_gaussian (r,noisesigma);
  }
  gsl_rng_free (r);
}

void addPoissonianNoise(NumVector<data_t>& data, data_t noisemean) {
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc (T);
  // seed the random number generator with the actual time in microsecs
  // otherwise possible noise correlations
  unsigned int seeder = time(NULL)*1000; 
  struct timeb tp;
  ftime(&tp);
  seeder += tp.millitm;
  gsl_rng_set(r,seeder);
  for (int i=0; i < data.size(); i++) {
    data(i) += noisemean + gsl_ran_gaussian (r,sqrt(fabs(noisemean+data(i))));
  }
  gsl_rng_free (r);
}

void convolveGaussian(const NumVector<data_t>& image, NumVector<data_t>& result, int width,int height) {
  NumVector<data_t> temp = image;
        
  int i,x,y,k0,k;
  int size=width*height;
  int noffset[4];

  noffset[0]=0; noffset[1]=1; noffset[2]=width; noffset[3]=width+1;

  if(result.size() != image.size()) result.resize(image.size());
  result.clear();

  // convolution (without the appropriate normalization to reduce the
  // number of necessary multiplications)
  for(y=0;y<height-1;y++) {
    k0=y*width;
    for(x=0;x<width-1;x++) {
      k=k0+x;
      // 0 is left out because of the memcpy earlier
      for(i=1;i<4;i++) temp[k]+=image[k+noffset[i]];
    }
  }
  k0=(height-1)*width;
  for(i=1;i<width-1;i++) {
    temp[k0+i]*=4.;
  }
  for(i=1;i<height-1;i++) {
    k=i*width+width-1;
    temp[k]*=4.;
  }
  for(y=1;y<height;y++) {
    k0=y*width;
    for(x=1;x<width;x++) {
      k=k0+x;
      for(i=0;i<4;i++) result[k]+=temp[k-noffset[i]];
    }
  }
  // fill up first row and column
  for(x=0;x<width;x++)
    result[x] += temp[x];
  for(y=1;y<height;y++) {
    k0=y*width;
    result[k0] += temp[k0];
  }

  // normalization
  result[0]/=4.; 
  result[size-1]/=16.;
  k0=(height-1)*width;
  for(i=1;i<width-1;i++) {
    result[i]/=4.; 
    result[k0+i]/=16.;
  }
  for(i=1;i<height-1;i++) {
    k=i*width;
    result[k]/=4.; 
    result[k+width-1]/=16.;
  }
  for(y=1;y<height-1;y++) {
    k0=y*width;
    for(x=1;x<width-1;x++) {
      k=k0+x;
      result[k]/=16.;
    }
  }
}

int makeColorMatrix(NumMatrix<unsigned int>& m, std::string colorscheme) {
  // convert colorscheme into appropriate parameters
  // default parameters are "SPECTRUM"
  char scheme = 0;
  int maxcolors = 768;
  char steps = 6;
  if (colorscheme.compare("RED")==0) {
    scheme = 1;
    maxcolors = 512;
    steps = 2;
  }
  if (colorscheme.compare("BLUE")==0) {
    scheme = 2;
    maxcolors = 512;
    steps = 2;
  }
  if (colorscheme.compare("GREEN")==0) {
    scheme = 3;
    maxcolors = 512;
    steps = 2;
  }
  if (colorscheme.compare("GRAY")==0) {
    scheme = 4;
    maxcolors = 256;
    steps = 1;
  }
  if (colorscheme.compare("WARM")==0) {
    scheme = 5;
    maxcolors = 512;
    steps = 4;
  }
  
  int C = maxcolors;
  int N = steps;
  int M = C/N;

  if (m.getRows() != C || m.getColumns() != 3)
    m = NumMatrix<unsigned char>(C,3);

  

  Point3D<int> start, diff, color;
  for (int k=0; k < C; k++) {
    int i = k/M;
    int j = k%M;
    switch (scheme) {
    case 0: // define "SPECTRUM" transitions
      switch (i) {
      case 0:
	if (j==0) {
	  start = Point3D<int>(0,0,0);
	  diff = Point3D<int>(0,0,255);
	  diff -= start;
	}
	break;
      case 1:
	if (j==0) {
	  start = Point3D<int>(0,0,255);
	  diff = Point3D<int>(0,255,255);
	  diff -= start;
	}
	break;
      case 2:
	if (j==0) {
	  start = Point3D<int>(0,255,255);
	  diff = Point3D<int>(0,255,0);
	  diff -= start;
	}
	break;
      case 3:
	if (j==0) {
	  start = Point3D<int>(0,255,0);
	  diff = Point3D<int>(255,255,0);
	  diff -= start;
	}
	break;
      case 4:
	if (j==0) {
	  start = Point3D<int>(255,255,0);
	  diff = Point3D<int>(255,0,0);
	  diff -= start;
	}
	break;
      case 5:
	if (j==0) {
	  start = Point3D<int>(255,0,0);
	  diff = Point3D<int>(255,255,255);
	  diff -= start;
	}
	break;
      }
      break;
    case 1: // "RED"
      switch (i) {
      case 0:
	if (j==0) {
	  start = Point3D<int>(0,0,0);
	  diff = Point3D<int>(255,0,0);
	  diff -= start;
	}
	break;
      case 1:
	if (j==0) {
	  start = Point3D<int>(255,0,0);
	  diff = Point3D<int>(255,255,255);
	  diff -= start;
	}
	break;
      }
      break;
    case 2: // "BLUE"
      switch (i) {
      case 0:
	if (j==0) {
	  start = Point3D<int>(0,0,0);
	  diff = Point3D<int>(64,64,255);
	  diff -= start;
	}
	break;
      case 1:
	if (j==0) {
	  start = Point3D<int>(64,64,255);
	  diff = Point3D<int>(255,255,255);
	  diff -= start;
	}
	break;
      }
      break;
    case 3: // "GREEN"
      switch (i) {
      case 0:
	if (j==0) {
	  start = Point3D<int>(0,0,0);
	  diff = Point3D<int>(0,255,0);
	}
	break;
      case 1:
	if (j==0) {
	  start = Point3D<int>(0,255,0);
	  diff = Point3D<int>(255,255,255);
	  diff -= start;
	}
	break;
      }
      break;
    case 4: // "GRAY"
      if (j==0) {
	start = Point3D<int>(0,0,0);
	diff = Point3D<int>(255,255,255);
	diff -= start;
      }
      break;
    case 5: // "WARM"
      switch (i) {
      case 0:
	if (j==0) {
	  start = Point3D<int>(0,0,0);
	  diff = Point3D<int>(192,0,0);
	  diff -= start;
	}
	break;
      case 1:
	if (j==0) {
	  start = Point3D<int>(192,0,0);
	  diff = Point3D<int>(255,127,0);
	  diff -= start;
	}
	break;
      case 2:
	if (j==0) {
	  start = Point3D<int>(255,127,0);
	  diff = Point3D<int>(255,255,0);
	  diff -= start;
	}
	break;
      case 3:
	if (j==0) {
	  start = Point3D<int>(255,255,0);
	  diff = Point3D<int>(255,255,255);
	  diff -= start;
	}
	break;
      }
      break;  
      break;
    }
    color = diff;
    color *= j;
    color /= M;
    color += start;
    boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<unsigned int> > (m,k) = color;
  }
  return maxcolors;
}

unsigned int getScaledValue(data_t value, int maxcolors, data_t min, data_t max, char scaling) {
  if (value < min) 
    value = min;
  else if (value > max)
    value = max;
  switch(scaling) {
  case 0: // linear
    return (unsigned int) floor((maxcolors-1)*(value-min)/(max-min));
    break;
  case 1: // square root
    return (unsigned int) floor((maxcolors-1)*sqrt((value-min)/(max-min)));
    break;
  case 2: // logarithmic
    data_t value_scaled = (value-min)/(max-min); // between 0 and 1
    data_t offset = 1e-2;
    data_t value_log = (log(value_scaled + offset)-log(offset)) /
      (log(1+offset) - log(offset));
    return (unsigned int) floor((maxcolors-1)*value_log);
    break;
  }
}
void writePPMImage(std::string filename, std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data) {
  FILE *file;
  long x,y,k;
    
  int width = grid.getSize(0);
  int height= grid.getSize(1);

  if((file=fopen(filename.c_str(),"wb"))==NULL) {
    printf("ERROR: writePPMImage: could not open file %s\n",filename.c_str());
    return;
  }
   
  // write the ppm header
  fprintf(file,"P6\n");
  fprintf(file,"%d %d\n",width,height);
  fprintf(file,"255\n");
  
  

  // convert scaling to scale char
  // default is "LINEAR"
  char scale = 0;
  if (scaling.compare("SQUARE_ROOT")==0)
    scale = 1;
  else if (scaling.compare("LOGARITHMIC")==0)
    scale = 2;

  // construct color matrix
  NumMatrix<unsigned int> rgbColors;
  int maxcolors = makeColorMatrix(rgbColors,colorscheme);

  // get a scaled value of data, ranging from 0 to maxcolor
  // read out the three entries of the colormatrix (r,g,b)
  // and insert it at the correct position k
  for(y=0;y<height;y++) {
    for(x=0;x<width;x++) {
      k=x+y*width;
      unsigned int scaled = getScaledValue(data(k),maxcolors,min,max,scale);
      fputc(rgbColors(scaled,0),file);
      fputc(rgbColors(scaled,1),file);
      fputc(rgbColors(scaled,2),file);
    }
  }
  fclose(file);
}

void makeRGBImage(NumMatrix<unsigned int>& rgbImage, std::string colorscheme, std::string scaling, data_t min, data_t max, const Grid& grid, const NumVector<data_t>& data) {
  unsigned int width = grid.getSize(0);
  unsigned int height= grid.getSize(1);
  
  if (rgbImage.getRows() != width*height || rgbImage.getColumns() != 3)
    rgbImage.resize(width*height,3);

  // convert scaling to scale char
  // default is "LINEAR"
  char scale = 0;
  if (scaling.compare("SQUARE_ROOT")==0)
    scale = 1;
  else if (scaling.compare("LOGARITHMIC")==0)
    scale = 2;

  // construct color matrix
  NumMatrix<unsigned int> rgbColors;
  int maxcolors = makeColorMatrix(rgbColors,colorscheme);

  // get a scaled value of data, ranging from 0 to maxcolor
  // read out the three entries of the colormatrix (r,g,b)
  // and insert it at the correct position k
  unsigned int k;
  for(unsigned int y=0;y<height;y++) {
    for(unsigned int x=0;x<width;x++) {
      k=x+y*width;
      unsigned int scaled = getScaledValue(data(k),maxcolors,min,max,scale);
      ublas::matrix_row<ublas::matrix<unsigned int> >(rgbImage,k) = 
	ublas::matrix_row<ublas::matrix<unsigned int> > (rgbColors,scaled);
    }
  }
}

void writeRGB2PPMImage (std::string filename, const Grid& grid, const NumMatrix<unsigned int>& rgbImage) {
  FILE *file;
  unsigned int width = grid.getSize(0);
  unsigned int height= grid.getSize(1);

  if((file=fopen(filename.c_str(),"wb"))==NULL) {
    printf("ERROR: writeRGB2PPMImage: could not open file %s\n",filename.c_str());
    return;
  }
   
  // write the ppm header...
  fprintf(file,"P6\n");
  fprintf(file,"%d %d\n",width,height);
  fprintf(file,"255\n");

  // and now the RGB entries of rgbImage successively
  unsigned int k;
  for(unsigned int y=0;y<height;y++) {
    for(unsigned int x=0;x<width;x++) {
      k=x+y*width;
      fputc(rgbImage(k,0),file);
      fputc(rgbImage(k,1),file);
      fputc(rgbImage(k,2),file);
    }
  }
  fclose(file);
}
