#include <shapelens/ShapeLens.h>
#include <shapelens/lensing/DEIMOS.h>
#include <shapelens/lensing/LensHelper.h>
#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>

using namespace shapelens;

int N = 2;
int C = 6;
int B = 10;
int L = 1014;
int D = L/B;

data_t scale_list[] = {
  1.5,
  2.0,
  2.25,
  2.5,
  2.75,
  3.0,
  3.5,
  4.0,
  //  5.0
};

// now extend to region around the object by
// typically objectsize/2, minimum 10 pixels
void addFrameBorder(data_t factor, int minsize, data_t& xmin, data_t& xmax, data_t& ymin, data_t& ymax) { 
  if (factor > 0) {
    data_t xrange, yrange, xborder, yborder;
    xrange = xmax - xmin;
    yrange = ymax - ymin;

    if (xrange < yrange) {
      yborder = std::max(yrange*factor, data_t(minsize));
      xborder = yborder + (yrange - xrange)/2;
    } else {
      xborder = std::max(xrange*factor, data_t(minsize));
      yborder = xborder + (xrange - yrange)/2;
    }
    xmin = xmin - xborder;
    xmax = xmax + xborder;
    ymin = ymin - yborder;
    ymax = ymax + yborder;
  }
}

std::map<Point<int>, DEIMOS::PSFMultiScale> getPSFMultiScaleMap(const std::string& starfile, const std::set<data_t>& scales, const CoordinateTransformation* wcs) {
  std::map<Point<int>, DEIMOS::PSFMultiScale> psfmap;
  fitsfile *fptr = IO::openFITSFile(starfile);
  Object star;
  std::complex<data_t> pos;
  int maxhdu, status = 0;
  fits_get_num_hdus(fptr, &maxhdu, &status);
  for (int i = 1; i <= maxhdu; i++) {
    // read star data and weight map
    IO::moveToFITSExtension(fptr, i);
    IO::readFITSImage(fptr, star);
    IO::readFITSKeyword(fptr, "POSITION", pos);
    star.centroid(0) = real(pos);
    star.centroid(1) = imag(pos);
    IO::readFITSKeyword(fptr, "POS_MIN", pos);
    star.centroid(0) -= real(pos);
    star.centroid(1) -= imag(pos);
    star.grid.setWCS(*wcs);
    i++;
    IO::moveToFITSExtension(fptr, i);
    IO::readFITSImage(fptr, star.weight);
    // define cell
    Point<int> P;
    P(0) = std::min(B-1, int(floor(star.centroid(0) + real(pos)))/D); // actual position
    P(1) = std::min(B-1, int(floor(star.centroid(1) + imag(pos)))/D);
    DEIMOS::PSFMultiScale& psfms = psfmap[P];
    // measure PSF moments at all scales
    for (std::set<data_t>::const_iterator iter = scales.begin(); iter != scales.end(); iter++) {
      ShapeLensConfig::USE_WCS = true;
      DEIMOS psf(star, N, C, *iter);
      ShapeLensConfig::USE_WCS = false;
      data_t flux = psf.mo(0,0);
      psf.mo /= flux; // flux normalization
      Moments& mo = psfms[*iter];
      if (mo.size() > 0)
	mo += psf.mo;
      else
	mo = psf.mo;
    }
  }
  
//   // go into all cells and check for missing data
//   Point<int> P;
//   for (P(0)=0; P(0) < B; P(0)++) {
//     for (P(1)=0; P(1) < B; P(1)++) {
//       DEIMOS::PSFMultiScale& psfms = psfmap[P];
//       std::cout << P(0) << "\t" << P(1);
//       for (std::set<data_t>::iterator iter = scales.begin(); iter != scales.end(); iter++) {
// 	Moments& mo = psfms[*iter];
// 	std::complex<data_t> eps = epsilon(mo);
// 	std::cout << "\t" << real(eps) << "\t" << imag(eps);
//       }
//       std::cout << std::endl;
//     }
//   }

  return psfmap;
}


int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Measure moments on all WFC3/IR exposures based on detection in the coadd", ' ', "0.3");
  TCLAP::ValueArg<std::string> catfile("c","catalog","Coadd catalog file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> segmapfile("s","segmentation","Coadd segmentation map", true, "","string",cmd);
  TCLAP::ValueArg<std::string> file("f","file","Coadd image file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> expfile("e","exposures","File with list of exposures to consider", true, "", "string", cmd);
  TCLAP::ValueArg<std::string> starfile("S","starfile","Star file for PSF measurement", true, "","string",cmd);
  TCLAP::ValueArg<std::string> outfile("o","output","Filename for shear catalog", true, "", "string", cmd);
  cmd.parse(argc,argv);


  // open coadd, catalog and segmap
  fitsfile* fptr = IO::openFITSFile(file.getValue());
  WCSTransformation wcs(fptr, false);
  IO::closeFITSFile(fptr);
  std::list<std::string> optional;
  optional.push_back("FWHM_IMAGE");
  optional.push_back("CLASS_STAR");
  Catalog cat(catfile.getValue(), optional);
  SegmentationMap segmap(segmapfile.getValue());

  //get pointers to exposures as vector
  // and set up the mapping from coadd pixels to exposure pixels
  std::vector<Image<data_t> > exposures;
  std::vector<Image<data_t> > weights;
  std::vector<CoordinateTransformation*> co2exp;
  std::vector<CoordinateTransformation*> wcsptr;
  {
    std::ifstream ifs;
    ifs.open(expfile.getValue().c_str());
    std::string line;
    Image<data_t> im(L,L);
    while (ifs >> line) {
      fptr = IO::openFITSFile(line);
      WCSTransformation wcs_exp(fptr, false);
      wcs_exp.invert(); // sky->xy
      WCSTransformation* tmp = new WCSTransformation(wcs);
      *tmp *= wcs_exp;
      co2exp.push_back(tmp);
      
      WCSTransformation* wcs_intermediate = new WCSTransformation(fptr);
      wcsptr.push_back(wcs_intermediate);
      IO::readFITSImage(fptr, im);
      exposures.push_back(im);
      IO::moveToFITSExtension(fptr, "WEIGHT");
      IO::readFITSImage(fptr, im);
      weights.push_back(im);
      IO::closeFITSFile(fptr);
    }
  }

  // extract objects in each exposure
  std::vector<std::list<Object> > meo;
  for (int i=0; i < exposures.size(); i++) {
    CoordinateTransformation* im2im = co2exp[i];
    std::list<Object> objs;
    for (Catalog::const_iterator iter = cat.begin(); iter != cat.end(); iter++) {
      const CatObject& co = iter->second;
      if (co.FLAGS == 0) { // select only "clean" object"
	Rectangle<data_t> bb;
	bb.ll(0) = co.XMIN;
	bb.ll(1) = co.YMIN;
	bb.tr(0) = co.XMAX;
	bb.tr(1) = co.YMAX;
	// increase box size by 15% each side (min. 10 pixels) 
	addFrameBorder(0.15, 10, bb.ll(0), bb.tr(0), bb.ll(1), bb.tr(1));
	bb.apply(*im2im);
	// check whether object is fully contained in image
	if (bb.ll(0) >= 0 && bb.ll(1) >= 0 && bb.tr(0) < L && bb.tr(1) < L) {
	  Object obj;
	  obj.id = iter->first;
	  Point<int> P1(int(floor(bb.ll(0))), int(floor(bb.ll(1)))),
	    P2(std::min(L-1,int(floor(bb.tr(0)))), std::min(L-1, int(floor(bb.tr(1)))));
	  exposures[i].slice(obj, P1, P2);
	  weights[i].slice(obj.weight, P1, P2);
	  //obj.segMap.grid = obj.grid;
	  //obj.segMap.resize(obj.grid.size());
	  // for each pixel in segmap: find out the value in coadd segmap
	  for (long j = 0; j < obj.size(); j++) {
	    Point<data_t> P = obj.grid(j);
	    im2im->inverse_transform(P);
	    long segmapval = segmap.get(P);
	    // obj.segMap(j) = segmapval;
	    // set weight of pixels from overlapping objects to 0
	    if (segmapval != 0 && segmapval != obj.id)
	      obj.weight(j) = 0;
	    else if (segmapval == 0)
	      obj.noise_rms = 1./sqrt(obj.weight(j));
	  }
	  obj.centroid = Point<data_t>(co.XCENTROID, co.YCENTROID);
	  im2im->transform(obj.centroid);
	  obj.grid.setWCS(*wcsptr[i]);
	  objs.push_back(obj);
	}
      }
    }
    meo.push_back(objs);
  }
  exposures.clear(); // not needed anymore
  weights.clear();


  // define set of scales to be used for the moment measurement
  std::set<data_t> scales;
  int scale_len = sizeof(scale_list) / sizeof(data_t);
  for(int i = 0; i < scale_len; i++)
    scales.insert(scale_list[i]);
  
  // for each exposure:
  // 1) construct PSF map from starfile with exposure WCS
  // 2) measure galactic moments with exposure WCS
  // 3) add galactic moments based on ID
  std::map<unsigned long, Moments> moments;
  for (int i=0; i < meo.size(); i++) {
    std::cout << "############" << std::endl;
    std::cout << "# Exposure " << i << std::endl;
    std::cout << "############" << std::endl;
    // get PSF moment maps
    std::map<Point<int>, DEIMOS::PSFMultiScale> psfmap = getPSFMultiScaleMap(starfile.getValue(), scales, wcsptr[i]);
    
    // iterate thru all (galaxy) objects
    const std::list<Object>& objects = meo[i];
    for (std::list<Object>::const_iterator iter = objects.begin(); 
	 iter != objects.end(); iter++) {
      const Object& obj = *iter;

      // get psf multiscale for this point
      Point<int> P;
      P(0) = std::min(B-1, int(floor(obj.centroid(0)))/D);
      P(1) = std::min(B-1, int(floor(obj.centroid(1)))/D);
      // measure deconvolved moments
      ShapeLensConfig::USE_WCS = true;
      DEIMOS::FIX_CENTROID = true;
      DEIMOS d(obj, psfmap[P], N, C, psfmap[P].getMaximumScale());
      ShapeLensConfig::USE_WCS = false;
      DEIMOS::FIX_CENTROID = false;
      std::cout << obj.id << "\t" << obj.centroid(0) << "\t" << obj.centroid(1) << "\t" << d.mo << "\t";
      std::cout << real(d.epsilon()) << "\t" << imag(d.epsilon());
      std::cout << "\t" << d.scale << "\t" << d.matching_scale << "\t" << d.SN[d.matching_scale] << "\t" << d.flags.to_string() << std::endl;
      
      // if measurement is good, add moments for this object based on id
      if (d.flags.none()) {
	std::map<unsigned long, Moments>::iterator miter = moments.find(obj.id);
	if (miter != moments.end()) {
	  miter->second += d.mo;
	  cat[obj.id].OPT["EXPCOUNT"] = boost::get<int>(cat[obj.id].OPT["EXPCOUNT"]) + 1;
	  cat[obj.id].OPT["TOTAL S/N"] = boost::get<data_t>(cat[obj.id].OPT["TOTAL S/N"]) + d.SN[d.matching_scale];
	}
        else {
	  moments[obj.id] = d.mo;
	  cat[obj.id].OPT["EXPCOUNT"] = 1;
	  cat[obj.id].OPT["TOTAL S/N"] = d.SN[d.matching_scale];
	}
      }
    }
  }

  // output moment/shear catalog
  std::ofstream ofs;
  ofs.open(outfile.getValue().c_str());
  for (std::map<unsigned long, Moments>::iterator miter = moments.begin();
       miter != moments.end(); miter++) {
    Moments& mo = miter->second;
    CatObject& co = cat[miter->first];
    std::complex<data_t> eps = shapelens::epsilon(mo);
    ofs << miter->first << "\t" << co.XCENTROID << "\t" << co.YCENTROID << "\t" << real(eps) << "\t" << imag(eps) << "\t" << boost::get<int>(co.OPT["EXPCOUNT"]) << "\t" << boost::get<data_t>(co.OPT["TOTAL S/N"]) << std::endl;
  }
  ofs.close();

  // clean up
  for (int i=0; i<co2exp.size(); i++) {
    delete co2exp[i];
    delete wcsptr[i];
  }
}
