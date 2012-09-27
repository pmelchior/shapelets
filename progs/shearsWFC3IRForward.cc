#include <shapelens/ShapeLens.h>
#include <shapelens/lensing/DEIMOSForward.h>
#include <shapelens/lensing/LensHelper.h>
#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace shapelens;

int N = 2;
int C = 6;
int B = 1;
int L = 1014;
int D = L/B;

data_t scale_list[] = {
  1.0,
  1.25,
  1.5,
  1.75,
  2.0,
  //2.25,
  2.5,
  //2.75,
  3.0,
  //3.5,
  //4.0,
  5.0
};

data_t fiducial_scale = 1.0;

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

int main(int argc, char* argv[]) {
  TCLAP::CmdLine cmd("Measure moments on all WFC3/IR exposures based on detection in the coadd", ' ', "0.4");
  TCLAP::ValueArg<std::string> catfile("c","catalog","Coadd catalog file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> segmapfile("s","segmentation","Coadd segmentation map", true, "","string",cmd);
  TCLAP::ValueArg<std::string> file("f","file","Coadd image file", true, "","string",cmd);
  TCLAP::ValueArg<std::string> expfile("e","exposures","File with list of exposures to consider", true, "", "string", cmd);
  TCLAP::ValueArg<std::string> starfile("S","starfile","Star file for PSF measurement", true, "","string",cmd);
  TCLAP::ValueArg<std::string> outfile("o","outfile","File to write shape catalog", false, "","string",cmd);
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

  // get pointers to exposures as vector
  // and set up the mapping from coadd pixels to exposure pixels
  std::vector<Image<data_t> > exposures;
  std::vector<Image<data_t> > weights;
  std::vector<data_t> noise_rms;
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
      data_t rms;
      IO::readFITSKeyword(fptr, "NOISERMS", rms);
      noise_rms.push_back(rms);
      IO::readFITSImage(fptr, im);
      weights.push_back(im);
      IO::closeFITSFile(fptr);
    }
  }
  int K = exposures.size();

  // define set of scales to be used for the moment measurement
  std::set<data_t> scales;
  int scale_len = sizeof(scale_list) / sizeof(data_t);
  for(int i = 0; i < scale_len; i++)
    scales.insert(scale_list[i]);
  
  // get PSF moments at each of scales from the starfile
  // use individual WCS for each exposure
  // -> PSF moments are different even though PSF image is identical
  Image<data_t> star(starfile.getValue());
  ShapeLensConfig::USE_WCS = true;
  int L_star = star.grid.getSize(0);
  int S_star = 4; // oversampling factor
  std::vector<DEIMOS::PSFMultiScale> mepsfs;
  for (int i=0; i < K; i++) {
    DEIMOS::PSFMultiScale psfs;
    Object psf = star;
    psf.centroid = Point<data_t>(L/2, L/2); // center of exposure
    psf.grid = Grid(L/2 - L_star/2, L/2 - L_star/2, L_star, L_star);

    // WCS needs to account for oversample and exposure information
    ScalarTransformation s(1./S_star);
    s *= *wcsptr[i];
    psf.grid.setWCS(s);

    for (std::set<data_t>::iterator iter = scales.begin(); iter != scales.end(); iter++) {
      DEIMOSElliptical d(psf, N, C, (*iter)*S_star);
      data_t flux = d.mo(0,0);
      d.mo /= flux;
      std::complex<data_t> eps = shapelens::epsilon(d.mo);
      data_t s = sqrt((d.mo(0,2) + d.mo(2,0))/d.mo(0,0));
      s /= d.scale_factor * S_star;
      //std::cout << i << "\t" << *iter << "\t" << d.mo << "\t" << eps << "\t" << s << std::endl;
      psfs.insert(*iter, d.mo);
    } 
    mepsfs.push_back(psfs);
  }


  // extract objects in each exposure
  // measure moments for all available exposures
  std::ofstream fout;
  std::ios_base::fmtflags ff, fd;
  fd = std::cout.flags();
  ff = std::ios_base::fixed;
  if (outfile.isSet()) {
    fout.open(outfile.getValue().c_str());
    fout << "# ID\tRA\t\tDEC\t\tX\tY\tEPS1\tEPS2\tN_EXP\tS_N\tFLAGS" << std::endl;
    fout.flags(ff);
  } else {
    std::cout << "# ID\tRA\t\tDEC\t\tX\tY\tEPS1\tEPS2\tN_EXP\tS_N\tFLAGS" << std::endl;
    std::cout.flags(ff);
  }
  std::map<unsigned long, Moments> moments;
  for (Catalog::const_iterator iter = cat.begin(); iter != cat.end(); iter++) {
    const CatObject& co = iter->second;
    if (1) {//co.FLAGS == 0) { // select only "clean" object"
      std::vector<Object> meo;
      std::vector<DEIMOS::PSFMultiScale> mepsf;
      Rectangle<data_t> bb;
      bb.ll(0) = co.XMIN;
      bb.ll(1) = co.YMIN;
      bb.tr(0) = co.XMAX;
      bb.tr(1) = co.YMAX;
      // increase box size by 5% each side (min. 10 pixels) 
      addFrameBorder(0.05, 10, bb.ll(0), bb.tr(0), bb.ll(1), bb.tr(1));
      for (int i=0; i < K; i++) {
	CoordinateTransformation* im2im = co2exp[i];
	Rectangle<data_t> bb_exp = bb;
	bb_exp.apply(*im2im);
	// check whether object is fully contained in image
	if (bb_exp.ll(0) >= 0 && bb_exp.ll(1) >= 0 && bb_exp.tr(0) < L && bb_exp.tr(1) < L) {
	  Object obj;
	  obj.id = iter->first;
	  Point<int> P1(int(floor(bb_exp.ll(0))), int(floor(bb_exp.ll(1)))),
	    P2(std::min(L-1,int(floor(bb_exp.tr(0)))), std::min(L-1, int(floor(bb_exp.tr(1)))));
	  exposures[i].slice(obj, P1, P2);
	  obj.noise_rms = noise_rms[i];
	  obj.segmentation.resize(obj.size());
	  obj.segmentation.grid = obj.grid;
	  weights[i].slice(obj.weight, P1, P2);
	  // for each pixel in segmap: find out the value in coadd segmap
	  for (long j = 0; j < obj.size(); j++) {
	    Point<data_t> P = obj.grid(j);
	    im2im->inverse_transform(P);
	    obj.segmentation(j) = segmap.get(P);
	    if (obj.weight(j) == 0)
	      obj.segmentation(j) = -1;
	  }
	  obj.centroid = Point<data_t>(co.XCENTROID, co.YCENTROID);
	  im2im->transform(obj.centroid);
	  obj.grid.setWCS(*wcsptr[i]);
	  obj.segmentation.grid.setWCS(*wcsptr[i]);
	  obj.weight.resize(0);
	  meo.push_back(obj);
	  mepsf.push_back(mepsfs[i]);
	}
      }
      if (meo.size() > 0) {
	// measure galactic moments with exposure WCS
	//ShapeLensConfig::VERBOSITY = true;
	DEIMOSForward d(meo, mepsf, N, C, fiducial_scale);
	std::complex<data_t> eps = shapelens::epsilon(d.mo);
	// filter nan's
	if (isnan(real(eps)) || isnan(imag(eps))) {
	  eps = std::complex<data_t>(0,0);
	  d.flags[2] = 1;
	}
	if (isnan(d.SN[fiducial_scale]) || d.SN[fiducial_scale] < 0 || d.SN[fiducial_scale] > 1e7) {
	  d.SN[fiducial_scale] = 0;
	  d.flags[3] = 1;
	}
	// get RA/Dec from centroid, 
	// need to add one since coordinates start with 1 in SExtractor
	Point<data_t> wcs_centroid(co.XCENTROID+1, co.YCENTROID+1);
	wcs.transform(wcs_centroid); // (ra/dec) in coadd
	if (fout.is_open())
	  fout << iter->first << "\t" << std::setprecision(6) << wcs_centroid(0) << "\t" << wcs_centroid(1) << "\t" << std::setprecision(1) << co.XCENTROID << "\t" << co.YCENTROID << "\t" << std::setprecision(3) << real(eps) << "\t" << imag(eps) << "\t" << meo.size() << "\t" << std::setprecision(0) << d.SN[fiducial_scale] << "\t" << d.flags.to_string() << std::endl;
	else
	  std::cout << iter->first << "\t" << std::setprecision(6) << wcs_centroid(0) << "\t" << wcs_centroid(1) << "\t" << std::setprecision(1) << co.XCENTROID << "\t" << co.YCENTROID << "\t" << std::setprecision(3)  << real(eps) << "\t" << imag(eps) << "\t" << meo.size() << "\t" << std::setprecision(0) << d.SN[fiducial_scale] << "\t" << d.flags.to_string() << std::endl;
      }
    }
  }

  // clean up
  exposures.clear();
  weights.clear();
  for (int i=0; i<co2exp.size(); i++) {
    delete co2exp[i];
    delete wcsptr[i];
  }
  if (fout.is_open())
    fout.close();
}
