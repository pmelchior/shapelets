#include <shapelens/lensing/DEIMOS.h>
#include <shapelens/utils/IO.h>

namespace shapelens {

  DEIMOS::DEIMOSWeightFunction::DEIMOSWeightFunction(data_t scale, const Point<data_t>& centroid, const complex<data_t>& eps) :
    GaussianWeightFunction(scale,centroid),
    T(0,eps) { 
    // no need to set centroid as it is contained in 
  }

  data_t DEIMOS::DEIMOSWeightFunction::operator()(const Point<data_t>& P_) const {
    Point<data_t> P = P_;
    T.inverse_transform(P);
    return GaussianWeightFunction::operator()(P);
  }

  DEIMOS::DEIMOS() : id(0), flags(0), scale(0), eps(0,0), C(0) {}

  DEIMOS::DEIMOS (Object& obj, int N, int C_, data_t scale_) :
    id(obj.id), flags(0), scale(scale_), eps(0,0), C(C_) {
    //DEIMOSWeightFunction w(scale, obj.centroid, eps);
    //mo = MomentsOrdered(obj,w,N);
    focus(obj,N);
  }

  // DEIMOS::DEIMOS(std::string filename) {
//     fitsfile* fptr = IO::openFITSFile(filename);
//     NumMatrix<data_t> M;
//     IO::readFITSImage(fptr,M);
//     int N = int(M.getRows()) - 1;
//     mo = MomentsOrdered(N);
//     for(int n1=0; n1 <= N; n1++)
//       for(int n2=0; n2 <= N-n1; n2++)
// 	mo(n1,n2) = M(n2,n1);
//     IO::readFITSKeyword(fptr,"ID",id);
//     IO::readFITSKeyword(fptr,"WIDTH",width);
//     int f;
//     IO::readFITSKeyword(fptr,"FLAGS",f);
//     flags = std::bitset<2>(f);
//     IO::closeFITSFile(fptr);
//   }

//   void DEIMOS::save(std::string filename) const {
//     fitsfile* fptr = IO::createFITSFile(filename);
//     int N = mo.getOrder();
//     NumMatrix<data_t> M(N+1,N+1);
//     for(int n1=0; n1 <= N; n1++)
//       for(int n2=0; n2 <= N-n1; n2++)
// 	M(n2,n1) = mo(n1,n2); // transpose to have correctly oriented image
//     IO::writeFITSImage(fptr,M);
//     IO::updateFITSKeyword(fptr,"ID",id);
//     IO::updateFITSKeyword(fptr,"WIDTH",width);
//     IO::updateFITSKeyword(fptr,"FLAGS",(int)flags.to_ulong());
//     IO::closeFITSFile(fptr);
//   }
  
  // determine optimal weighting parameters, centroid, and ellipticity
  void DEIMOS::focus(Object& obj, int N) {
    int iter = 0;
    Point<data_t> centroid_shift;
    while (true) {
      iter++;
      DEIMOSWeightFunction w(scale, obj.centroid, eps);
      mo = MomentsOrdered(obj, w, N + C);
      flags.reset(0);
      flags.reset(1);

      // do initial centroiding and ellipticty estimates
      // from weighted moments: more stable
      centroid_shift(0) = mo(1,0)/mo(0,0);
      centroid_shift(1) = mo(0,1)/mo(0,0);
      obj.centroid += centroid_shift;
      eps = epsilon();
      data_t trQ = mo(2,0) + mo(0,2);
      data_t trQs = trQ*(1-gsl_pow_2(abs(eps)));
      data_t newscale = sqrt(trQs/mo(0,0));

      // deweight now and check whether |ellipticity| < 1
      // if so: use deweighted values, otherwise weighted ones
      MomentsOrdered mo_tmp = mo;
      deweight();
      if (abs(epsilon()) < 0.9 && mo(2,0) > 0 && mo(0,2) > 0) {
	// use new centroid correction, remove old one
	obj.centroid(0) += mo(1,0)/mo(0,0) - centroid_shift(0);
	obj.centroid(1) += mo(0,1)/mo(0,0) - centroid_shift(1);
	eps = epsilon();
	trQ = mo(2,0) + mo(0,2);
	trQs = trQ*(1-gsl_pow_2(abs(eps)));
	scale = sqrt(trQs/mo(0,0));
      } else {
	mo = mo_tmp;
	mo.setOrder(N);
	flags.set(0);
	flags.set(1);
	scale = newscale;
      }
      if (iter >= 20) // stop after 20 iterations
	break;
    }
  }

  void DEIMOS::deweight() {
    // check if already deweighted
    if (!flags.test(0)) {

      data_t c1 = gsl_pow_2(1-real(eps)) + gsl_pow_2(imag(eps));
      data_t c2 = gsl_pow_2(1+real(eps)) + gsl_pow_2(imag(eps));
      data_t e2 = imag(eps);

      data_t s2 = gsl_pow_2(scale);
      data_t s4 = gsl_pow_4(scale);
      data_t s6 = gsl_pow_6(scale);

      int N = mo.getOrder();
      for (int n=0; n <= N; n++) {
	for (int m=0; m <= n; m++) {
	  if (C >= 2 && n <= N - 2) { 
	    mo(m,n-m) += (c1/2*mo(m+2,n-m) - 2*e2*mo(m+1,n-m+1) +
			  c2/2*mo(m,n-m+2))/s2;
	  }
	  if (C >= 4 && n <= N - 4) {
	    mo(m,n-m) += (c1*c1/8*mo(m+4,n-m) - 
			  c1*e2*mo(m+3,n-m+1)  +
			  (c1*c2/4 + 2*e2*e2)*mo(m+2,n-m+2) -
			  c2*e2*mo(m+1,n-m+3) +
			  c2*c2/8*mo(m,n-m+4))/s4;
	  }
	  if (C >= 6 && n <= N - 6) {
	    mo(m,n-m) += (c1*c1*c1/48*mo(m+6,n-m) -
			  c1*c1*e2/4*mo(m+5,n-m+1) +
			  (c1*c1*c2/16 + c1*e2*e2)*mo(m+4,n-m+2) -
			  (c1*c2*e2/2 + 4*e2*e2*e2/3)*mo(m+3,n-m+3) +
			  (c1*c2*c2/16 + c2*e2*e2)*mo(m+2,n-m+4) -
			  c2*c2*e2/4*mo(m+1,n-m+5) +
			  c2*c2*c2/48*mo(m,n-m+6))/s6;
	  }
	}
      }
	
      if (C > 0) {
 	mo.setOrder(N-C);
      }
    
      // note in flags
      flags.set(0);
    }
  }

  // helper class
  class BinoMo {
  public:
    BinoMo(const MomentsOrdered& mo_) : mo(mo_) {
    }
    data_t operator()(int i, int j) const {
      return binomial(i+j,i) * mo(i,j);
    }
  private:
    const MomentsOrdered& mo;
    NumMatrix<unsigned int> binom;
    // CAUTION: this overflows for n > 10!!!
    unsigned long factorial(int n) const {
      unsigned long f = 1;
      for (int m=2; m <= n; m++)
	f *= m;
      return f;
    }
    unsigned long binomial(int n, int m) const {
      return factorial(n)/(factorial(m)*factorial(n-m));
    }
  };
  
  void DEIMOS::deconvolve(const DEIMOS& psf) {
    // set up easy access
    BinoMo b(mo), c(psf.mo);
    int Nmin = std::min(psf.mo.getOrder(),mo.getOrder());
    // no change for monopole: m(0,0) = b(0,0);
    if (Nmin > 0) {
      mo(0,1) -= b(0,0)*c(0,1);
      mo(1,0) -= b(0,0)*c(1,0);
      if (Nmin > 1) {
	mo(0,2) -= b(0,0)*c(0,2) + 2*b(0,1)*c(0,1);
	mo(1,1) -= 0.5*(b(0,0)*c(1,1) + 2*(b(0,1)*c(1,0) + b(1,0)*c(0,1)));
	mo(2,0) -= b(0,0)*c(2,0) + 2*b(1,0)*c(1,0);
	if (Nmin > 2) {
	  mo(0,3) -= b(0,0)*c(0,3) + 3*b(0,1)*c(0,2) + 3*b(0,2)*c(0,1);
	  mo(1,2) -= 1./3*(b(0,0)*c(1,2) + 3*(b(0,1)*c(1,1) + b(1,0)*c(0,2)) + 3*(b(0,2)*c(1,0) + b(1,1)*c(0,1)));
	  mo(2,1) -= 1./3*(b(0,0)*c(2,1) + 3*(b(0,1)*c(2,0) + b(1,0)*c(1,1)) + 3*(b(1,1)*c(1,0) + b(2,0)*c(0,1)));
	  mo(3,0) -= b(0,0)*c(3,0) + 3*b(1,0)*c(2,0) + 3*b(2,0)*c(1,0);
	  if (Nmin > 3) {
	    mo(0,4) -= b(0,0)*c(0,4) + 4*(b(0,1)*c(0,3) + b(0,3)*c(0,1)) + 6*b(0,2)*c(0,2);
	    mo(1,3) -= 1./4*(b(0,0)*c(1,3) + 4*(b(0,1)*c(1,2) + b(1,0)*c(0,3)) + 6*(b(0,2)*c(1,1) + b(1,1)*c(0,2)) + 4*(b(0,3)*c(1,0) + b(1,2)*c(0,1)));
	    mo(2,2) -= 1./6*(b(0,0)*c(2,2) + 4*(b(0,1)*c(2,1) + b(1,0)*c(1,2)) + 6*(b(0,2)*c(2,0) + b(1,1)*c(1,1) + b(2,0)*c(0,2)) + 4*(b(1,2)*c(1,0) + b(2,1)*c(0,1)));
	    mo(3,1) -= 1./4*(b(0,0)*c(3,1) + 4*(b(0,1)*c(3,0) + b(1,0)*c(2,1)) + 6*(b(1,1)*c(2,0) + b(2,0)*c(1,1)) + 4*(b(2,1)*c(1,0) + b(3,0)*c(0,1)));
	    mo(4,0) -= b(0,0)*c(4,0) + 4*(b(1,0)*c(3,0) + b(3,0)*c(1,0)) + 6*b(2,0)*c(2,0);
	  }
	}
      }
    }
    // correction for psf monopole different from unity
    mo /= c(0,0);

    // note in flags
    flags.set(1);

  }

  complex<data_t> DEIMOS::epsilon() const {
    if (mo.getOrder() >= 2) {
      complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
      e/= (complex<data_t>(mo(2,0) + mo(0,2)) + 2.*sqrt(complex<data_t>(mo(0,2)*mo(2,0) - gsl_pow_2(mo(1,1)))));
      return e;
    } else
      return complex<data_t>(0,0);
  }
  
  complex<data_t> DEIMOS::chi() const {
    if (mo.getOrder() >= 2) {
      complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
      e/= mo(2,0) + mo(0,2);
      return e;
    } else
      return complex<data_t>(0,0);
  }
 
  // HOLICS equations from OMU 2007
  complex<data_t> DEIMOS::zeta() const {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4); // eq. 25
      complex<data_t> zeta(mo(3,0) + mo(1,2),    // eq. 26
			   mo(2,1) + mo(0,3));
      zeta /= xi;
      return zeta;
    } else
      return complex<data_t>(0,0);
  }

  complex<data_t> DEIMOS::delta() const {
    if (mo.getOrder() >= 4) {
      data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4);  // eq. 25
      complex<data_t> delta(mo(3,0) - 3*mo(1,2),  // eq. 27
			    3*mo(2,1) - mo(0,3));
      delta /= xi;
      return delta;
    } else
      return complex<data_t>(0,0);
  }

  
//   // #### DEIMOSList ####
//   DEIMOSList::DEIMOSList() : std::vector<boost::shared_ptr<DEIMOS> >() {}

//   DEIMOSList::DEIMOSList(SQLiteDB& sql, std::string table, std::string where) :
//     std::vector<boost::shared_ptr<DEIMOS> > () {
    
//   }
//   void DEIMOSList::save(SQLiteDB& sql, std::string table) const {
//     // drop table: faster than deleting its entries
//     std::string query = "DROP TABLE IF EXISTS " + table + ";";
//     sql.query(query);

//     // create it
//     query = "CREATE TABLE " + table + "(";
//     query += "`id` int NOT NULL PRIMARY KEY,";
//     query += "`width` float NOT NULL,";
//     query += "`flags` int NOT NULL,";
//     query += "`mo` blob NOT NULL);";
//     sql.query(query); 
    
//     // create prepared statement
//     sqlite3_stmt *stmt;
//     query = "INSERT INTO `" + table + "` VALUES (?,?,?,?);";
//     sql.checkRC(sqlite3_prepare_v2(sql.db, query.c_str(), query.size(), &stmt, NULL));
//     for(DEIMOSList::const_iterator iter = DEIMOSList::begin(); iter != DEIMOSList::end(); iter++) {
//       const DEIMOS& d = *(*iter);
//       sql.checkRC(sqlite3_bind_int(stmt,1,d.id));
//       sql.checkRC(sqlite3_bind_double(stmt,2,d.width));
//       sql.checkRC(sqlite3_bind_int(stmt,3,d.flags.to_ulong()));
//       sql.checkRC(sqlite3_bind_blob(stmt,4,d.mo.c_array(),d.mo.size()*sizeof(data_t),SQLITE_STATIC));
//       if(sqlite3_step(stmt)!=SQLITE_DONE)
// 	throw std::runtime_error("ShapeletObjectDB: insertion failed: " + std::string(sqlite3_errmsg(sql.db)));
//       sql.checkRC(sqlite3_reset(stmt));
//     }
//     sql.checkRC(sqlite3_finalize(stmt));
//   }
  

} // end namespace
