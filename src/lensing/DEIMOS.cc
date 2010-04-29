#include <shapelens/lensing/DEIMOS.h>
#include <shapelens/utils/IO.h>

namespace shapelens {
  DEIMOS::DEIMOS() : N(0), mo(), id(0), width(0), flags(0) {}

  DEIMOS::DEIMOS(std::string filename) {
    fitsfile* fptr = IO::openFITSFile(filename);
    NumMatrix<data_t> M;
    IO::readFITSImage(fptr,M);
    int N = int(M.getRows()) - 1;
    mo = MomentsOrdered(N);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	mo(n1,n2) = M(n1,n2);
    IO::readFITSKeyword(fptr,"ID",id);
    IO::readFITSKeyword(fptr,"WIDTH",width);
    int f;
    IO::readFITSKeyword(fptr,"FLAGS",f);
    flags = std::bitset<2>(f);
    IO::closeFITSFile(fptr);
  }

  DEIMOS::DEIMOS (const Object& obj, unsigned int N_) :
    N(N_), mo(obj,N_), id(obj.id), width(obj.w.getScale()), flags(0) {
  }

  void DEIMOS::save(std::string filename) const {
    fitsfile* fptr = IO::createFITSFile(filename);
    int N = mo.getOrder();
    NumMatrix<data_t> M(N+1,N+1);
    for(int n1=0; n1 <= N; n1++)
      for(int n2=0; n2 <= N-n1; n2++)
	M(n1,n2) = mo(n1,n2);
    IO::writeFITSImage(fptr,M);
    IO::updateFITSKeyword(fptr,"ID",id);
    IO::updateFITSKeyword(fptr,"WIDTH",width);
    IO::updateFITSKeyword(fptr,"FLAGS",(int)flags.to_ulong());
    IO::closeFITSFile(fptr);
  }

  void DEIMOS::deweight(unsigned int C) {
    Point<data_t> zero(0,0);
    WeightFunction w(width,zero);
    w.setDerivative(-1);
    data_t w_0 = w(zero);
    w.setDerivative(-2);
    data_t w__0 = w(zero);
    w.setDerivative(-3);
    data_t w___0 = w(zero);
    for (int n=0; n <= N; n++) {
      for (int m=0; m <= n; m++) {
	if (C >= 2 && n <= N - 2)
	  mo(m,n-m) -= w_0*(mo(m+2,n-m) + mo(m,n-m+2));
	if (C >= 4 && n <= N - 4)
	  mo(m,n-m) += w__0/2*(mo(m+4,n-m) + 2*mo(m+2,n-m+2) + mo(m,n-m+4));
	if (C >= 6 && n <= N - 6)
	  mo(m,n-m) -= 2*w___0/15*(mo(m+6,n-m) + 3*mo(m+4,n-m+2) + 3*mo(m+2,n-m+4) + mo(m,n-m+6));
      }
    }
    // note in flags
    flags.set(0);
  }

  // helper cclass
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
  
  void DEIMOS::deconvolve(const DEIMOS& psf, unsigned int D) {
    // set up easy access
    BinoMo b(mo), c(psf.mo);
    // no change for monopole: m(0,0) = b(0,0);
    if (D > 0) {
      mo(0,1) -= b(0,0)*c(0,1);
      mo(1,0) -= b(0,0)*c(1,0);
      if (D > 1) {
	mo(0,2) -= b(0,0)*c(0,2) + 2*b(0,1)*c(0,1);
	mo(1,1) -= 0.5*(b(0,0)*c(1,1) + 2*(b(0,1)*c(1,0) + b(1,0)*c(0,1)));
	mo(2,0) -= b(0,0)*c(2,0) + 2*b(1,0)*c(1,0);
	if (D > 2) {
	  mo(0,3) -= b(0,0)*c(0,3) + 3*b(0,1)*c(0,2) + 3*b(0,2)*c(0,1);
	  mo(1,2) -= 1./3*(b(0,0)*c(1,2) + 3*(b(0,1)*c(1,1) + b(1,0)*c(0,2)) + 3*(b(0,2)*c(1,0) + b(1,1)*c(0,1)));
	  mo(2,1) -= 1./3*(b(0,0)*c(2,1) + 3*(b(0,1)*c(2,0) + b(1,0)*c(1,1)) + 3*(b(1,1)*c(1,0) + b(2,0)*c(0,1)));
	  mo(3,0) -= b(0,0)*c(3,0) + 3*b(1,0)*c(2,0) + 3*b(2,0)*c(1,0);
	  if (D > 3) {
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
    complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
    e/= (complex<data_t>(mo(2,0) + mo(0,2)) + 2.*sqrt(complex<data_t>(mo(0,2)*mo(2,0) - gsl_pow_2(mo(1,1)))));
    return e;
  }
  
  complex<data_t> DEIMOS::chi() const {
    complex<data_t> e(mo(2,0) - mo(0,2),2*mo(1,1));
    e/= mo(2,0) + mo(0,2);
    return e;
  }
 
  // HOLICS equations from OMU 2007
  complex<data_t> DEIMOS::zeta() const {
    data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4); // eq. 25
    complex<data_t> zeta(mo(3,0) + mo(1,2),    // eq. 26
			 mo(2,1) + mo(0,3));
    zeta /= xi;
    return zeta;
  }

  complex<data_t> DEIMOS::delta() const {
    data_t xi = mo(4,0) + 2*mo(2,2) + mo(0,4);  // eq. 25
    complex<data_t> delta(mo(3,0) - 3*mo(1,2),  // eq. 27
			  3*mo(2,1) - mo(0,3));
    delta /= xi;
    return delta;
  }

  
  // #### DEIMOSList ####
  DEIMOSList::DEIMOSList() : std::vector<boost::shared_ptr<DEIMOS> >() {}

  DEIMOSList::DEIMOSList(SQLiteDB& sql, std::string table, std::string where) :
    std::vector<boost::shared_ptr<DEIMOS> > () {
    
  }
  void DEIMOSList::save(SQLiteDB& sql, std::string table) const {
    // drop table: faster than deleting its entries
    std::string query = "DROP TABLE IF EXISTS " + table + ";";
    sql.query(query);

    // create it
    query = "CREATE TABLE " + table + "(";
    query += "`id` int NOT NULL PRIMARY KEY,";
    query += "`width` float NOT NULL,";
    query += "`flags` int NOT NULL,";
    query += "`mo` blob NOT NULL);";
    sql.query(query); 
    
    // create prepared statement
    sqlite3_stmt *stmt;
    query = "INSERT INTO `" + table + "` VALUES (?,?,?,?);";
    sql.checkRC(sqlite3_prepare_v2(sql.db, query.c_str(), query.size(), &stmt, NULL));
    for(DEIMOSList::const_iterator iter = DEIMOSList::begin(); iter != DEIMOSList::end(); iter++) {
      const DEIMOS& d = *(*iter);
      sql.checkRC(sqlite3_bind_int(stmt,1,d.id));
      sql.checkRC(sqlite3_bind_double(stmt,2,d.width));
      sql.checkRC(sqlite3_bind_int(stmt,3,d.flags.to_ulong()));
      sql.checkRC(sqlite3_bind_blob(stmt,4,d.mo.c_array(),d.mo.size()*sizeof(data_t),SQLITE_STATIC));
      if(sqlite3_step(stmt)!=SQLITE_DONE)
	throw std::runtime_error("ShapeletObjectDB: insertion failed: " + std::string(sqlite3_errmsg(sql.db)));
      sql.checkRC(sqlite3_reset(stmt));
    }
    sql.checkRC(sqlite3_finalize(stmt));
  }
  

} // end namespace
