#include "../../include/shapelets/ShapeletObjectList.h"
#include <fstream>
#include <list>
#include <cmath>
#include <cstring>

using namespace std;

namespace shapelens {

  ShapeletObjectList::ShapeletObjectList() : vector<boost::shared_ptr<ShapeletObject> >() {
  }

  ShapeletObjectList::ShapeletObjectList(string filename) : vector<boost::shared_ptr<ShapeletObject> >() {
    // open file with list of ShapeletObjects
    ifstream listfile (filename.c_str());
    if (listfile.fail()) {
      cout << "ShapeletObjectList: list file does not exists!" << endl;
      terminate();
    }
    // read in list file
    string sifname;
    while(getline(listfile, sifname)) {
      boost::shared_ptr<ShapeletObject> safePointer (new ShapeletObject(sifname));
      ShapeletObjectList::push_back(safePointer);
    }
  }

#ifdef HAS_SQLiteDB
  ShapeletObjectList::ShapeletObjectList(SQLiteDB& sqlite, std::string table, std::string where, bool loadHistory, bool loadCovariance) : vector<boost::shared_ptr<ShapeletObject> >() {

    //send SQL query
    string query = "SELECT id,nmax,beta,chi2,flags,min_x,min_y,size_x,size_y,centroid_x,centroid_y,basefile,prop,";
    if (loadHistory)
      query += "history,";
    else
      query += "NULL,";
    query += "coeffs,";
    if (loadCovariance)
      query += "cov";
    else
    query += "NULL";
  
    query +=" FROM " + table;

    if (where.size() > 0)
      query += " WHERE " + where;
    query += ";";

    // create prepared statement
    sqlite3_stmt *stmt;
    sqlite3_prepare_v2(sqlite.db, query.c_str(), query.size(), &stmt, NULL);

    // get objects from DB
    ShapeletObject tmp;
    CoefficientVector<data_t>& cv = tmp.coeffs;
    NumMatrix<data_t>& cov = tmp.cov;
    int rc;
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
      tmp.history.clear();
      tmp.id = sqlite3_column_int(stmt,0);
      cv.setNMax(sqlite3_column_int(stmt,1));
      tmp.setBeta(sqlite3_column_double(stmt,2));
      tmp.chisquare = sqlite3_column_double(stmt,3);
      tmp.flags = bitset<16>(sqlite3_column_int(stmt,4));
      tmp.model.grid = Grid(sqlite3_column_int(stmt,5),
			    sqlite3_column_int(stmt,6),
			    sqlite3_column_int(stmt,7),
			    sqlite3_column_int(stmt,8));
      tmp.xcentroid(0) = sqlite3_column_double(stmt,9);
      tmp.xcentroid(1) = sqlite3_column_double(stmt,10);
      tmp.basefilename = string(reinterpret_cast<const char*>(sqlite3_column_text(stmt,11)));

      // set properties
      tmp.prop.clear();
      if (sqlite3_column_bytes(stmt,12)) {
	istringstream is(reinterpret_cast<const char*>(sqlite3_column_text(stmt,12)));
	tmp.prop.read(is);
      }

      // copy history
      if (sqlite3_column_bytes(stmt,13)) {
	tmp.history.setSilent();
	tmp.history << sqlite3_column_text(stmt,13);
	tmp.history.unsetSilent();
      }
	
      // load coeffs and cov
      memcpy(cv.c_array(),
	     reinterpret_cast<const data_t*>(sqlite3_column_blob(stmt, 14)),
	     sqlite3_column_bytes(stmt, 14));

      // check whether covariance matrix is stored
      if (sqlite3_column_bytes(stmt, 15)) {
	unsigned int nCoeffs = cv.getNCoeffs();
	cov = NumMatrix<data_t>(nCoeffs,nCoeffs);
	// cov is stored in either symmetric-packed format or compressed
	// into one number (first diagonal element)
	if (sqlite3_column_bytes(stmt, 15) == sizeof(data_t)) {// single element
	  data_t sigma = *reinterpret_cast<const data_t*>(sqlite3_column_blob(stmt, 15));
	  for (unsigned int i=0; i < nCoeffs; i++)
	    cov(i,i) = sigma;
	} else {
	  // as symmetric matrix is stored in packed format
	  // we have to unpack it here again:
	  // readout in row-wise fashion
	  int done = 0;
	  for (unsigned int i = 0; i < nCoeffs; i++) {
	    memcpy(cov.c_array() + i*(nCoeffs+1),
		   reinterpret_cast<const data_t*>(sqlite3_column_blob(stmt, 15)) + done,
		   (nCoeffs-i)*sizeof(data_t));
	    done += nCoeffs-i;
	  }
	  // reconstruct the upper left from the lower right side of cov
	  for (unsigned int i = 1; i < nCoeffs; i++)
	    for (unsigned int j = 0; j < i; j++)
	      cov(i,j) = cov(j,i);
	}
      }

      // push copy of tmp to sl
      boost::shared_ptr<ShapeletObject> safePointer (new ShapeletObject(tmp));
      ShapeletObjectList::push_back(safePointer);
    }
    sqlite3_finalize(stmt);
  }

  void ShapeletObjectList::save(SQLiteDB& sqlite, std::string table) {
    // create table (if it doesn't exist)
    createTable(sqlite,table);

    // create prepared statement
    sqlite3_stmt *stmt;
    std::string query = "INSERT INTO `" + table + "` VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
    sqlite.checkRC(sqlite3_prepare_v2(sqlite.db, query.c_str(), query.size(), &stmt, NULL));
    for(ShapeletObjectList::const_iterator iter = ShapeletObjectList::begin(); iter != ShapeletObjectList::end(); iter++) {
      const ShapeletObject& sobj = *(*iter);
      sqlite.checkRC(sqlite3_bind_int(stmt,1,sobj.getObjectID()));
      sqlite.checkRC(sqlite3_bind_int(stmt,2,sobj.getNMax()));
      sqlite.checkRC(sqlite3_bind_double(stmt,3,sobj.getBeta()));
      sqlite.checkRC(sqlite3_bind_double(stmt,4,sobj.getChiSquare()));
      sqlite.checkRC(sqlite3_bind_int(stmt,5,sobj.getFlags().to_ulong()));
      const Grid& grid = sobj.getGrid();
      sqlite.checkRC(sqlite3_bind_int(stmt,6,grid.getStartPosition(0)));
      sqlite.checkRC(sqlite3_bind_int(stmt,7,grid.getStartPosition(1)));
      sqlite.checkRC(sqlite3_bind_int(stmt,8,grid.getSize(0)));
      sqlite.checkRC(sqlite3_bind_int(stmt,9,grid.getSize(1)));
      const Point<data_t>& centroid = sobj.getCentroid();
      sqlite.checkRC(sqlite3_bind_double(stmt,10,centroid(0)));
      sqlite.checkRC(sqlite3_bind_double(stmt,11,centroid(1)));
      sqlite.checkRC(sqlite3_bind_text(stmt,12,sobj.getBaseFilename().c_str(), sobj.getBaseFilename().size(),SQLITE_STATIC));

      if (sobj.prop.size() > 0) {
	std::ostringstream pis;
	sobj.prop.write(pis,"\\n");
	sqlite.checkRC(sqlite3_bind_text(stmt,13,pis.str().c_str(),pis.str().size(),SQLITE_STATIC));
      } else
	sqlite.checkRC(sqlite3_bind_null(stmt,13));
    
      std::string hist = sobj.getHistory();
      if (hist.size() > 0) {
	sqlite.checkRC(sqlite3_bind_text(stmt,14,hist.c_str(),hist.size(),SQLITE_STATIC));
      } else
	sqlite.checkRC(sqlite3_bind_null(stmt,14));
    
      const CoefficientVector<data_t>& coeffs = sobj.getCoeffs();
      int nCoeffs = coeffs.size();
      sqlite.checkRC(sqlite3_bind_blob(stmt,15,coeffs.c_array(),nCoeffs*sizeof(data_t),SQLITE_STATIC));

      // same for cov if not empty
      // only store the neccessary entries of the symmetric matrix
      // by copying them row-wise (and shrinking the number of columns);
      // in case of Gaussian errors with diagonal matrix, 
      // store only first diagonal element
      const NumMatrix<data_t>& cov = sobj.getCovarianceMatrix();
      if (cov.getRows() > 0) {
	if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN")
	  sqlite.checkRC(sqlite3_bind_blob(stmt,16,cov.c_array(),sizeof(data_t),SQLITE_STATIC));
	else {
	  unsigned long covlen = ((nCoeffs + 1)*nCoeffs)/2 * sizeof(data_t);
	  data_t* covpacked = (data_t*) malloc(covlen);
	  int done = 0;
	  for (unsigned int i = 0; i < nCoeffs; i++) {
	    memcpy(covpacked + done, cov.c_array() + i*(nCoeffs+1),(nCoeffs-i)*sizeof(data_t));
	    done += nCoeffs-i;
	  }
	  sqlite.checkRC(sqlite3_bind_blob(stmt,16,cov.c_array(),covlen*sizeof(data_t),SQLITE_STATIC));
	  free(covpacked);
	}
      } else
	sqlite.checkRC(sqlite3_bind_null(stmt,16));

      if(sqlite3_step(stmt)!=SQLITE_DONE)
	throw std::runtime_error("ShapeletObjectDB: insertion failed: " + string(sqlite3_errmsg(sqlite.db)));
      sqlite.checkRC(sqlite3_reset(stmt));
    }
    sqlite.checkRC(sqlite3_finalize(stmt));
  }

  void ShapeletObjectList::createTable(SQLiteDB& sqlite, std::string table) {
    string query = "CREATE TABLE IF NOT EXISTS `" + table + "` (";
    query += "`id` int(10) NOT NULL default '0' PRIMARY KEY,";
    query += "`nmax` tinyint(3) NOT NULL default '0',";
    query += "`beta` float NOT NULL default '0',";
    query += "`chi2` float NOT NULL default '0',";
    query += "`flags` smallint(5) NOT NULL default '0',";
    query += "`min_x` int(10) NOT NULL default '0',";
    query += "`min_y` int(10) NOT NULL default '0',";
    query += "`size_x` smallint(5) NOT NULL default '0',";
    query += "`size_y` smallint(5) NOT NULL default '0',";
    query += "`centroid_x` float NOT NULL default '0',";
    query += "`centroid_y` float NOT NULL default '0',";
    query += "`basefile` varchar(255) default NULL,";
    query += "`prop` text,";
    query += "`history` text,";
    query += "`coeffs` blob NOT NULL,";
    query += "`cov` mediumblob)";
    sqlite.query(query);
  }
#endif // HAS_SQLiteDB

  void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta) {
    average(mean,std_mean,beta,NULL);
  }

  void ShapeletObjectList::average(CoefficientVector<data_t>& mean, CoefficientVector<data_t>& std_mean, data_t& beta, data_t (* weightFunction) (ShapeletObject&, void*), void* p) {
    // set up two empty matrices for mean and std_mean
    // as they are easier to resize
    NumMatrix<data_t> meanMatrix(1,1), stdMatrix(1,1);
    beta = 0;
    data_t sum_weights = 0, sum_weights2 = 0;
    int nmax=0, n1, n2;

    // go through all ShapeletObjects
    for (ShapeletObjectList::iterator iter = ShapeletObjectList::begin(); iter != ShapeletObjectList::end(); iter++) {
      const CoefficientVector<data_t>& coeffs = (*iter)->getCoeffs();
      const IndexVector& nVector = coeffs.getIndexVector();
      data_t weight = 1;
      if (weightFunction != NULL)
	weight = (*weightFunction)(*(*iter), p);
      // if new coeff matrix is bigger than current matrices:
      // expand them
      if (coeffs.getNMax() > nmax) {
	nmax = coeffs.getNMax();
	meanMatrix.resize_clear(nmax+1, nmax+1);
	stdMatrix.resize_clear(nmax+1,nmax+1);
      }
      // go thru all coeffs
      for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
	n1 = nVector.getState1(i);
	n2 = nVector.getState2(i);
	meanMatrix(n1,n2) += weight*coeffs(i);
	stdMatrix(n1,n2) += weight*coeffs(i)*coeffs(i);
      }
      beta += weight*((*iter)->getBeta());
      sum_weights += weight;
      sum_weights2 += weight*weight;
    }
    // compute average of beta and all coeffs
    beta /= sum_weights;
  
    // set up mean and std_mean with appropriate nmax
    mean.setNMax(nmax);
    std_mean.setNMax(nmax);
    const IndexVector& nVector = mean.getIndexVector();
    for (unsigned int i=0; i < nVector.getNCoeffs(); i++) {
      n1 = nVector.getState1(i);
      n2 = nVector.getState2(i);
      std_mean(i) = sqrt((stdMatrix(n1,n2)*sum_weights - meanMatrix(n1,n2)*meanMatrix(n1,n2)) / (sum_weights*sum_weights - sum_weights2));
      mean(i) = meanMatrix(n1,n2) / sum_weights;
    }
  }

  ShapeletObjectList ShapeletObjectList::select(bool (* selectionFunction) (ShapeletObject&, void*), void* p) {
    ShapeletObjectList selection;
    for(vector<boost::shared_ptr<ShapeletObject> >::iterator iter = vector<boost::shared_ptr<ShapeletObject> >::begin(); iter != vector<boost::shared_ptr<ShapeletObject> >::end(); iter++)
      if ((*selectionFunction)(*(*iter),p))
	selection.push_back(*iter);				
    return selection;
  }		

} // end namespace
