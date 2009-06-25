#ifdef SHAPELETDB

#include "../../include/shapelets/ShapeletObjectDB.h"
#include <boost/lexical_cast.hpp>
#include <iomanip>

using namespace shapelens;
using namespace std;

#if SHAPELETDB==MySQL

ShapeletObjectDB::ShapeletObjectDB() {
  db.connect("SHAPELENSDBCONF");
  exists = false;
  loadHistory = loadCovariance = true;
}

void ShapeletObjectDB::selectTable(std::string t) {
  table = t;
}

void ShapeletObjectDB::selectDatabase(std::string d) {
  db.selectDatabase(d);
}

void ShapeletObjectDB::save(const ShapeletObjectList& sl) {
  // check whether table exists (if not already tested)
  if (!exists) checkForTable();
  // if it does not exist: create it
  if (!exists) createTable();
  for(ShapeletObjectList::const_iterator iter = sl.begin(); iter != sl.end(); iter++) {
    save(*(*iter));
  }
}

void ShapeletObjectDB::save(const ShapeletObject& sobj) {
  if (!exists) checkForTable();
  if (!exists) createTable();
  ostringstream query;
  //query << setprecision(23); // maximum precision of FLOAT in MySQL
  query << "INSERT INTO `" << table << "` VALUES (" << sobj.getObjectID() << ",";
  query << sobj.getNMax() << "," << sobj.getBeta() << ",";
  query << sobj.getChiSquare() << "," << sobj.getFlags().to_ulong() << ",";
  const Grid& grid = sobj.getGrid();
  query << grid.getStartPosition(0) << "," << grid.getStartPosition(1) << ",";
  query << grid.getSize(0) << "," << grid.getSize(1) << ",";
  const Point2D<data_t>& centroid = sobj.getCentroid();
  query << centroid(0) << "," << centroid(1) << ",";
  query << "'" << sobj.getBaseFilename() << "','";
  // write sobj.prop to query with "\n" lineseparator (see below)
  sobj.prop.write(query,"\\n");
  query << "','";
  // for MySQL we have to convert the newline character (ASCII 10) to the string "\n"
  // therefore we have to copy the history and replace the newlines...
  string hist = sobj.getHistory();
  string newline = "\\n";
  string::size_type loc = hist.find(char(10),0);
  while( loc != string::npos ) {
    hist.replace(loc,1,newline,0,3);
    loc = hist.find(char(10),loc);
  }
  // ... and ' (ASCII 39) if they occur
  loc = hist.find(char(39),0);
  while( loc != string::npos ) {
    hist.erase(loc,1);
    loc = hist.find(char(39),loc);
  }
  query << hist << "',";

  // coefficients (and covariance matrix) have to be converted to hexadecimal values
  const CoefficientVector<data_t>& coeffs = sobj.getCoeffs();
  unsigned long nCoeffs = coeffs.getNCoeffs();
  unsigned long coefflen = nCoeffs * sizeof(data_t);
  //char coeffhex[2*coefflen + 1];
  char* coeffhex = (char*) malloc(2*coefflen + 1);
  mysql_hex_string(coeffhex,reinterpret_cast<const char*>(coeffs.c_array()),coefflen);
  query << "0x" << string(coeffhex) << ",";
  free(coeffhex);

  // same for cov if not empty
  // only store the neccessary entries of the symmetric matrix
  // by copying them row-wise (and shrinking the number of columns);
  // in case of Gaussian errors with diagonal matrix, store only first diagonal
  // element
  const NumMatrix<data_t>& cov = sobj.getCovarianceMatrix();
  if (cov.getRows() > 0) {
    bool diag = true;
    if (ShapeLensConfig::NOISEMODEL == "GAUSSIAN") {
    // check whether cov is diagonal by running along first rows
      for (unsigned int i=0; i < 4; i++) {
	for (unsigned int j=i+1; j < sobj.coeffs.size(); j++) {
	  if (sobj.cov(i,j) != 0) {
	    diag = false;
	    break;
	  }
	}
      }
      // covlen  = 1 here:
      if (diag) {
	char* covhex = (char*) malloc(2*sizeof(data_t) + 1);
	mysql_hex_string(covhex, reinterpret_cast<const char*>(&cov(0,0)),sizeof(data_t));
	query << "0x" << string(covhex) << ")";
	free(covhex);
      }
    }
    if (ShapeLensConfig::NOISEMODEL != "GAUSSIAN" || diag == false) {
      unsigned long covlen = ((nCoeffs + 1)*nCoeffs)/2 * sizeof(data_t);
      char* covhex = (char*) malloc(2*covlen + 1);
      int done = 0;
      for (unsigned int i = 0; i < nCoeffs; i++) {
	mysql_hex_string(covhex + 2*done*sizeof(data_t), reinterpret_cast<const char*>(cov.c_array() + i*(nCoeffs+1)),(nCoeffs-i)*sizeof(data_t));
	done += nCoeffs-i;
      }
      query << "0x" << string(covhex) << ")";
      free(covhex);
    }
  }
  else 
    query << "'')";

  // send SQL query
  db.query(query.str());
}

ShapeletObjectList ShapeletObjectDB::load(std::string where_clause, std::string join_clause) {
  //send SQL query
  string query = "SELECT " + ref("id") + "," + ref("nmax") + "," + ref("beta") + "," +  ref("chi2") + "," + ref("flags") + "," + ref("min_x") + "," + ref("min_y") + "," + ref("size_x") + "," + ref("size_y") + "," + ref("centroid_x") + "," + ref("centroid_y") + "," + ref("basefile") + "," + ref("prop") + ",";
  if (loadHistory)
    query += ref("history") + ",";
    else
      query += "NULL,";
  query += ref("coeffs") + ",";
  if (loadCovariance)
    query += ref("cov") + ",OCTET_LENGTH(" +ref("cov") + ")";
  else
    query += "NULL,0";
  
  query +=" FROM `" + db.getDatabaseName() + "`.`" + table + "`";
  if (join_clause.size() > 0)
    query += " JOIN " + join_clause;

  if (where_clause.size() > 0)
    query += " WHERE " + where_clause;

  DBResult res = db.query(query);
  MYSQL_ROW row;

  // legacy: check whether name, tag is saved or prop
  MYSQL_FIELD* fields = res.getFields();
  bool legacy = false;
  int Nleg = 12; // last field index for legacy block
  for (unsigned int i = 0; i < res.getFieldCount(); i++) {
    std::string f(fields[i].name);
    if (f=="name") {
      legacy = true;
      Nleg = 14;
      break;
    }
  }

  // store results in list
  ShapeletObjectList sl;
  ShapeletObject tmp;
  CoefficientVector<data_t>& cv = tmp.coeffs;
  NumMatrix<data_t>& cov = tmp.cov;

  while (row = res.getRow()) {
    tmp.history.clear();
    tmp.id = boost::lexical_cast<unsigned long>(row[0]);
    cv.setNMax(boost::lexical_cast<unsigned int>(row[1]));
    tmp.setBeta(boost::lexical_cast<data_t>(row[2]));
    tmp.chisquare = boost::lexical_cast<data_t>(row[3]);
    tmp.flags = bitset<16>(boost::lexical_cast<unsigned long>(row[4]));
    tmp.model.grid = Grid(boost::lexical_cast<int>(row[5]),boost::lexical_cast<int>(row[6]),boost::lexical_cast<int>(row[7]),boost::lexical_cast<int>(row[8]));
    tmp.xcentroid(0) = boost::lexical_cast<data_t>(row[9]);
    tmp.xcentroid(1) = boost::lexical_cast<data_t>(row[10]);
    tmp.basefilename = string(row[11]);
    if (legacy) {
      tmp.prop["name"] = string(row[12]);
      tmp.prop["classifier"] = boost::lexical_cast<data_t>(row[13]);
      tmp.prop["tag"] = boost::lexical_cast<data_t>(row[14]);
    } else {
      istringstream is(row[12]);
      tmp.prop.read(is);
    }
    // copy history
    tmp.history.setSilent();
    tmp.history << row[Nleg+1];
    tmp.history.unsetSilent();
    // load coeffs and cov: binary string has to be interpreted as data_t
    unsigned int nCoeffs = cv.getNCoeffs();
    memcpy(cv.c_array(),reinterpret_cast<data_t*>(row[Nleg+2]),nCoeffs*sizeof(data_t));
    // check whether covariance matrix is stored
    if (row[Nleg+3] != NULL) {
      cov = NumMatrix<data_t>(nCoeffs,nCoeffs); // quickly creates a matrix with zeros
      // cov is stored in either symmetric-packed format or compressed
      // into one number (first diagonal element)
      if (atoi(row[Nleg+4]) == sizeof(data_t)) {// single element
	data_t sigma = *reinterpret_cast<data_t*>(row[Nleg+3]);
	for (unsigned int i=0; i < nCoeffs; i++)
	  cov(i,i) = sigma;
      } else {
	// as symmetric matrix is stored in packed format
	// we have to unpack it here again:
	// readout in row-wise fashion
	int done = 0;
	for (unsigned int i = 0; i < nCoeffs; i++) {
	  memcpy(cov.c_array() + i*(nCoeffs+1),reinterpret_cast<data_t*>(row[Nleg+3]) + done,(nCoeffs-i)*sizeof(data_t));
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
    sl.push_back(safePointer);
  }
  return sl;
}

void ShapeletObjectDB::clear() {
  if (!exists) checkForTable();
  if (exists) {
    string query = "TRUNCATE `" + table + "`";
    db.query(query);
  }
}

void ShapeletObjectDB::checkForTable() {
  string query = "SHOW tables";
  DBResult result = db.query(query);

  exists = false;
  while (MYSQL_ROW row = result.getRow()) {
    string current(row[0]);
    if (table == current) {
      exists = true;
      break;
    }
  }
}

DBResult ShapeletObjectDB::query(std::string query) {
  return db.query(query);
}

void ShapeletObjectDB::createTable() {
  string query = "CREATE TABLE `" + table + "` (";
  query += "`id` int(10) unsigned NOT NULL default '0' COMMENT 'object id',";
  query += "`nmax` tinyint(3) unsigned NOT NULL default '0' COMMENT 'shapelet order',";
  query += "`beta` float unsigned NOT NULL default '0' COMMENT 'shapelet scale',";
  query += "`chi2` float unsigned NOT NULL default '0' COMMENT 'goodness of fit',";
  query += "`flags` smallint(5) unsigned NOT NULL default '0' COMMENT 'segmentation and decomposition flags',";
  query += "`min_x` int(10) unsigned NOT NULL default '0' COMMENT 'min(X) in image pixels',";
  query += "`min_y` int(10) unsigned NOT NULL default '0' COMMENT 'min(Y) in image pixels',";
  query += "`size_x` smallint(5) unsigned NOT NULL default '0' COMMENT 'image dimension in x',";
  query += "`size_y` smallint(5) unsigned NOT NULL default '0' COMMENT 'image dimension in y',";
  query += "`centroid_x` float unsigned NOT NULL default '0' COMMENT 'x-position of centroid in image coordinates',";
  query += "`centroid_y` float unsigned NOT NULL default '0' COMMENT 'y-Position of centroid in image coordinates',";
  query += "`basefile` varchar(255) default NULL COMMENT 'object source file',";
  query += "`prop` text COMMENT 'arbitrary object properties',";
  query += "`history` text COMMENT 'object history',";
  query += "`coeffs` blob NOT NULL COMMENT 'shapelet coefficients',";
  query += "`cov` mediumblob COMMENT 'covariance matrix of shapelet coefficients',";
  query += "PRIMARY KEY  (`id`))";
  query += "ENGINE=MyISAM DEFAULT CHARSET=latin1";

  db.query(query);
  exists = true;
}

void ShapeletObjectDB::useHistory(bool use) {
  loadHistory = use;
}

void ShapeletObjectDB::useCovariance(bool use) {
  loadCovariance = use;
}

std::string ShapeletObjectDB::ref(std::string col) {
  return "`"+db.getDatabaseName()+ "`.`" + table + "`.`" + col + "`";
}

#endif // MySQL

#endif // SHAPELETDB


