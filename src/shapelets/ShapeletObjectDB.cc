#ifdef SHAPELETDB

#include <shapelets/ShapeletObjectDB.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fenv.h>

using namespace std;

#if SHAPELETDB==MySQL

ShapeletObjectDB::ShapeletObjectDB() {
  // open file
  ifstream connfile (getenv("SHAPELENSDBCONF"));
  if (connfile.fail()) {
    cerr << "ShapeletObjectDB: connection file could not be opened!" << endl;
    terminate();
  }
  // read in config file
  string line;
  string host, user, password, database;
  while(getline(connfile, line)) {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    // split entries at tabs
    boost::char_separator<char> sep("\t");
    Tok tok(line, sep);
    // first of all we copy the token into string vector
    // though this is not too sophisticated it does not hurt and is more 
    // convenient
    std::vector<std::string> column;
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
      column.push_back(*tok_iter);
    // exclude empty and comment lines
    if (column.size() >= 2 && column[0] != "#") {
      if (column[0] == "HOST")
	host = column[1];
      if (column[0] == "USER")
	user = column[1];
      if (column[0] == "PASSWORD")
	password = column[1];
      if (column[0] == "DATABASE")
	database = column[1];
      if (column[0] == "TABLE")
	table = column[1];
    }
  }
  checkConnectionDetails(host,user,password,database,table);
  db.connect(host,user,password,database);
  exists = false;
}

void ShapeletObjectDB::checkConnectionDetails(std::string host, std::string user, std::string password, std::string database, std::string table) {
  // all keywords (apart from password) must be specified
  if (host.size() == 0) {
    cerr << "ShapeletObjectDB: HOST keyword not specified!";
    terminate();
  }
  if (user.size() == 0) {
    cerr << "ShapeletObjectDB: USER keyword not specified!";
    terminate();
  }
}

void ShapeletObjectDB::selectTable(std::string t) {
  table = t;
}

void ShapeletObjectDB::selectDatabase(std::string database) {
  db.query("USE `" + database + "`");
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
  query << setprecision(23); // maximum precision of FLOAT in MySQL
  query << "INSERT INTO `" << table << "` VALUES (" << sobj.getObjectID() << ",";
  //query << "INSERT INTO `" << table << "` VALUES ('',";
  query << sobj.getNMax() << "," << sobj.getBeta() << ",";
  query << sobj.getChiSquare() << "," << sobj.getFlags().to_ulong() << ",";
  const Grid& grid = sobj.getGrid();
  query << grid.getStartPosition(0) << "," << grid.getStartPosition(1) << ",";
  query << grid.getSize(0) << "," << grid.getSize(1) << ",";
  const Point2D<data_t>& centroid = sobj.getCentroid();
  query << centroid(0) << "," << centroid(1) << ",";
  query << "'" << sobj.getBaseFilename() << "','" << sobj.getName() << "',";
  query << sobj.getObjectClassifier() << "," << sobj.getTag() << ",'";
  // for MySQL we have to convert the newline character (ASCII 10) to the string "\n"
  // therefore we have to copy the history and replace the newlines
  string hist = sobj.getHistory();
  string newline = "\\n";
  string::size_type loc = hist.find(char(10),0);
  while( loc != string::npos ) {
    hist.replace(loc,1,newline,0,3);
    loc = hist.find(char(10),loc);
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

ShapeletObjectList ShapeletObjectDB::load(std::string where_clause) {
  //send SQL query
  string query = "SELECT *, OCTET_LENGTH(`cov`) FROM `" + table + "`";
  if (where_clause.size() > 0) {
    query += " WHERE " + where_clause;
  }
  DBResult res = db.query(query);
  MYSQL_ROW row;

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
    // this assumes integer grids
    tmp.model.grid = Grid(boost::lexical_cast<grid_t>(row[5]),boost::lexical_cast<grid_t>(row[6]),boost::lexical_cast<grid_t>(row[7]),boost::lexical_cast<grid_t>(row[8]));
    tmp.xcentroid(0) = boost::lexical_cast<data_t>(row[9]);
    tmp.xcentroid(1) = boost::lexical_cast<data_t>(row[10]);
    tmp.basefilename = string(row[11]);
    tmp.name = string(row[12]);
    tmp.classifier = boost::lexical_cast<data_t>(row[13]);
    tmp.tag = boost::lexical_cast<data_t>(row[14]);
    // copy history
    tmp.history.setSilent();
    tmp.history << row[15];
    tmp.history.unsetSilent();
    // load coeffs and cov: binary string has to be interpreted as data_t
    unsigned int nCoeffs = cv.getNCoeffs();
    memcpy(cv.c_array(),reinterpret_cast<data_t*>(row[16]),nCoeffs*sizeof(data_t));
    // check whether covariance matrix is stored
    if (row[17] != NULL) {
      cov = NumMatrix<data_t>(nCoeffs,nCoeffs); // quickly creates a matrix with zeros
      // cov is stored in either symmetric-packed format or compressed
      // into one number (first diagonal element)
      if (atoi(row[18]) == sizeof(data_t)) {// single element
	data_t sigma = *reinterpret_cast<data_t*>(row[17]);
	for (unsigned int i=0; i < nCoeffs; i++)
	  cov(i,i) = sigma;
      } else {
	// as symmetric matrix is stored in packed format
	// we have to unpack it here again:
	// readout in row-wise fashion
	int done = 0;
	for (unsigned int i = 0; i < nCoeffs; i++) {
	  memcpy(cov.c_array() + i*(nCoeffs+1),reinterpret_cast<data_t*>(row[17]) + done,(nCoeffs-i)*sizeof(data_t));
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
  query += "`name` varchar(255) default NULL COMMENT 'arbitrary object name',";
  query += "`classifier` float default NULL COMMENT 'object classification number',";
  query += "`tag` float default NULL COMMENT 'arbitrary object tag',";
  query += "`history` text COMMENT 'object history',";
  query += "`coeffs` blob NOT NULL COMMENT 'shapelet coefficients',";
  query += "`cov` mediumblob COMMENT 'covariance matrix of shapelet coefficients',";
  query += "PRIMARY KEY  (`id`))";
  query += "ENGINE=MyISAM DEFAULT CHARSET=latin1";

  db.query(query);
  exists = true;
}

#endif // MySQL


// #if SHAPELETDB=="Berkeley"
// ShapeletObjectDB::ShapeletObjectDB(std::string dbfilename):
//   db(NULL,0), changed(0) {
//   try {
//     // Redirect debugging information to std::cerr
//     db.set_error_stream(&std::cerr);
//     // Open the database
//     db.open(NULL, dbfilename.c_str(), NULL, DB_BTREE, DB_CREATE, 0);
//     // Open cursor
//     db.cursor(NULL, &cursorp, 0);
//   }
//   // DbException is not a subclass of std::exception, so we
//   // need to catch them both.
//   catch(DbException &e) {
//     std::cerr << "ShapeletObjectDB: Error opening database: " << dbfilename << std::endl;
//     std::cerr << e.what() << std::endl;
//   }
//   catch(std::exception &e) {
//     std::cerr << "ShapeletObjectDB: Error opening database: " << dbfilename << std::endl;
//     std::cerr << e.what() << std::endl;
//   }
// }

// ShapeletObjectDB::~ShapeletObjectDB() {
//   // Close the db
//   try {
//     cursorp->close();
//     db.close(0);
//   } 
//   catch(DbException &e) {
//     std::cerr << "ShapeletObjectDB: Error closing database" << std::endl;
//     std::cerr << e.what() << std::endl;
//   }
//   catch(std::exception &e) {
//     std::cerr << "ShapeletObjectDB: Error closing database" << std::endl;
//     std::cerr << e.what() << std::endl;
//   }
// }

// ShapeletObjectList ShapeletObjectDB::all() {
//   ShapeletObjectList sl;
//   try {
//     // Iterate over the database, from the first record to the last
//     Dbt key, data;
//     int ret;
//     while ((ret = cursorp->get(&key, &data, DB_NEXT)) == 0 ) {
//       std::cout << *(u_int32_t*)key.get_data() << "\t";
//       // load all data into ShapeletObject
//       fillShapeletObject(data);
//       // construct new ShapeletObject entity on heap and append it to list
//       boost::shared_ptr<ShapeletObject> safePointer (new ShapeletObject(tmpsobj));
//       sl.push_back(safePointer);
//     }
//     db.cursor(NULL, &cursorp, 0);
//   } catch(DbException &e) {
//     db.err(e.get_errno(), "ShapeletObjectDB: Error in all()");
//     throw e;
//   } catch(std::exception &e) {
//     throw e;
//   }
//   return sl;
// } 

// ShapeletObject ShapeletObjectDB::get(unsigned int i) {
//   try {
//     // Store index i in key, but leave data empty (only search keys)
//     Dbt key(&i,sizeof(u_int32_t));
//     Dbt data;
//     int ret = cursorp->get(&key, &data, DB_SET);
//     // match found
//     if (!ret) {
//       // load all data into ShapeletObject
//       fillShapeletObject(data);
//       return tmpsobj;
//     }
//     else {
//       std::cerr << "ShapeletObjectDB: Key " << i << " does not exist." << std::endl;
//       // FIXME: what to return here?
//     }
//     db.cursor(NULL, &cursorp, 0);
//   } catch(DbException &e) {
//     db.err(e.get_errno(), "ShapeletObjectDB: Error in get()");
//     throw e;
//   } catch(std::exception &e) {
//     throw e;
//   }
// }

// void ShapeletObjectDB::put(unsigned int i, const ShapeletObject& sobj) {
//   try {
//     Dbt key(&i, sizeof(u_int32_t));
//     Dbt data;
//     // compute lengths of elements to store in data
//     datalengths dl;
//     setLengths(dl,sobj);
//     // allocate and fill buffer for data with values from sobj
//     u_int8_t buffer[dl.totallength];
//     fillBuffer(dl,sobj,buffer);
//     data.set_size(dl.totallength);
//     data.set_data(buffer);
//     // since keys are unique, not check is necessary for existing entries
//     int ret = cursorp->put(&key, &data, DB_KEYFIRST); 
//     db.cursor(NULL, &cursorp, 0);
//   }
//   catch(DbException &e) {
//     db.err(e.get_errno(), "ShapeletObjectDB: Error in put()");
//   } catch(std::exception &e) {
//     db.errx("ShapeletObjectDB: Error %s", e.what());
//   }
// }

// void ShapeletObjectDB::fillShapeletObject(Dbt& data) {
//   // first get header to infer sizes of variable elements in data
//   int totallength = data.get_size();
//   sdbheader header;
//   u_int8_t *p = (u_int8_t*) data.get_data();
//   memcpy(&header,p,sizeof(sdbheader));
//   // then read strings and coeffs (errors)
//   p += sizeof(sdbheader);
//   tmpsobj.basefilename = std::string((char*)p,header.basefilelength);
//   p += header.basefilelength;
//   tmpsobj.history.clear();
//   tmpsobj.history << std::string((char*)p,header.historylength);
//   p += header.historylength;
//   // determine size of coefficient vector, read it from db and convert it to matrix
//   unsigned int nCoeffs = (header.nmax+1)*(header.nmax+2)/2;
//   NumVector<double> cV(nCoeffs);
//   memcpy(cV.c_array(),p,nCoeffs*sizeof(double));
//   CoefficientVector<double> coeffVector(cV);
//   coeffVector.fillCoeffMatrix(tmpsobj.coeffs);
//   // if there is something left, its the coeff errors
//   if (sizeof(sdbheader) + header.basefilelength + header.historylength + 2*nCoeffs*sizeof(double) == totallength) {
//     p += nCoeffs*sizeof(double);
//     memcpy(cV.c_array(),p,nCoeffs*sizeof(double));
//     coeffVector = CoefficientVector<double>(cV);
//     coeffVector.fillCoeffMatrix(tmpsobj.errors);
//   }
//   std::cout << tmpsobj.errors << std::endl;
// }
// void ShapeletObjectDB::fillBuffer(const datalengths& dl, const ShapeletObject& sobj, u_int8_t* p) {

//   // set up header part of data
//   const Grid& grid = sobj.getGrid();
//   const Point2D<data_t>& centroid = sobj.getCentroid();
//   sdbheader header = {
//     1,
//     sobj.getObjectID(),
//     sobj.getNMax(),
//     sobj.getFlags().to_ulong(),
//     sobj.getBeta(),
//     sobj.getDecompositionChiSquare(),
//     sobj.getRegularizationR(),
//     sobj.getTag(),
//     sobj.getObjectClassifier(),
//     grid.getStartPosition(0),
//     grid.getStopPosition(0),
//     grid.getStartPosition(1),
//     grid.getStopPosition(1),
//     centroid(0),
//     centroid(1),
//     dl.basefilelength,
//     dl.historylength
//   };
  
//   // fill header, strings and coeffs (errors) into buffer
//   memcpy(p,&header,dl.headerlength);
//   p += dl.headerlength;
//   memcpy(p,sobj.getBaseFilename().c_str(),header.basefilelength);
//   p += header.basefilelength;
//   memcpy(p,sobj.getHistory().c_str(),header.historylength);
//   p += header.historylength;
//   CoefficientVector<double> coeffVector(sobj.getCoeffs());
//   memcpy(p,coeffVector.c_array(),dl.coefflength);
//   if (dl.errors) {
//      coeffVector = CoefficientVector<double>(sobj.getDecompositionErrors());
//      p += dl.coefflength;
//      memcpy(p,coeffVector.c_array(),dl.coefflength);
//   }
// }

// void ShapeletObjectDB::setLengths(datalengths& dl, const ShapeletObject& sobj) {
//   unsigned int nmax = sobj.getNMax();
//   dl.coefflength = (nmax+1)*(nmax+2)/2 * sizeof(double);
//   dl.errors = (bool) (sobj.getDecompositionErrors()).getRows();
//   dl.headerlength = sizeof(sdbheader);
//   dl.basefilelength = sobj.getBaseFilename().size() * sizeof(char);
//   dl.historylength = sobj.getHistory().size() * sizeof(char);
//   dl.totallength = dl.headerlength + dl.basefilelength + dl.historylength + dl.coefflength;
//   if (dl.errors)
//     dl.totallength += dl.coefflength;
// }

// #endif // Berkeley

#endif // SHAPELETDB


