#ifdef SHAPELETDB
#ifndef SHAPELENS_SHAPELETOBJECTDB_H
#define SHAPELENS_SHAPELETOBJECTDB_H

#include <string>
#include "../Typedef.h"
#include "ShapeletObject.h"
#include "ShapeletObjectList.h"

#if SHAPELETDB==MySQL
#include "../utils/MySQLDB.h"
#endif

namespace shapelens {

/// Database backend class.
/// This class provides an interface to a database server/file (later called \p DB).
/// Single ShapeletObject entities or whole ShapeletObjectList entities can be stored
/// in and retrieved from a \p DB.\n\n
/// Depending on the permissions granted to the \p DB user who opens the connection,
/// some of the calls above may be forbidden.\n\n
/// The database table layout is as follows:
/// \code
/// CREATE TABLE  `DATABASE`.`TABLE` (
///  `id` int(10) unsigned NOT NULL default '0' COMMENT 'object id',
///  `nmax` tinyint(3) unsigned NOT NULL default '0' COMMENT 'shapelet order',
///  `beta` float unsigned NOT NULL default '0' COMMENT 'shapelet scale',
///  `chi2` float unsigned NOT NULL default '0' COMMENT 'goodness of fit',
///  `flags` smallint(5) unsigned NOT NULL default '0' COMMENT 'segmentation and decomposition flags',
///  `min_x` int(10) unsigned NOT NULL default '0' COMMENT 'min(X) in image pixels',
///  `min_y` int(10) unsigned NOT NULL default '0' COMMENT 'min(Y) in image pixels',
///  `size_x` smallint(5) unsigned NOT NULL default '0' COMMENT 'image dimension in x',
///  `size_y` smallint(5) unsigned NOT NULL default '0' COMMENT 'image dimension in y',
///  `centroid_x` float unsigned NOT NULL default '0' COMMENT 'x-position of centroid in image coordinates',
///  `centroid_y` float unsigned NOT NULL default '0' COMMENT 'y-Position of centroid in image coordinates',
///  `basefile` varchar(255) default NULL COMMENT 'object source file',
///  `prop` text COMMENT 'arbitrary object properties',
///  `history` text COMMENT 'object history',
///  `coeffs` blob NOT NULL COMMENT 'shapelet coefficients',
///  `cov` mediumblob COMMENT 'covariance matrix of shapelet coefficients',
///  PRIMARY KEY  (`id`)
///) ENGINE=MyISAM DEFAULT CHARSET=latin1
/// \endcode
///
/// This means, that essentially all properties of a ShapeletObject entity are
/// stored in the \p DB. As a peculiarity, the column \p id, which
/// is populated with ShapeletObject::getObjectID(), has to have unique entries.
///
/// Currently, \p MySQL is the only supported \p DB driver, but extensions to
/// \p BerkeleyDB or \p SQLite are straightforward.
/// 
/// \todo
/// - allow queries without history and/or covariance matrix to reduce bandwith
/// - implement additional 'book-keeping' table to store ShapeLensConfig parameters,
///   user, creating time etc.

class ShapeletObjectDB {
 public:
  /// Default constructor.
  /// The required connection details are stored in a file whose
  /// path has to be set in the environment variable \p SHAPELENSDBCONF.\n
  /// see MySQLDB::connect() for details.
  
  ShapeletObjectDB();
  /// Change the table within the database.
  void selectTable(std::string table);
  /// Change the selected database.
  void selectDatabase(std::string database);
  /// Load all ShapeletObject entities which match the \p WHERE statement.
  /// \p where_clause is the string which specifies the selection within the
  /// \p DB table. E.g.
  /// \code
  /// ShapeletObjectDB db(connection_file);
  /// ShapeletObjectList sl = db.load("`chi2` < '1' AND `nmax` >= '8'");
  /// \endcode
  /// would open a \p DB connection and load all objects with \f$\chi^2 < 1\f$ and
  /// \f$n_{max}\geq 8\f$.
  ShapeletObjectList load(std::string where_clause = "", std::string join_clause = "");
  /// Save all ShapeletObject entities in \p sl to the current \p table.
  /// This call requires the permission to execute a \p INSERT statement.
  void save(const ShapeletObjectList& sl);
  /// Save \p sobj to the current \p table.
  /// This call requires the permission to execute a \p INSERT statement.
  void save(const ShapeletObject& sobj);
  /// Clear all entries from the current \p table.
  /// This call requires the permission to execute a \p TRUNCATE statement.\n\n
  /// \b CAUTION: This call is irreversible.
  void clear();
  /// Submit query to \p DB.
  /// The call requires permissions to execute the specified statement.
  DBResult query(std::string query);
  /// Load history for ShapeletObject during load().
  /// Using <tt>use==false</tt> lowers DB traffic and memory consumption.
  void useHistory(bool use);
  /// Load coefficient covariance matrix for ShapeletObject during load().
  /// Using <tt>use==false</tt> lowers DB traffic and memory consumption.
  void useCovariance(bool use);

 private:
  std::string table;
  void createTable();
  void checkForTable();
  std::string ref(std::string col);
  bool exists, loadHistory, loadCovariance;
#if SHAPELETDB==MySQL
  MySQLDB db;
#endif
};
} // end namespace
#endif  // SHAPELETOBJECTDB_H
#endif  // SHAPELETDB
