#ifndef SHAPELENS_MYSQLDB_H
#define SHAPELENS_MYSQLDB_H

#if SHAPELETDB==MySQL

#include <string>
#include <mysql/mysql.h>

namespace shapelens {  
/// MySQLDB result class.
/// This class creates and frees the memory associated with a
/// \p MYSQL_RES instance.
/// Single rows of the result set are obtained from getRow();
class DBResult{
 public:
  /// Constructor from DB connection \p conn.
  DBResult(MYSQL* conn);
  /// Get number of rows in DBResult.
  unsigned int getRowCount();
  /// Get number of fields in the result.
  unsigned int getFieldCount();
  /// Return field informations for result.
  MYSQL_FIELD* getFields();
  /// Get a row from the result set.
  /// Returns \p NULL if there are no (more) rows.
  MYSQL_ROW getRow();
  /// Default destructor.
  ~DBResult();
 private:
  MYSQL_RES *res;
};

/// DB accessor class for MySQL.
class MySQLDB {
 public:
  /// Constructor.
  MySQLDB();
  /// Destructor.
  ~MySQLDB();
  /// Connect to DB with given details.
  void connect(std::string host, std::string user, std::string password, std::string schema = "");
  /// Connect to DB using a connection file.
  /// The file is specified by its path which has to be set in the 
  /// environment variable \p envvar.
  /// The file contains the following connection details
  /// and the appropriate values in format of one keyword/value pair per line,
  /// separated by one or more tab characters:
  /// \code
  /// HOST      localhost
  /// USER      username
  /// PASSWORD  pass
  /// DATABASE  shapelets
  /// \endcode
  ///
  /// \p HOST and \p USER are mandatory.\n\n
  /// \b CAUTION: If the password is written in this file, make sure to change 
  /// the file permissions such that \p username is the only user who can read 
  /// this file. 
  void connect(std::string envvar);
  /// Get name of active database.
  const std::string& getDatabaseName();
  /// Change the selected database.
  void selectDatabase(std::string database);
  /// Send query to DB.
  DBResult query(std::string query);
 private:
  MYSQL* conn;
  std::string database;
};
} // end namespace
#endif // SHAPELETDB
#endif // MYSQLDB_H
