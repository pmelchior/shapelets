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
  /// Send query to DB.
  DBResult query(std::string query);
 private:
  MYSQL* conn;
};
} // end namespace
#endif // SHAPELETDB
#endif // MYSQLDB_H
