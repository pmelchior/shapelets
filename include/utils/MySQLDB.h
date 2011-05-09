#ifndef SHAPELENS_MYSQLDB_H
#define SHAPELENS_MYSQLDB_H

#ifdef HAS_MySQLDB

#include <string>
#include <mysql/mysql.h>

namespace shapelens { 
 
/// DB accessor class for MySQL.
 class MySQLDB {
 public:
  /// Constructor.
  MySQLDB();
  /// Destructor.
  ~MySQLDB();
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
  /// Execute a \p query using \p callback.
  void exec(std::string query, DBCallback& dbc);
  /// Retrieve data from \p query.
  MySQLResult query(std::string query);
  /// Type of DBHandle (returns 1).
  int getDBType();
 private:
  MYSQL* conn;
  std::string database;
  void connect(std::string host, std::string user, std::string password, std::string schema = "");
  
  class MySQLResult {
  public:
    MySQLResult(MYSQL* conn);
    ~MySQLResult();
    unsigned int getRowCount();
    unsigned int getFieldCount();
    char* getFieldName(unsigned int i);
    char** getRow();
    void free();
    friend class MySQLDB;
  private:
    MYSQL_RES *res;
    MYSQL_FIELD* fields;
  };
  };

} // end namespace
#endif // HAS_MYSQLDB
#endif // MYSQLDB_H
