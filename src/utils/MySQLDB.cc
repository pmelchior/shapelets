#ifdef SHAPELETDB
#if SHAPELETDB==MySQL

#include "../../include/utils/MySQLDB.h"
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <fenv.h>
#include <stdexcept>

using namespace shapelens;

DBResult::DBResult(MYSQL* conn) {
  res = mysql_store_result(conn);
}

MYSQL_ROW DBResult::getRow() {
  return mysql_fetch_row(res);
}

unsigned int DBResult::getRowCount() {
  return mysql_num_rows(res);
}

unsigned int DBResult::getFieldCount() {
  return mysql_num_fields(res);
}

MYSQL_FIELD* DBResult::getFields() {
  return mysql_fetch_fields(res);
}

DBResult::~DBResult() {
  mysql_free_result(res); 
}

MySQLDB::MySQLDB() {
  conn = mysql_init(NULL);
}

void MySQLDB::connect(std::string envvar) {
  // open file
  std::ifstream connfile (getenv(envvar.c_str()));
  if (connfile.fail()) {
    std::cerr << "MySQLDB: connection file could not be opened!" << std::endl;
    std::terminate();
  }
  // read in config file
  std::string line;
  std::string host, user, password;
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
    }
  }

  if (host.size() == 0)
    throw std::invalid_argument("MySQLDB: HOST keyword not specified!");
  if (user.size() == 0)
    throw std::invalid_argument("MySQLDB: USER keyword not specified!");

  connect(host,user,password,database);
}

void MySQLDB::connect(std::string host, std::string user, std::string password, std::string schema) {
  // try to connect to DB server
  // if schema is set, use it; otherwise use NULL
  if (schema.size() != 0) {
    if (!mysql_real_connect(conn, host.c_str(), user.c_str(), password.c_str(), schema.c_str(), 0, NULL, 0))
      throw std::invalid_argument(mysql_error(conn));
  } else {
    if (!mysql_real_connect(conn, host.c_str(), user.c_str(), password.c_str(), NULL, 0, NULL, 0))
      throw std::invalid_argument(mysql_error(conn));
  }
}

const std::string& MySQLDB::getDatabaseName() {
  return database;
}

void MySQLDB::selectDatabase(std::string d) {
  database = d;
  if (mysql_select_db(conn,database.c_str()))
    throw std::invalid_argument(mysql_error(conn));
}

MySQLDB::~MySQLDB() {
  mysql_close(conn);
}

DBResult MySQLDB::query(std::string query) {
  if (mysql_real_query(conn, query.c_str(),query.size()))
    throw std::invalid_argument(mysql_error(conn));
  return DBResult(conn);
}

#endif // SHAPELETDB==MySQL
#endif // SHAPELETDB
