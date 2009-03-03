#include <utils/MySQLDB.h>
#include <iostream>

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

void MySQLDB::connect(std::string host, std::string user, std::string password, std::string schema) {
  // try to connect to DB server
  // if schema is set, use it; otherwise use NULL
  if (schema.size() != 0) {
    if (!mysql_real_connect(conn, host.c_str(), user.c_str(), password.c_str(), schema.c_str(), 0, NULL, 0)) {
      std::cerr << mysql_error(conn) << std::endl;
      std::terminate();
    }
  } else {
    if (!mysql_real_connect(conn, host.c_str(), user.c_str(), password.c_str(), NULL, 0, NULL, 0)) {
      std::cerr << mysql_error(conn) << std::endl;
      std::terminate();
    }
  }
}

MySQLDB::~MySQLDB() {
  mysql_close(conn);
}

DBResult MySQLDB::query(std::string query) {
  if (mysql_real_query(conn, query.c_str(),query.size())) {
    std::cerr <<  mysql_error(conn) << std::endl;
    std::terminate();
  }
  return DBResult(conn);
}
