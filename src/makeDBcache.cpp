#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <pgsql/libpq-fe.h>

struct db_entry {
  std::string schema_name;
  std::string table_name;
  std::string owner;
  std::string host;
  int port;

  db_entry(std::string _schema_name, std::string _table_name,
           std::string _owner, std::string _host, int _port) {
    schema_name = _schema_name;
    table_name = _table_name;
    owner = _owner;
    host = _host;
    port = _port;
  }
};

std::unordered_map<int, std::shared_ptr<struct db_entry>> db_index;
void make_gaia_db(int hpx, std::shared_ptr<struct db_entry> entry);

void load_db_index(std::string filename) {
  std::ifstream index_file(filename);

  for (std::string line; getline(index_file, line);) {
    // std::cout << line << std::endl;

    std::vector<std::string> columns;
    boost::split(columns, line, [](char c) { return c == '|'; });

    // extract the database connection info based on the healpix index
    if (columns.size() == 6) {
      int hpx;

      sscanf(columns[1].c_str(), "gaia_source_%d", &hpx);

      std::shared_ptr<struct db_entry> entry(
          new db_entry(columns[0], columns[1], columns[3], columns[4],
                       std::stoi(columns[5])));

      db_index.insert(std::make_pair(hpx, entry));
    }
  }

  std::cout << "PostgreSQL HEALPix index contains " << db_index.size()
            << " entries." << std::endl;
}

int main(int argc, char *argv[]) {
  // load the db healpix index file
  load_db_index("gaiadr2-table.dat");

  // go through the GAIA db once and make a cache
  for (auto &it : db_index) {
    int hpx = it.first;
    auto entry = it.second;

    make_gaia_db(hpx, entry);

    break;
  }
}

void make_gaia_db(int hpx, std::shared_ptr<struct db_entry> entry) {
  std::cout << "schema: " << entry->schema_name
            << " table: " << entry->table_name << " owner: " << entry->owner
            << " host:port: " << entry->host << ":" << entry->port << std::endl;

  std::string conn_str = "dbname=gaiadr2 host=" + entry->host +
                         " port=" + std::to_string(entry->port) +
                         " user=" + entry->owner + " password=jvo!";

  PGconn *gaia_db = PQconnectdb(conn_str.c_str());
  uint64_t count = 0;

  if (PQstatus(gaia_db) != CONNECTION_OK) {
    fprintf(stderr, "PostgreSQL connection failed: %s\n",
            PQerrorMessage(gaia_db));

    PQfinish(gaia_db);
    gaia_db = NULL;
    return;
  } else
    printf("PostgreSQL connection successful.\n");

  // clean up the db connection
  if (gaia_db != NULL)
    PQfinish(gaia_db);
}