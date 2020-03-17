#include <fstream>
#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

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

void load_db_index(std::string filename) {
  std::ifstream index_file(filename);

  for (std::string line; getline(index_file, line);) {
    std::cout << line << std::endl;

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
}