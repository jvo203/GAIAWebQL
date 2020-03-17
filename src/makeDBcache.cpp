#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <pgsql/libpq-fe.h>

#include "PJMCoords.h"

#define HTMfunc_Pr 3.1415926535897932385E0 / 180.0

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
  OmniCoords coords;
  coords.change_sol_pos(8.3, 0.027);

  std::cout << "schema.table: " << entry->schema_name
            << "." << entry->table_name << " owner: " << entry->owner
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

  // alter table add column
  std::string sql = "alter table " + entry->schema_name + "." + entry->table_name +
  " add column _x real,"
  " add column _y real,"
  " add column _z real,"
  " add column _r real,"
  " add column _phi real,"
  " add column _m_g real;";

  PGresult* res = PQexec(gaia_db, sql.c_str());
  ExecStatusType stat = PQresultStatus(res);

  if (stat != PGRES_COMMAND_OK)
    std::cout << sql << " : " << PQresStatus(stat) << std::endl;  

  if(res != NULL)
    PQclear(res);

  // create search indices
  std::string columns[6] = {"_x", "_y", "_z", "_r", "_phi", "_m_g"};

  for(int i=0;i<6;i++) {
    std::string column = columns[i];

    //sql = "create index on " + entry->schema_name + "." + entry->table_name + " (" + column + ");";
    sql = "create index " + entry->table_name + "_" + column + "_idx on " + entry->schema_name + "." + entry->table_name + " (" + column + ");";    

    res = PQexec(gaia_db, sql.c_str());
    ExecStatusType stat = PQresultStatus(res);

    if (stat != PGRES_COMMAND_OK)
      std::cout << sql << " : " << PQresStatus(stat) << std::endl;  

    if(res != NULL)
      PQclear(res);
  }

  // iterate through all the rows
  sql = "select "
                      "ra,dec,phot_g_mean_mag,bp_rp,parallax,pmra,pmdec,radial_"
                      "velocity,source_id from " +
                      entry->schema_name + "." + entry->table_name +
                      " where parallax > 0 and ra is not null and dec is not "
                      "null and phot_g_mean_mag is not null and bp_rp is not "
                      "null and pmra is not null and pmdec is not null and "
                      "radial_velocity is not null;";

  if (PQsendQuery(gaia_db, sql.c_str())) {
      if (PQsetSingleRowMode(gaia_db)) {
        PGresult *res = NULL;

        while ((res = PQgetResult(gaia_db)) != NULL) {
          if (PQresultStatus(res) == PGRES_SINGLE_TUPLE) {
            count++;

            std::stringstream res_str;
            res_str << count << ":\t";

            int nRows = PQntuples(res);
            int nFields = PQnfields(res);

            /*for (int i = 0; i < nFields; i++)
                res_str << PQfname(res, i) << "\t";
            res_str << std::endl;*/

            for (int i = 0; i < nRows; i++) {
              if (nFields >= 8) {
                /*char *e;
                errno = 0;*/

                std::size_t pos;
                long long source_id;
                double ra, dec, phot_g_mean_mag, bp_rp, parallax, pmra, pmdec,
                    radial_velocity;

                bool valid_data = true;

                try {
                  ra = std::stod(std::string(PQgetvalue(res, i, 0)), &pos);

                  dec = std::stod(std::string(PQgetvalue(res, i, 1)), &pos);

                  phot_g_mean_mag =
                      std::stod(std::string(PQgetvalue(res, i, 2)), &pos);

                  bp_rp = std::stod(std::string(PQgetvalue(res, i, 3)), &pos);

                  parallax =
                      std::stod(std::string(PQgetvalue(res, i, 4)), &pos);

                  pmra = std::stod(std::string(PQgetvalue(res, i, 5)), &pos);

                  pmdec = std::stod(std::string(PQgetvalue(res, i, 6)), &pos);

                  radial_velocity =
                      std::stod(std::string(PQgetvalue(res, i, 7)), &pos);

                  source_id = std::stoll(std::string(PQgetvalue(res, i, 8)), &pos);
                } catch (const std::out_of_range &err) {
                  printf("(%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s)\n",
                         PQgetvalue(res, i, 0), PQgetvalue(res, i, 1),
                         PQgetvalue(res, i, 2), PQgetvalue(res, i, 3),
                         PQgetvalue(res, i, 4), PQgetvalue(res, i, 5),
                         PQgetvalue(res, i, 6), PQgetvalue(res, i, 7),
                         PQgetvalue(res, i, 8));
                  valid_data = false;
                } catch (const std::invalid_argument &err) {
                  printf("(%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s)\n",
                         PQgetvalue(res, i, 0), PQgetvalue(res, i, 1),
                         PQgetvalue(res, i, 2), PQgetvalue(res, i, 3),
                         PQgetvalue(res, i, 4), PQgetvalue(res, i, 5),
                         PQgetvalue(res, i, 6), PQgetvalue(res, i, 7),
                         PQgetvalue(res, i, 7));
                  valid_data = false;
                }

                if (valid_data) {
                  std::cout << source_id << "¥t";
                  double M_G =
                      phot_g_mean_mag + 5.0 + 5.0 * log10(parallax / 1000.0);

                  double alpha = ra * HTMfunc_Pr;  //[rad]
                  double delta = dec * HTMfunc_Pr; //[rad]
                  double d = 1000.0 / parallax;    // distance [parsec, pc]

                  vec6 sHEQ, sGCA, sGCY;
                  sHEQ[0] = d / 1000.0;      //[kpc]
                  sHEQ[1] = ra;              //[deg]
                  sHEQ[2] = dec;             //[deg]
                  sHEQ[3] = radial_velocity; //[km/s]
                  sHEQ[4] = pmra;
                  sHEQ[5] = pmdec;

                  sGCA = coords.GCAfromHEQ(sHEQ);
                  sGCY = coords.GCYfromHEQ(sHEQ);

                  double X = sGCA[0]; //[kpc]
                  double Y = sGCA[1]; //[kpc]
                  double Z = sGCA[2]; //[kpc]

                  double R = sGCY[0];    //[kpc]
                  double Phi = sGCY[2];  //[rad]
                  double VR = sGCY[3];   //[km/s]
                  double VZ = sGCY[4];   //[km/s]
                  double VPhi = sGCY[5]; //[km/s]

                  // update the db row

                }
              }
            }
          };

          PQclear(res);
        };
      } else
        std::cout << "error setting PQsetSingleRowMode.\n";
    } else
      std::cout << "PQsendQuery error.\n";

  // clean up the db connection
  if (gaia_db != NULL)
    PQfinish(gaia_db);
}