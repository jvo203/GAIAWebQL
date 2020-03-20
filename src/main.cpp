#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_SUB 0

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

static const double halfpi = 1.570796326794896619231321691639751442099;
#define HTMfunc_Pr 3.1415926535897932385E0 / 180.0

/*#include <math.h>

const double R[3][3] = {{-0.0548755604162154, -0.8734370902348850,
-0.4838350155487132}, {+0.4941094278755837, -0.4448296299600112,
+0.7469822444972189}, {-0.8676661490190047, -0.1980763734312015,
+0.4559837761750669}}; const double _theta = 0.003253017785496385; const double
H[3][3] = {{cos(_theta), 0.0, sin(_theta)}, {0.0, 1.0, 0.0}, {-sin(_theta), 0.0,
cos(_theta)}}; const double dGC = 8300.0;*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define SERVER_PORT 8081
#define SERVER_STRING                                                          \
  "GAIAWebQL v" STR(VERSION_MAJOR) "." STR(VERSION_MINOR) "." STR(VERSION_SUB)
#define VERSION_STRING "SV2020-03-20.0"

#include <pwd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// base64 encoding with SSL
#include <openssl/bio.h>
#include <openssl/buffer.h>
#include <openssl/evp.h>
#include <openssl/hmac.h>
#include <openssl/sha.h>

// CERN ROOT
#include <TROOT.h>

#include "PJMCoords.h"

char *base64(const unsigned char *input, int length);

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_map>
using std::chrono::steady_clock;

#include <curl/curl.h>
#include <pgsql/libpq-fe.h>
//#include "healpix_base.h"
#include "json.h"

#include "SeedHist2D.hpp"
#include "SeedHist3D.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/atomic.hpp>
#include <boost/lockfree/stack.hpp>

struct search_criteria {
  double X_min, X_max;
  double Y_min, Y_max;
  double Z_min, Z_max;
  double R_min, R_max;
  double Phi_min, Phi_max;
  double M_G_min, M_G_max;
  double parallax_error_min;

  search_criteria() {
    X_min = NAN;
    X_max = NAN;
    Y_min = NAN;
    Y_max = NAN;
    Z_min = NAN;
    Z_max = NAN;
    R_min = NAN;
    R_max = NAN;
    Phi_min = NAN;
    Phi_max = NAN;
    M_G_min = NAN;
    M_G_max = NAN;
    parallax_error_min = NAN;
  }
};

std::string make_uuid(std::shared_ptr<struct search_criteria> search,
                      std::string where) {
  std::ostringstream uri;

  if (!std::isnan(search->X_min))
    uri << search->X_min;

  if (!std::isnan(search->X_max))
    uri << search->X_max;

  if (!std::isnan(search->Y_min))
    uri << search->Y_min;

  if (!std::isnan(search->Y_max))
    uri << search->Y_max;

  if (!std::isnan(search->Z_min))
    uri << search->Z_min;

  if (!std::isnan(search->Z_max))
    uri << search->Z_max;

  if (!std::isnan(search->R_min))
    uri << search->R_min;

  if (!std::isnan(search->R_max))
    uri << search->R_max;

  if (!std::isnan(search->Phi_min))
    uri << search->Phi_min;

  if (!std::isnan(search->Phi_max))
    uri << search->Phi_max;

  if (!std::isnan(search->M_G_min))
    uri << search->M_G_min;

  if (!std::isnan(search->M_G_max))
    uri << search->M_G_max;

  if (!std::isnan(search->parallax_error_min))
    uri << search->parallax_error_min;

  if (where != "")
    uri << where;

  std::size_t id = std::hash<std::string>{}(uri.str());

  std::stringstream sstream;
  sstream << std::hex << id;
  std::string uuid = sstream.str();

  return uuid;
}

void print_search_criteria(std::shared_ptr<struct search_criteria> search) {
  printf("X: [%f - %f]\n", search->X_min, search->X_max);
  printf("Y: [%f - %f]\n", search->Y_min, search->Y_max);
  printf("Z: [%f - %f]\n", search->Z_min, search->Z_max);
  printf("R: [%f - %f]\n", search->R_min, search->R_max);
  printf("Φ: [%f - %f]\n", search->Phi_min, search->Phi_max);
  printf("M_G: [%f - %f]\n", search->M_G_min, search->M_G_max);
  printf("parallax_error: ≥%f\n", search->parallax_error_min);
}

struct plot_data {
  double bp_rp;
  double M_G;
  double X;
  double Y;
  double Z;
  double R;
  double VR;
  double VPhi;
  double VZ;
};

struct global_search {
  std::vector<struct plot_data> queue;
  std::mutex queue_mtx;
};

#include <omp.h>

#include <nghttp2/asio_http2_server.h>

using namespace nghttp2::asio_http2;
using namespace nghttp2::asio_http2::server;

http2 *http2_server;
std::string docs_root = "htdocs";
std::string home_dir;

void http_ok(const response *res) {
  res->write_head(200);
  res->end("OK");
}

void http_created(const response *res) {
  res->write_head(201);
  res->end("Created");
}

void http_accepted(const response *res) {
  res->write_head(202);
  res->end("Accepted");
}

void http_not_found(const response *res) {
  res->write_head(404);
  res->end("Not Found");
}

void http_no_content(const response *res) {
  res->write_head(204);
  res->end("No Content");
}

void http_not_implemented(const response *res) {
  res->write_head(501);
  res->end("Not Implemented");
}

void http_service_unavailable(const response *res) {
  res->write_head(503);
  res->end("Service Unavailable");
}

void http_internal_server_error(const response *res) {
  res->write_head(500);
  res->end("Internal Server Error");
}

void serve_file(const request *req, const response *res, std::string uri) {
  // a safety check against directory traversal attacks
  if (!check_path(uri))
    return http_not_found(res);

  // check if a resource exists
  std::string path = docs_root + uri;

  if (std::filesystem::exists(path)) {
    // detect mime-types
    header_map mime;

    size_t pos = uri.find_last_of(".");

    if (pos != std::string::npos) {
      std::string ext = uri.substr(pos + 1, std::string::npos);

      if (ext == "htm" || ext == "html")
        mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                         {"text/html", false}));

      if (ext == "txt")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"text/plain", false}));

      if (ext == "js")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"application/javascript", false}));

      if (ext == "ico")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"image/x-icon", false}));

      if (ext == "png")
        mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                         {"image/png", false}));

      if (ext == "gif")
        mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                         {"image/gif", false}));

      if (ext == "webp")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"image/webp", false}));

      if (ext == "jpg" || ext == "jpeg")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"image/jpeg", false}));

      if (ext == "mp4")
        mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                         {"video/mp4", false}));

      if (ext == "css")
        mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                         {"text/css", false}));

      if (ext == "pdf")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"application/pdf", false}));

      if (ext == "svg")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"image/svg+xml", false}));

      if (ext == "root")
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"application/octet-stream", false}));
    }

    // check for compression
    header_map headers = req->header();

    auto it = headers.find("accept-encoding");
    if (it != headers.end()) {
      auto value = it->second.value;
      // std::cout << "Supported compression: " << value << std::endl;

      // prefer brotli due to smaller file sizes
      size_t pos = value.find("br"); // brotli or gzip

      if (pos != std::string::npos) {
        if (std::filesystem::exists(path + ".br")) {
          path += ".br";
          // append the compression mime
          mime.insert(std::pair<std::string, header_value>("Content-Encoding",
                                                           {"br", false}));
        }
      } else {
        // fallback to gzip
        size_t pos = value.find("gzip");

        if (pos != std::string::npos) {
          if (std::filesystem::exists(path + ".gz")) {
            path += ".gz";
            // append the compression mime
            mime.insert(std::pair<std::string, header_value>("Content-Encoding",
                                                             {"gzip", false}));
          }
        }
      }
    }

    res->write_head(200, mime);
    res->end(file_generator(path));
  } else {
    std::cout << "[HTTPServer]: " << path << " not found." << std::endl;
    http_not_found(res);
  }
}

int server_port = SERVER_PORT;

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

struct db_search_job {
  // steady_clock::time_point timestamp = steady_clock::now();
  time_t timestamp = std::time(nullptr);
  std::vector<int> completed;
  std::mutex completed_mtx;
  bool abort = false;
};

std::unordered_map<std::string, std::shared_ptr<struct db_search_job>> requests;
std::mutex requests_mtx;

struct gaia_plots {
  std::string hr;
  std::string xy;
  std::string rz;
  std::string xyvr;
  std::string xyvphi;
  std::string xyvz;
  std::string rzvr;
  std::string rzvphi;
  std::string rzvz;
};

std::unordered_map<std::string, std::shared_ptr<struct gaia_plots>> results;
std::mutex results_mtx;

struct gaia_hist {
  // 2D
  SeedH2 _hr;
  SeedH2 _xy;
  SeedH2 _rz;

  // 3D
  SeedH3 _xyvr;
  SeedH3 _xyvphi;
  SeedH3 _xyvz;
  SeedH3 _rzvr;
  SeedH3 _rzvphi;
  SeedH3 _rzvz;

  gaia_hist() {
    // 3D
    _hr = SeedH2(true);
    _xy = SeedH2(false);
    _rz = SeedH2(false);

    // 3D
    _xyvr = SeedH3(false);
    _xyvphi = SeedH3(false);
    _xyvz = SeedH3(false);
    _rzvr = SeedH3(false);
    _rzvphi = SeedH3(false);
    _rzvz = SeedH3(false);
  }
};

std::mutex root_mtx;

inline const char *check_null(const char *str) {
  if (str != nullptr)
    return str;
  else
    return "\"\"";
};

void signalHandler(int signum) {
  std::cout << "Interrupt signal (" << signum << ") received.\n";

  curl_global_cleanup();

  // cleanup and close up stuff here
  // terminate program

  std::cout << "remaining requests: " << requests.size()
            << "\tremaining results: " << results.size() << std::endl;

  http2_server->stop();
}

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

bool dataset_exists(std::string uuid) {
  std::string dir = docs_root + "/gaiawebql/DATA/" + uuid;

  if (std::filesystem::exists(dir))
    return true;
  else
    return false;
}

bool search_gaia_db(int hpx, std::shared_ptr<struct db_entry> entry, std::string uuid,
                    std::shared_ptr<struct search_criteria> search, std::string where,
                    std::shared_ptr<OmniCoords> coords,
                    boost::lockfree::stack<struct plot_data> &plot_stack,
                    struct global_search& queue
                    /*,
                    struct gaia_hist &hist*/) {
  std::stringstream msg;
  bool abort_search = false;

  msg << uuid << ":\t" << entry->schema_name << "/" << entry->table_name << "/"
      << entry->owner << "/" << entry->host << ":" << entry->port << "\t";

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
  } else {
    printf("PostgreSQL connection successful.\n");

    // std::vector<struct plot_data> local_queue;

    // perform a whole-data search
    // and ra is not null and dec is not null and phot_g_mean_mag is not null
    // and bp_rp is not null and pmra is not null and pmdec is not null and
    // radial_velocity is not null
    std::string sql = "select "
                      "ra,dec,phot_g_mean_mag,bp_rp,parallax,pmra,pmdec,radial_"
                      "velocity from " +
                      entry->schema_name + "." + entry->table_name +
                      " where parallax > 0 and ra is not null and dec is not "
                      "null and phot_g_mean_mag is not null and bp_rp is not "
                      "null and pmra is not null and pmdec is not null and "
                      "radial_velocity is not null"; // and parallax_over_error
                                                     // > 10;"; //" limit 1;";

    // add optional search conditions
    if (where != "")
      sql.append(" and " + where);

    // parallax error
    if (!std::isnan(search->parallax_error_min))
      sql.append(" and parallax_over_error >= " +
                 std::to_string(search->parallax_error_min));

    // validate data against the search critera
    if (!std::isnan(search->X_min))
      sql.append(" and _x >= " + std::to_string(search->X_min));

    if (!std::isnan(search->X_max))
      sql.append(" and _x <= " + std::to_string(search->X_max));

    if (!std::isnan(search->Y_min))
      sql.append(" and _y >= " + std::to_string(search->Y_min));

    if (!std::isnan(search->Y_max))
      sql.append(" and _y <= " + std::to_string(search->Y_max));

    if (!std::isnan(search->Z_min))
      sql.append(" and _z >= " + std::to_string(search->Z_min));

    if (!std::isnan(search->Z_max))
      sql.append(" and _z <= " + std::to_string(search->Z_max));

    if (!std::isnan(search->R_min))
      sql.append(" and _r >= " + std::to_string(search->R_min));

    if (!std::isnan(search->R_max))
      sql.append(" and _r <= " + std::to_string(search->R_max));

    if (!std::isnan(search->Phi_min))
      sql.append(" and _phi >= " +
                 std::to_string(search->Phi_min * HTMfunc_Pr));

    if (!std::isnan(search->Phi_max))
      sql.append(" and _phi <= " +
                 std::to_string(search->Phi_max * HTMfunc_Pr));

    if (!std::isnan(search->M_G_min))
      sql.append(" and _m_g >= " + std::to_string(search->M_G_min));

    if (!std::isnan(search->M_G_max))
      sql.append(" and _m_g <= " + std::to_string(search->M_G_max));

    sql.append(";"); // finish the sql

    // std::cout << sql << std::endl;

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
                } catch (const std::out_of_range &err) {
                  printf("(%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s)\n",
                         PQgetvalue(res, i, 0), PQgetvalue(res, i, 1),
                         PQgetvalue(res, i, 2), PQgetvalue(res, i, 3),
                         PQgetvalue(res, i, 4), PQgetvalue(res, i, 5),
                         PQgetvalue(res, i, 6), PQgetvalue(res, i, 7));
                  valid_data = false;
                } catch (const std::invalid_argument &err) {
                  printf("(%s) (%s) (%s) (%s) (%s) (%s) (%s) (%s)\n",
                         PQgetvalue(res, i, 0), PQgetvalue(res, i, 1),
                         PQgetvalue(res, i, 2), PQgetvalue(res, i, 3),
                         PQgetvalue(res, i, 4), PQgetvalue(res, i, 5),
                         PQgetvalue(res, i, 6), PQgetvalue(res, i, 7));
                  valid_data = false;
                }

                if (valid_data) {
                  double M_G =
                      phot_g_mean_mag + 5.0 + 5.0 * log10(parallax / 1000.0);

                  // res_str << "ra: " << ra << "\t dec: " << dec <<
                  // "phot_g_mean_mag = " << phot_g_mean_mag << "\tbp_rp = " <<
                  // bp_rp << "\tparallax = " << parallax << "\tM_G = " << M_G
                  // << std::endl;

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

                  sGCA = coords->GCAfromHEQ(sHEQ);
                  sGCY = coords->GCYfromHEQ(sHEQ);

                  double X = sGCA[0]; //[kpc]
                  double Y = sGCA[1]; //[kpc]
                  double Z = sGCA[2]; //[kpc]

                  double R = sGCY[0];    //[kpc]
                  double Phi = sGCY[2];  //[rad]
                  double VR = sGCY[3];   //[km/s]
                  double VZ = sGCY[4];   //[km/s]
                  double VPhi = sGCY[5]; //[km/s]

                  bool data_ok = true;

                  // validate data against the search critera
                  /*if (!std::isnan(search->X_min))
                    if (X < search->X_min)
                      data_ok = false;

                  if (!std::isnan(search->X_max))
                    if (X > search->X_max)
                      data_ok = false;

                  if (!std::isnan(search->Y_min))
                    if (Y < search->Y_min)
                      data_ok = false;

                  if (!std::isnan(search->Y_max))
                    if (Y > search->Y_max)
                      data_ok = false;

                  if (!std::isnan(search->Z_min))
                    if (Z < search->Z_min)
                      data_ok = false;

                  if (!std::isnan(search->Z_max))
                    if (Z > search->Z_max)
                      data_ok = false;

                  if (!std::isnan(search->R_min))
                    if (R < search->R_min)
                      data_ok = false;

                  if (!std::isnan(search->R_max))
                    if (R > search->R_max)
                      data_ok = false;

                  if (!std::isnan(search->Phi_min))
                    if (Phi < search->Phi_min * HTMfunc_Pr)
                      data_ok = false;

                  if (!std::isnan(search->Phi_max))
                    if (Phi > search->Phi_max * HTMfunc_Pr)
                      data_ok = false;

                  if (!std::isnan(search->M_G_min))
                    if (M_G < search->M_G_min)
                      data_ok = false;

                  if (!std::isnan(search->M_G_max))
                    if (M_G > search->M_G_max)
                      data_ok = false;*/

                  // the data is OK, add values to histograms
                  if (data_ok) {
                    struct plot_data data = {bp_rp, M_G, X,    Y, Z,
                                             R,     VR,  VPhi, VZ};

                    // local_queue.push_back(data);

                    while (!plot_stack.push(data))
                      ;

                    // printf("added struct plot_data\n");
                  }
                }
              }
            }
          };

          PQclear(res);
        };
      } else
        std::cout << "error setting PQsetSingleRowMode.\n";

      /*if (local_queue.size() > 0) {
        // append local queue data to the mutex-protected global queue
        std::lock_guard<std::mutex> lock(queue.queue_mtx);
        queue.queue.reserve(queue.queue.size() +
                            local_queue.size()); // replaced reserve with resize
        std::copy(local_queue.begin(), local_queue.end(), queue.queue.end());
        queue.queue.resize(queue.queue.size() + local_queue.size());

        std::cout << entry->schema_name << "/" << entry->table_name << "/"
                  << entry->owner << "/" << entry->host << ":" << entry->port
                  << "\t"
                  << " copied " << local_queue.size()
                  << " records; queue.queue.length(): " << queue.queue.size()
                  << std::endl;
        ;
      }*/
    } else
      std::cout << "PQsendQuery error.\n";
  }

  if (gaia_db != NULL)
    PQfinish(gaia_db);

  try {
    auto entry = requests.at(uuid);
    std::lock_guard<std::mutex> lock(entry->completed_mtx);
    entry->completed.push_back(hpx);
    abort_search = entry->abort;
  } catch (const std::out_of_range &err) {
    printf("no entry found for a job request %s\n", uuid.c_str());
    abort_search = true;
  }

  msg << "processed " << count << " records." << std::endl;
  std::cout << msg.str();

  return abort_search;
}

void execute_gaia(const response *res,
                  std::shared_ptr<struct search_criteria> search,
                  std::string where, std::string uuid, bool offline) {

  bool underway = false;

  // find out if the uuid has already been processed
  bool exists = dataset_exists(uuid);

  // check if processing a request is already under way
  {
    std::lock_guard<std::mutex> req_lock(requests_mtx);
    if (requests.find(uuid) != requests.end()) {
      underway = true;

      // respond with the dataset id
      if (offline) {
        std::string html = uuid;

        header_map mime;
        mime.insert(std::pair<std::string, header_value>(
            "Content-Type", {"text/plain", false}));
        res->write_head(200, mime);
        res->end(html);
        return;
      }

      /*std::string html = uuid + " is being processed. Please check later.";
      res->write_head(202);
      res->end(html);
      return;*/
    }
  }

  if (!exists && !underway) {
    {
      std::lock_guard<std::mutex> req_lock(requests_mtx);
      requests.insert(
          std::make_pair(uuid, std::make_shared<struct db_search_job>()));
    }

    {
      std::lock_guard<std::mutex> res_lock(results_mtx);
      results.insert(
          std::make_pair(uuid, std::make_shared<struct gaia_plots>()));
    }

    std::thread([=]() {
      printf("starting a GAIA database search thread, job id:%s\n",
             uuid.c_str());

      int max_threads = omp_get_max_threads();
      // std::vector<std::shared_ptr<struct gaia_hist>> thread_hist;
      // //(max_threads);
      std::vector<std::shared_ptr<OmniCoords>> thread_coords(max_threads);

      struct gaia_hist global_hist;
      char name[255];

      // 2D histograms
      global_hist._hr.set_title(uuid + "::HR", "Hertzsprung-Russell diagram",
                                "(BP-RP) [mag]", "M_{G} [mag]");
      global_hist._xy.set_title(uuid + "::XY", "X-Y", "X [kpc]", "Y [kpc]");
      global_hist._rz.set_title(uuid + "::RZ", "R-Z", "R [kpc]", "Z [kpc]");

      // 3D histograms
      global_hist._xyvr.set_title(uuid + "::XYVR", "X-Y-#bar{V}_{R}", "X [kpc]",
                                  "Y [kpc]", "#bar{V}_{R} [km/s]");
      global_hist._xyvr.set_z_axis("#bar{V}_{R}", "[km/s]");

      global_hist._xyvphi.set_title(uuid + "::XYVPhi", "X-Y-#bar{V}_{#Phi}",
                                    "X [kpc]", "Y [kpc]",
                                    "#bar{V}_{#Phi} [km/s]");
      global_hist._xyvphi.set_z_axis("#bar{V}_{#Phi}", "[km/s]");

      global_hist._xyvz.set_title(uuid + "::XYVZ", "X-Y-#bar{V}_{Z}", "X [kpc]",
                                  "Y [kpc]", "#bar{V}_{Z} [km/s]");
      global_hist._xyvz.set_z_axis("#bar{V}_{Z}", "[km/s]");

      global_hist._rzvr.set_title(uuid + "::RZVR", "R-Z-#bar{V}_{R}", "R [kpc]",
                                  "Z [kpc]", "#bar{V}_{R} [km/s]");
      global_hist._rzvr.set_z_axis("#bar{V}_{R}", "[km/s]");

      global_hist._rzvphi.set_title(uuid + "::RZVPhi", "R-Z-#bar{V}_{#Phi}",
                                    "R [kpc]", "Z [kpc]",
                                    "#bar{V}_{#Phi} [km/s]");
      global_hist._rzvphi.set_z_axis("#bar{V}_{#Phi}", "[km/s]");

      global_hist._rzvz.set_title(uuid + "::RZVZ", "R-Z-#bar{V}_{Z}", "R [kpc]",
                                  "Z [kpc]", "#bar{V}_{Z} [km/s]");
      global_hist._rzvz.set_z_axis("#bar{V}_{Z}", "[km/s]");

      for (int i = 0; i < max_threads; i++) {
        thread_coords[i] = std::make_shared<OmniCoords>(OmniCoords());
        thread_coords[i]->change_sol_pos(8.3, 0.027);

        char name[255];
        sprintf(name, "%s/%d", uuid.c_str(), (i + 1));

        // struct gaia_hist hist;
        // hist.hr_hist = new TH2D(name, "M_{G} vs. BP-RP", 600, -1.0, 5.0, 600,
        // -5.0, 15.0);
        /*hist->SetCanExtend(TH1::kAllAxes);
        hist->SetBit(TH1::kAutoBinPTwo);*/
        // thread_hist[i] = hist;
        // thread_hist.push_back(std::make_shared<struct
        // gaia_hist>(std::move(hist)));
      }

      struct global_search queue {};

      boost::lockfree::stack<struct plot_data> plot_stack(100000);
      boost::atomic<bool> search_done(false);

      std::thread plot_thread([&]() {
        std::cout << "starting a histogram thread function.\n";

        unsigned long counter = 0;
        struct plot_data data {};

        auto plotter = [&](struct plot_data data) {
          // 2D
          global_hist._hr.update(data.bp_rp, data.M_G);
          global_hist._xy.update(data.X, data.Y);
          global_hist._rz.update(data.R, data.Z);

          // 3D
          global_hist._xyvr.update(data.X, data.Y, data.VR);
          global_hist._xyvphi.update(data.X, data.Y, data.VPhi);
          global_hist._xyvz.update(data.X, data.Y, data.VZ);
          global_hist._rzvr.update(data.R, data.Z, data.VR);
          global_hist._rzvphi.update(data.R, data.Z, data.VPhi);
          global_hist._rzvz.update(data.R, data.Z, data.VZ);
        };

        // older Boost on py1 does not have consume_all!!!

        while (!search_done) {
          // counter += plot_stack.consume_all(plotter);
          while (plot_stack.pop(data)) {
            plotter(data);
            counter++;
          }
        }

        // counter += plot_stack.consume_all(plotter);
        while (plot_stack.pop(data)) {
          plotter(data);
          counter++;
        }

        global_hist._hr.flush();
        global_hist._xy.flush();
        global_hist._rz.flush();

        global_hist._xyvr.flush();
        global_hist._xyvphi.flush();
        global_hist._xyvz.flush();

        global_hist._rzvr.flush();
        global_hist._rzvphi.flush();
        global_hist._rzvz.flush();

        std::cout << "a histogram thread ended after processing " << counter
                  << " values.\n";
      });

      boost::atomic<bool> search_aborted(false);

#pragma omp parallel shared(global_hist)
      {
#pragma omp single
        {
          for (auto &it : db_index)
            if (!search_aborted) {
#pragma omp task
              {
                int hpx = it.first;
                auto entry = it.second;

                // for (size_t i = 0; i < pixels->size(); i++)
                int tid = omp_get_thread_num();
                bool abort_search = false;

                printf("parametric search through hpx %d\n", hpx);

                abort_search = search_gaia_db(hpx, entry, uuid, search, where,
                                              thread_coords[tid], plot_stack,
                                              queue); //, global_hist);

                if (abort_search) {
                  search_aborted = true;
                  printf("cancelling a parallel OpenMP search loop.\n");
                }
              }
            }
        }
      }

      printf("OpenMP parallel for loop done.\n");

      search_done = true;
      plot_thread.join();
      std::cout << "a global queue length: " << queue.queue.size()
                << "\tabort status: " << search_aborted << std::endl;

      // if (!search_aborted)
      {
        /*
#pragma omp parallel
        {
#pragma omp single
          {
            // save as JSON too for plotly.js to use

#pragma omp task
            global_hist._hr.save(uuid, docs_root, "hr");

#pragma omp task
            global_hist._xy.save(uuid, docs_root, "xy");

#pragma omp task
            global_hist._rz.save(uuid, docs_root, "rz");
          }
        }*/

        // save the histograms to disk
        // cannot be done in parallel due to the global ROOT lock
        // 2D
        global_hist._hr.export_root(uuid, docs_root, "hr");
        global_hist._xy.export_root(uuid, docs_root, "xy");
        global_hist._rz.export_root(uuid, docs_root, "rz");

        // 3D
        global_hist._xyvr.export_root(uuid, docs_root, "xyvr");
        global_hist._xyvphi.export_root(uuid, docs_root, "xyvphi");
        global_hist._xyvz.export_root(uuid, docs_root, "xyvz");
        global_hist._rzvr.export_root(uuid, docs_root, "rzvr");
        global_hist._rzvphi.export_root(uuid, docs_root, "rzvphi");
        global_hist._rzvz.export_root(uuid, docs_root, "rzvz");

        // rename the temporary dir to just "DATA/uuid"
        std::string tmp = docs_root + "/gaiawebql/DATA/" + uuid + ".tmp";
        std::string dir = docs_root + "/gaiawebql/DATA/" + uuid;
        rename(tmp.c_str(), dir.c_str());

        // the H-R diagram
        {
            /*ReverseYData(global_hist._hr.hist);
            global_hist._hr.hist->SetStats(false);
            global_hist._hr.hist->GetXaxis()->SetTitle("(BP-RP) [mag]");
            global_hist._hr.hist->GetYaxis()->SetTitle("M_{G} [mag]");*/
            // global_hist->GetZaxis()->SetTitle("star density");
        }

        // the X-Y plot
        {
            /*global_hist._xy.hist->SetStats(false);
            global_hist._xy.hist->GetXaxis()->SetTitle("X [kpc]");
            global_hist._xy.hist->GetYaxis()->SetTitle("Y [kpc]");
            SetView2D(global_hist._xy.hist);*/
        }

        // the R-Z plot
        {
            /*global_hist._rz.hist->SetStats(false);
            global_hist._rz.hist->GetXaxis()->SetTitle("R [kpc]");
            global_hist._rz.hist->GetYaxis()->SetTitle("Z [kpc]");
            SetView2D(global_hist._rz.hist);*/
        }

        // the X-Y-VR plot
        {
            /*global_hist.XYVR->SetStats(false);
            global_hist.XYVR->GetXaxis()->SetTitle("X [kpc]");
            global_hist.XYVR->GetYaxis()->SetTitle("Y [kpc]");
            global_hist.XYVR->GetZaxis()->SetTitle("V_{R} [km/s]");
            SetView3D(global_hist.XYVR);*/
        }

        // the X-Y-VPhi plot
        {
            /*global_hist.XYVPhi->SetStats(false);
            global_hist.XYVPhi->GetXaxis()->SetTitle("X [kpc]");
            global_hist.XYVPhi->GetYaxis()->SetTitle("Y [kpc]");
            global_hist.XYVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
            SetView3D(global_hist.XYVPhi);*/
        }

        // the X-Y-VZ plot
        {
            /*global_hist.XYVZ->SetStats(false);
            global_hist.XYVZ->GetXaxis()->SetTitle("X [kpc]");
            global_hist.XYVZ->GetYaxis()->SetTitle("Y [kpc]");
            global_hist.XYVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
            SetView3D(global_hist.XYVZ);*/
        }

        // the R-Z-VR plot
        {
            /*global_hist.RZVR->SetStats(false);
            global_hist.RZVR->GetXaxis()->SetTitle("R [kpc]");
            global_hist.RZVR->GetYaxis()->SetTitle("Z [kpc]");
            global_hist.RZVR->GetZaxis()->SetTitle("V_{R} [km/s]");
            SetView3D(global_hist.RZVR);*/
        }

        // the R-Z-VPhi plot
        {
            /*global_hist.RZVPhi->SetStats(false);
            global_hist.RZVPhi->GetXaxis()->SetTitle("R [kpc]");
            global_hist.RZVPhi->GetYaxis()->SetTitle("Z [kpc]");
            global_hist.RZVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
            SetView3D(global_hist.RZVPhi);*/
        }

        // the R-Z-VZ plot
        {
          /*global_hist.RZVZ->SetStats(false);
          global_hist.RZVZ->GetXaxis()->SetTitle("R [kpc]");
          global_hist.RZVZ->GetYaxis()->SetTitle("Z [kpc]");
          global_hist.RZVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
          SetView3D(global_hist.RZVZ);*/
        }
      }

      {
        // purge the results
        std::lock_guard<std::mutex> res_lock(results_mtx);
        results.erase(uuid);
      }

      /*delete global_hist.XYVR;
      delete global_hist.XYVPhi;
      delete global_hist.XYVZ;
      delete global_hist.RZVR;
      delete global_hist.RZVPhi;
      delete global_hist.RZVZ;*/

      {
        // purge the requests
        std::lock_guard<std::mutex> lock(requests_mtx);
        requests.erase(uuid);
      }
    }).detach();
  }

  // respond with the dataset id
  if (offline) {
    std::string html = uuid;

    header_map mime;
    mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                     {"text/plain", false}));
    res->write_head(200, mime);
    res->end(html);
    return;
  }

  // the result-set already exists, prepare a full response

  std::string html =
      "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n";
  html.append(
      "<link href=\"https://fonts.googleapis.com/css?family=Inconsolata\" "
      "rel=\"stylesheet\"/>\n");
  html.append(
      "<link href=\"https://fonts.googleapis.com/css?family=Material+Icons\" "
      "rel=\"stylesheet\"/>\n");
  html.append("<script "
              "src=\"//cdnjs.cloudflare.com/ajax/libs/numeral.js/2.0.6/"
              "numeral.min.js\"></script>\n");

  // plotly
  /*html.append(
      "<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>");*/

  // bootstrap
  html.append(
      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, "
      "user-scalable=no, minimum-scale=1, maximum-scale=1\">\n");
  html.append("<link rel=\"stylesheet\" "
              "href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/"
              "bootstrap.min.css\">\n");
  html.append("<script "
              "src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/"
              "jquery.min.js\"></script>\n");
  html.append("<script "
              "src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/"
              "bootstrap.min.js\"></script>\n");

  // CERN JSROOT
  /*std::string url =
      "https://root.cern/js/latest/scripts/JSRootCore.min.js?more2d&3d&io&mathjax&onload=";*/

  // an official development version
  /*std::string url = "https://root.cern/js/dev/scripts/"
                    "JSRootCore.min.js?more2d&3d&io&mathjax&onload=";*/

  // an early development version
  /*std::string url = "https://jsroot.gsi.de/dev/scripts/"
                    "JSRootCore.min.js?more2d&3d&io&mathjax&onload=";*/

  // serve JSROOT from a local copy
  //std::string url = "scripts/JSRootCore.js?more2d&3d&io&mathjax&onload=";

  // serve JSROOT from CDN (jsdelivr)
  std::string url = "https://cdn.jsdelivr.net/gh/jvo203/GAIAWebQL/htdocs/gaiawebql/scripts/JSRootCore.js?more2d&3d&io&onload=";

  if (!exists)
    html.append("<script defer type=\"text/javascript\" "
                "src=\"" +
                url + "onloaded\"></script>");
  else
    html.append("<script defer type=\"text/javascript\" "
                "src=\"" +
                url + "fetch_plots\"></script>"); // was fetch_plots()

  // GAIAWebQL main JavaScript + CSS
  html.append("<script src=\"gaiawebql.js?" VERSION_STRING "\"></script>\n");
  html.append("<link rel=\"stylesheet\" href=\"gaiawebql.css?" VERSION_STRING
              "\"/>\n");

  // HTML content
  html.append("<title>GAIA DR2 WebQL</title></head><body>\n");
  html.append("<div id='session-data' style='width: 0; height: 0;' "
              "data-search-count='" +
              std::to_string(db_index.size()) + "' ");
  html.append("data-uuid='" + uuid + "' ");
  html.append("data-where='" + where + "' ");
  html.append("data-server-version='" + std::string(VERSION_STRING) +
              "' data-server-string='" + SERVER_STRING + "'></div>\n");

  // add the main HTML container element
  html.append("<div id="
              "main"
              " class="
              "container"
              ">\n");

  // html.append("<h1>GAIA DR2 WebQL</h1>");
  if (!exists) {
    html.append(
        "<h3 id=\"processing\">processing request; #HEALPix search pixels: " +
        std::to_string(db_index.size()) + "</h3>");
    html.append("<h3 id="
                "completed"
                "></h3>");
    html.append("<div class=\"progress\">");
    html.append(
        "<div id=\"progress-bar\" class=\"progress-bar progress-bar-info "
        "progress-bar-striped\" role=\"progressbar\" aria-valuenow=0 "
        "aria-valuemin=0 aria-valuemax=" +
        std::to_string(db_index.size()) +
        " style=\"width:0%\">0/0</div></div>");
    html.append("</div>");    
  }
  else
    html.append("<div><h3 id="
              "fetching"
              ">Fetching the plots, please wait...</h3></div>");

  html.append("<div><h3 id="
              "no-data"
              "></h3></div>");

  html.append("<div class=\"plots\">");

  html.append("<div id=\"hr\" style=\"width: 1000px; height: 600px\"></div>");

  html.append("<div id=\"mg\"></div><hr>");

  html.append(
      "<div id=\"xy\" style=\"width: 1000px; height: 600px\"></div><hr>");

  html.append(
      "<div id=\"rz\" style=\"width: 1000px; height: 600px\"></div><hr>");

  html.append("<div id=\"xyvr\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"xyvr_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"xyvr_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div><hr>");

  html.append("<div id=\"xyvphi\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"xyvphi_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"xyvphi_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div><hr>");

  html.append("<div id=\"xyvz\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"xyvz_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"xyvz_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div><hr>");

  html.append("<div id=\"rzvr\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"rzvr_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"rzvr_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div><hr>");

  html.append("<div id=\"rzvphi\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"rzvphi_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"rzvphi_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div><hr>");

  html.append("<div id=\"rzvz\" style=\"width: 1200px; height: "
              "500px\"><table><tr><td><div "
              "id=\"rzvz_mean\" style=\"width: 600px; height: "
              "500px\"></div></td><td><div "
              "id=\"rzvz_error\" style=\"width: 600px; height: "
              "500px\"></div></td></tr></table></div>");

  // end of plots
  html.append("</div>");

  /*html.append("<hr><h3>TEST</h3>");
  html.append(
      "<div id=\"large\" style=\"width: 1000px; height: 600px\"></div>");
  html.append("<div id=\"small\" style=\"width: 500px; height:
  300px\"></div>");*/

  // end of the main div
  html.append("</div>");

  html.append("<script>main();</script>");

  if (!exists) {
    html.append("<script>poll_progress();</script>");
    // html.append("<script>check_progress();</script>");
    // html.append("<script>console.log('this is a test');</script>");
  }

  html.append("</body></html>");

  header_map mime;
  mime.insert(std::pair<std::string, header_value>("Content-Type",
                                                   {"text/html", false}));
  res->write_head(200, mime);
  res->end(html);
}

char *base64(const unsigned char *input, int length) {
  BIO *bmem, *b64;
  BUF_MEM *bptr;

  b64 = BIO_new(BIO_f_base64());
  bmem = BIO_new(BIO_s_mem());
  b64 = BIO_push(b64, bmem);
  BIO_write(b64, input, length);
  BIO_flush(b64);
  BIO_get_mem_ptr(b64, &bptr);

  char *buff = (char *)malloc(bptr->length);
  memcpy(buff, bptr->data, bptr->length - 1);
  buff[bptr->length - 1] = 0;

  BIO_free_all(b64);

  return buff;
}

/*void SetView2D(TH2 *h) {
  auto mean_x = h->GetMean(1);
  auto std_x = h->GetStdDev(1);
  auto mean_y = h->GetMean(2);
  auto std_y = h->GetStdDev(2);

  h->GetXaxis()->SetRangeUser(mean_x - 10.0 * std_x, mean_x + 10.0 * std_x);
  h->GetYaxis()->SetRangeUser(mean_y - 10.0 * std_y, mean_y + 10.0 * std_y);
}

void SetView3D(TH3 *h) {
  auto mean_x = h->GetMean(1);
  auto std_x = h->GetStdDev(1);
  auto mean_y = h->GetMean(2);
  auto std_y = h->GetStdDev(2);
  auto mean_z = h->GetMean(3);
  auto std_z = h->GetStdDev(3);

  h->GetXaxis()->SetRangeUser(mean_x - 10.0 * std_x, mean_x + 10.0 * std_x);
  h->GetYaxis()->SetRangeUser(mean_y - 10.0 * std_y, mean_y + 10.0 * std_y);
  h->GetZaxis()->SetRangeUser(mean_z - 10.0 * std_z, mean_z + 10.0 * std_z);
}*/

int main(int argc, char *argv[]) {
  curl_global_init(CURL_GLOBAL_ALL);
  ROOT::EnableThreadSafety();

  // load the db healpix index file
  load_db_index("gaiadr2-table.dat");

  // parse local command-line options
  if (argc > 2) {
    for (int i = 1; i < argc - 1; i++) {
      const char *key = argv[i];
      const char *value = argv[i + 1];

      if (!strcmp(key, "--port"))
        server_port = atoi(value);
    }
  }

  std::cout << SERVER_STRING << " (" << VERSION_STRING << ")" << std::endl;
  std::cout << "Browser URL: https://localhost:" << server_port << std::endl;

  boost::system::error_code ec;
  boost::asio::ssl::context tls(boost::asio::ssl::context::sslv23);

  tls.use_private_key_file("ssl/server.key", boost::asio::ssl::context::pem);
  tls.use_certificate_chain_file("ssl/server.crt");

  if (configure_tls_context_easy(ec, tls)) {
    std::cerr << "error: " << ec.message() << std::endl;
  }

  int no_threads = MAX(std::thread::hardware_concurrency() / 2, 1);

  http2_server = new http2();
  http2_server->num_threads(no_threads);

  http2_server->handle("/", [](const request &req, const response &res) {
    boost::system::error_code ec;

    auto uri = req.uri().path;

    if (uri == "/") {
      auto push = res.push(ec, "GET", "/favicon.ico");
      serve_file(&req, push, "/favicon.ico");

      push = res.push(ec, "GET", "/index.css");
      serve_file(&req, push, "/index.css");

      push = res.push(ec, "GET", "/index.js");
      serve_file(&req, push, "/index.js");

      push = res.push(ec, "GET", "/paper_texture_08.png");
      serve_file(&req, push, "/paper_texture_08.png");

      serve_file(&req, &res, "/index.html");
    } else {
      std::cout << uri << std::endl;

      // status
      if (uri.find("status/") != std::string::npos) {
        size_t pos = uri.find_last_of("/");

        if (pos != std::string::npos) {
          std::string uuid = uri.substr(pos + 1, std::string::npos);

          // process the response
          std::cout << "status(" << uuid << ")" << std::endl;

          // find out if the uuid has already been processed
          std::string dir = docs_root + "/gaiawebql/DATA/" + uuid;

          if (std::filesystem::exists(dir))
            return http_ok(&res);
          else
            return http_accepted(&res);

        } else
          return http_not_found(&res);
      }

      // progress
      if (uri.find("progress/") != std::string::npos) {
        size_t pos = uri.find_last_of("/");

        if (pos != std::string::npos) {
          std::string uuid = uri.substr(pos + 1, std::string::npos);

          // process the response
          std::cout << "progress(" << uuid << ")" << std::endl;

          try {
            auto entry = requests.at(uuid);
            std::lock_guard<std::mutex> lock(entry->completed_mtx);
            // std::shared_lock<std::shared_mutex> lock(fits->progress_mtx);

            size_t len = entry->completed.size();

            // send a progress notification
            std::ostringstream json;
            json << "{ \"type\" : \"progress\",  \"completed\" : " << len
                 << ", \"total\" : " << db_index.size()
                 << ", \"exists\" : false, \"elapsed\" : "
                 << (std::time(nullptr) - entry->timestamp) << " }";

            header_map mime;
            mime.insert(std::pair<std::string, header_value>(
                "Content-Type", {"application/json", false}));
            res.write_head(200, mime);
            res.end(json.str());
            return;
          } catch (const std::out_of_range &err) {
            printf("no entry found for a job request %s\n", uuid.c_str());

            // find out if the uuid has already been processed
            if (dataset_exists(uuid)) {
              // send a progress notification
              std::ostringstream json;
              json << "{ \"type\" : \"progress\",  \"completed\" : "
                   << db_index.size() << ", \"total\" : " << db_index.size()
                   << ", \"exists\" : true, \"elapsed\" : 0 }";

              header_map mime;
              mime.insert(std::pair<std::string, header_value>(
                  "Content-Type", {"application/json", false}));
              res.write_head(200, mime);
              res.end(json.str());
              return;
            } else
              return http_no_content(&res);
          }
        } else
          return http_not_found(&res);
      }

      // GAIAWebQL entry
      if (uri.find("GAIAWebQL.html") != std::string::npos) {
        auto push = res.push(ec, "GET", "/favicon.ico");
        serve_file(&req, push, "/favicon.ico");

        /*push = res.push(ec, "GET", "/gaiawebql/gaiawebql.js");
        serve_file(&req, push, "/gaiawebql/gaiawebql.js");*/

        push = res.push(ec, "GET", "/gaiawebql/gaiawebql.css?" VERSION_STRING);
        serve_file(&req, push, "/gaiawebql/gaiawebql.css");

        push = res.push(ec, "GET", "/gaiawebql/gaiawebql.js?" VERSION_STRING);
        serve_file(&req, push, "/gaiawebql/gaiawebql.js");

        auto uri = req.uri();
        auto query = percent_decode(uri.raw_query);
        std::cout << "query: (" << query << ")" << std::endl;

        bool offline = false;
        std::string where;
        auto search = std::make_shared<struct search_criteria>();
        bool valid_params = false;

        std::vector<std::string> params;
        boost::split(params, query, [](char c) { return c == '&'; });

        for (auto const &s : params) {
          // check the value-less 'offline' flag
          if (s == "offline")
            offline = true;

          // find '='
          size_t pos = s.find("=");

          if (pos != std::string::npos) {
            std::string key = s.substr(0, pos);
            std::string value = s.substr(pos + 1, std::string::npos);

            if (key == "xmin") {
              char *e;
              errno = 0;

              search->X_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->X_min = NAN;
              else
                valid_params = true;
            }

            if (key == "xmax") {
              char *e;
              errno = 0;

              search->X_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->X_max = NAN;
              else
                valid_params = true;
            }

            if (key == "ymin") {
              char *e;
              errno = 0;

              search->Y_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Y_min = NAN;
              else
                valid_params = true;
            }

            if (key == "ymax") {
              char *e;
              errno = 0;

              search->Y_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Y_max = NAN;
              else
                valid_params = true;
            }

            if (key == "zmin") {
              char *e;
              errno = 0;

              search->Z_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Z_min = NAN;
              else
                valid_params = true;
            }

            if (key == "zmax") {
              char *e;
              errno = 0;

              search->Z_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Z_max = NAN;
              else
                valid_params = true;
            }

            if (key == "rmin") {
              char *e;
              errno = 0;

              search->R_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->R_min = NAN;
              else
                valid_params = true;
            }

            if (key == "rmax") {
              char *e;
              errno = 0;

              search->R_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->R_max = NAN;
              else
                valid_params = true;
            }

            if (key == "phimin") {
              char *e;
              errno = 0;

              search->Phi_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Phi_min = NAN;
              else
                valid_params = true;
            }

            if (key == "phimax") {
              char *e;
              errno = 0;

              search->Phi_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->Phi_max = NAN;
              else
                valid_params = true;
            }

            if (key == "mgmin") {
              char *e;
              errno = 0;

              search->M_G_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->M_G_min = NAN;
              else
                valid_params = true;
            }

            if (key == "mgmax") {
              char *e;
              errno = 0;

              search->M_G_max = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->M_G_max = NAN;
              else
                valid_params = true;
            }

            if (key == "parallax_over_error") {
              char *e;
              errno = 0;

              search->parallax_error_min = std::strtod(value.c_str(), &e);
              if (*e != '\0' || // error, we didn't consume the entire
                                // string
                  errno != 0)   // error, overflow or underflow
                search->parallax_error_min = NAN;
              else
                valid_params = true;
            }

            if (key == "where") {
              CURL *curl = curl_easy_init();

              char *str =
                  curl_easy_unescape(curl, value.c_str(), value.length(), NULL);
              where = std::string(str);
              curl_free(str);

              curl_easy_cleanup(curl);
            }
          }
        }

        if (where != "")
          std::cout << " where " << where;
        std::cout << std::endl;

        if (valid_params) {
          print_search_criteria(search);

          std::string uuid = make_uuid(search, where);

          return execute_gaia(&res, search, where, uuid, offline);
        } else {
          return http_not_found(&res);
        }
      };

      // by default try to serve a file
      serve_file(&req, &res, uri);
    }
  });

  signal(SIGPIPE, SIG_IGN); // ignore SIGPIPE
  signal(SIGINT, signalHandler);
  signal(SIGTERM, signalHandler);

  if (http2_server->listen_and_serve(ec, tls, "0.0.0.0",
                                     std::to_string(server_port), true)) {
    std::cerr << "error: " << ec.message() << std::endl;
  }

  http2_server->join();
  delete http2_server;

  std::cout << "GAIAWebQL shutdown completed." << std::endl;
}
