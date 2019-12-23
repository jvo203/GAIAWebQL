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
#define VERSION_STRING "SV2019-12-11.0"

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
#include <uWS/uWS.h>

typedef uWS::WebSocket<uWS::SERVER> uWS_SERVER_CONNECTION;
std::unordered_map<std::string, uWS_SERVER_CONNECTION *> progress_list;
std::mutex progress_mtx;

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
  SeedH2 _hr;
  SeedH2 _xy;
  SeedH2 _rz;

  /*TH3 *XYVR;
  TH3 *XYVPhi;
  TH3 *XYVZ;
  TH3 *RZVR;
  TH3 *RZVPhi;
  TH3 *RZVZ;*/
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

  std::cout << "GAIAWebQL shutdown completed." << std::endl;

  exit(signum);
}

// resource not found
void http_not_found(uWS::HttpResponse *res) {
  const std::string not_found =
      "HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n";
  res->write(not_found.data(), not_found.length());
}

// server error
void http_internal_server_error(uWS::HttpResponse *res) {
  const std::string server_error =
      "HTTP/1.1 500 Internal Server Error\r\nContent-Length: 0\r\n\r\n";
  res->write(server_error.data(), server_error.length());
}

// request accepted but not ready yet
void http_accepted(uWS::HttpResponse *res) {
  const std::string accepted =
      "HTTP/1.1 202 Accepted\r\nContent-Length: 0\r\n\r\n";
  res->write(accepted.data(), accepted.length());
}

// functionality not implemented/not available
void http_not_implemented(uWS::HttpResponse *res) {
  const std::string not_implemented =
      "HTTP/1.1 501 Not Implemented\r\nContent-Length: 0\r\n\r\n";
  res->write(not_implemented.data(), not_implemented.length());
}

void write_status(uWS::HttpResponse *res, int code, std::string message) {
  std::string status =
      "HTTP/1.1 " + std::to_string(code) + " " + message + "\r\n";
  res->write(status.data(), status.length());
}

void write_content_length(uWS::HttpResponse *res, size_t length) {
  std::string content_length =
      "Content-Length: " + std::to_string(length) + "\r\n";
  res->write(content_length.data(), content_length.length());
}

void write_content_type(uWS::HttpResponse *res, std::string mime) {
  std::string content_type = "Content-Type: " + mime + "\r\n";
  res->write(content_type.data(), content_type.length());
}

void write_key_value(uWS::HttpResponse *res, std::string key,
                     std::string value) {
  std::string content_type = key + ": " + value + "\r\n";
  res->write(content_type.data(), content_type.length());
}

void serve_file(uWS::HttpResponse *res, std::string uri) {
  std::string resource = /*"htdocs" +*/ uri;

  // strip '?' from the requested file name
  size_t pos = resource.find("?");

  if (pos != std::string::npos)
    resource = resource.substr(0, pos);

  std::cout << "serving " << resource << std::endl;

  // mmap a disk resource
  int fd = -1;
  void *buffer = NULL;

  struct stat64 st;
  stat64(resource.c_str(), &st);
  long size = st.st_size;

  fd = open(resource.c_str(), O_RDONLY);

  if (fd != -1) {
    buffer = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);

    if (buffer != NULL) {
      write_status(res, 200, "OK");
      write_content_length(res, size);

      // detect mime-types
      size_t pos = resource.find_last_of(".");

      if (pos != std::string::npos) {
        std::string ext = resource.substr(pos + 1, std::string::npos);

        if (ext == "htm" || ext == "html")
          write_content_type(res, "text/html");

        if (ext == "txt")
          write_content_type(res, "text/plain");

        if (ext == "js")
          write_content_type(res, "application/javascript");

        if (ext == "ico")
          write_content_type(res, "image/x-icon");

        if (ext == "png")
          write_content_type(res, "image/png");

        if (ext == "gif")
          write_content_type(res, "image/gif");

        if (ext == "webp")
          write_content_type(res, "image/webp");

        if (ext == "jpg" || ext == "jpeg")
          write_content_type(res, "image/jpeg");

        if (ext == "bpg")
          write_content_type(res, "image/bpg");

        if (ext == "mp4")
          write_content_type(res, "video/mp4");

        if (ext == "hevc")
          write_content_type(res, "video/hevc");

        if (ext == "css")
          write_content_type(res, "text/css");

        if (ext == "pdf")
          write_content_type(res, "application/pdf");

        if (ext == "svg")
          write_content_type(res, "image/svg+xml");

        if (ext == "wasm")
          write_content_type(res, "application/wasm");
      }

      res->write("\r\n", 2);
      res->write((const char *)buffer, size);
      res->write("\r\n\r\n", 4);
    } else {
      perror("error mapping a file");
      http_not_found(res);
    }

    if (munmap(buffer, size) == -1)
      perror("un-mapping error");

    close(fd);
  } else
    http_not_found(res);
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
      sql.append(" and parallax_over_error > " +
                 std::to_string(search->parallax_error_min));

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

                /*double ra = std::strtod(PQgetvalue(res, i, 0), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  ra = NAN;

                double dec = std::strtod(PQgetvalue(res, i, 1), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  dec = NAN;

                double phot_g_mean_mag = std::strtod(PQgetvalue(res, i, 2), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  phot_g_mean_mag = NAN;

                double bp_rp = std::strtod(PQgetvalue(res, i, 3), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  bp_rp = NAN;

                double parallax = std::strtod(PQgetvalue(res, i, 4), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  parallax = NAN;

                double pmra = std::strtod(PQgetvalue(res, i, 5), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  pmra = NAN;

                double pmdec = std::strtod(PQgetvalue(res, i, 6), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  pmdec = NAN;

                double radial_velocity = std::strtod(PQgetvalue(res, i, 7), &e);
                if (*e != '\0' || // error, we didn't consume the entire string
                    errno != 0)   // error, overflow or underflow
                  radial_velocity = NAN;*/

                /*if (!std::isnan(ra) && !std::isnan(dec) &&
                    !std::isnan(phot_g_mean_mag) && !std::isnan(bp_rp) &&
                    !std::isnan(parallax) && !std::isnan(pmra) &&
                    !std::isnan(pmdec) && !std::isnan(radial_velocity))*/
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
                  if (!std::isnan(search->X_min))
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
                      data_ok = false;

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
    size_t len = entry->completed.size();
    abort_search = entry->abort;

    try {
      std::lock_guard<std::mutex> lock(progress_mtx);
      auto ws = progress_list.at(uuid);

      // send completed via websockets
      std::ostringstream json;
      json << "{ \"type\" : \"progress\",  \"completed\" : " << len
           << ", \"elapsed\" : " << (std::time(nullptr) - entry->timestamp)
           << " }";
      //((steady_clock::now() - entry->timestamp).count()) *
      // steady_clock::period::num /
      // static_cast<double>(steady_clock::period::den)
      ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
    } catch (const std::out_of_range &err) {
      printf("no websocket connection found for a job request %s\n",
             uuid.c_str());
    }
  } catch (const std::out_of_range &err) {
    printf("no entry found for a job request %s\n", uuid.c_str());
    abort_search = true;
  }

  msg << "processed " << count << " records." << std::endl;
  std::cout << msg.str();

  return abort_search;
}

void execute_gaia(uWS::HttpResponse *res,
                  std::shared_ptr<struct search_criteria> search,
                  std::string where, std::string uuid) {

  // check if processing a request is already under way
  {
    std::lock_guard<std::mutex> req_lock(requests_mtx);
    if (requests.find(uuid) != requests.end()) {
      // respond with the dataset id
      std::string html = uuid + " is being processed. Please check later.";
      size_t size = html.length();
      write_status(res, 202, "Accepted");
      write_content_length(res, size);
      write_content_type(res, "text/plain");
      res->write("\r\n", 2);
      res->write((const char *)html.data(), size);
      res->write("\r\n\r\n", 4);
      return;
    }
  }

  // find out if the uuid has already been processed
  std::string dir = "DATA/" + uuid;
  bool exists = false;

  if (std::filesystem::exists(dir))
    exists = true;

  if (!exists) {
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

      struct gaia_hist global_hist {};
      char name[255];

      global_hist._hr.set_title("Hertzsprung-Russell diagram");
      global_hist._xy.set_title("X-Y");
      global_hist._rz.set_title("R-Z");

      // 3D histograms
      /*sprintf(name, "%s/XYVR", uuid.c_str());
      global_hist.XYVR = new TH3D(name, "X-Y-V_{R}", 100, -0.1, 0.1, 100, -0.1,
                                  0.1, 100, -0.1, 0.1);
      global_hist.XYVR->SetCanExtend(TH1::kAllAxes);

      sprintf(name, "%s/XYVPhi", uuid.c_str());
      global_hist.XYVPhi = new TH3D(name, "X-Y-V_{\\Phi}", 100, -0.1, 0.1, 100,
                                    -0.1, 0.1, 100, -0.1, 0.1);
      global_hist.XYVPhi->SetCanExtend(TH1::kAllAxes);

      sprintf(name, "%s/XYVZ", uuid.c_str());
      global_hist.XYVZ = new TH3D(name, "X-Y-V_{Z}", 100, -0.1, 0.1, 100, -0.1,
                                  0.1, 100, -0.1, 0.1);
      global_hist.XYVZ->SetCanExtend(TH1::kAllAxes);

      sprintf(name, "%s/RZVR", uuid.c_str());
      global_hist.RZVR = new TH3D(name, "R-Z-V_{R}", 100, -0.1, 0.1, 100, -0.1,
                                  0.1, 100, -0.1, 0.1);
      global_hist.RZVR->SetCanExtend(TH1::kAllAxes);

      sprintf(name, "%s/RZVPhi", uuid.c_str());
      global_hist.RZVPhi = new TH3D(name, "R-Z-V_{\\Phi}", 100, -0.1, 0.1, 100,
                                    -0.1, 0.1, 100, -0.1, 0.1);
      global_hist.RZVPhi->SetCanExtend(TH1::kAllAxes);

      sprintf(name, "%s/RZVZ", uuid.c_str());
      global_hist.RZVZ = new TH3D(name, "R-Z-V_{Z}", 100, -0.1, 0.1, 100, -0.1,
                                  0.1, 100, -0.1, 0.1);
      global_hist.RZVZ->SetCanExtend(TH1::kAllAxes);*/

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
          // a custom solution
          global_hist._hr.update(data.bp_rp, data.M_G);
          global_hist._xy.update(data.X, data.Y);
          global_hist._rz.update(data.R, data.Z);

          // CERN ROOT
          /*global_hist.XYVR->Fill(data.X, data.Y, data.VR);
          global_hist.XYVPhi->Fill(data.X, data.Y, data.VPhi);
          global_hist.XYVZ->Fill(data.X, data.Y, data.VZ);
          global_hist.RZVR->Fill(data.R, data.Z, data.VR);
          global_hist.RZVPhi->Fill(data.R, data.Z, data.VPhi);
          global_hist.RZVZ->Fill(data.R, data.Z, data.VZ);*/
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
      std::cout << "a global queue length: " << queue.queue.size() << std::endl;

      if (!search_aborted) {
        // save the histograms to disk
        global_hist._hr.save(uuid, "hr");
        global_hist._xy.save(uuid, "xy");
        global_hist._rz.save(uuid, "rz");

        // rename the temporary dir to just "DATA/uuid"
        std::string tmp = "DATA/" + uuid + ".tmp";
        std::string dir = "DATA/" + uuid;
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
      } else {
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

      std::lock_guard<std::mutex> lock(requests_mtx);
      requests.erase(uuid);
    }).detach();

    // respond with the dataset id
    std::string html = uuid;
    size_t size = html.length();
    write_status(res, 200, "OK");
    write_content_length(res, size);
    write_content_type(res, "text/plain");
    res->write("\r\n", 2);
    res->write((const char *)html.data(), size);
    res->write("\r\n\r\n", 4);
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
  html.append("<script "
              "src=\"https://cdn.jsdelivr.net/gh/jvo203/fits_web_ql/htdocs/"
              "fitswebql/ra_dec_conversion.js\"></script>\n");

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

  // CERN ROOT JS
  html.append("<script type=\"text/javascript\" "
              "src=\"https://root.cern/js/latest/scripts/"
              "JSRootCore.min.js?hist&onload=main\"></script>");

  // GAIAWebQL main JavaScript + CSS
  html.append("<script src=\"gaiawebql.js?" VERSION_STRING
              "\" defer></script>\n");
  // html.append("<link rel=\"stylesheet\" href=\"gaiawebql.css?" VERSION_STRING
  // "\"/>\n");

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
  html.append("<h1>GAIA DR2 WebQL</h1>");
  html.append("</div");
  html.append("</body></html>");

  size_t size = html.length();
  write_status(res, 200, "OK");
  write_content_length(res, size);
  write_content_type(res, "text/html");
  res->write("\r\n", 2);
  res->write((const char *)html.data(), size);
  res->write("\r\n\r\n", 4);
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

  // load the db healpix index file
  load_db_index("gaiadr2-table.dat");

  // register signal SIGINT and signal handler
  signal(SIGINT, signalHandler);

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
  std::cout << "Browser URL: http://localhost:" << server_port << std::endl;

  int no_threads = MIN(MAX(std::thread::hardware_concurrency() / 2, 1), 4);
  std::vector<std::thread *> threads(no_threads);
  std::transform(
      threads.begin(), threads.end(), threads.begin(), [](std::thread *t) {
        return new std::thread([]() {
          uWS::Hub h;

          h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req,
                             char *data, size_t, size_t) {
            std::string uri = req.getUrl().toString();

            std::cout << "HTTP request for " << uri << std::endl;

            // root
            if (uri == "/")
              return serve_file(res, "htdocs/index.html");

            // serve plot data
            if (uri.find("/gaiawebql/DATA") != std::string::npos) {
              size_t pos = uri.find("DATA");
              if (pos != std::string::npos) {
                std::string file = uri.substr(pos, std::string::npos);
                return serve_file(res, file);
              }
            }

            // GAIAWebQL entry
            if (uri.find("GAIAWebQL.html") != std::string::npos) {
              // get a position of '?'
              size_t pos = uri.find("?");

              if (pos != std::string::npos) {
                double ra = NAN;
                double dec = NAN;
                double radius = NAN;
                std::string where;

                auto search = std::make_shared<struct search_criteria>();
                bool valid_params = false;

                // using std::string for now as std::string_view is broken
                // in the Intel C++ compiler v19 Update 1
                //(works OK with v19 Update 2)
                // LLVM CLANG works OK with std::string_view

                std::string query = uri.substr(pos + 1, std::string::npos);
                std::cout << "query: (" << query << ")" << std::endl;

                std::vector<std::string> params;
                boost::split(params, query, [](char c) { return c == '&'; });

                for (auto const &s : params) {
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

                      search->parallax_error_min =
                          std::strtod(value.c_str(), &e);
                      if (*e != '\0' || // error, we didn't consume the entire
                                        // string
                          errno != 0)   // error, overflow or underflow
                        search->parallax_error_min = NAN;
                      else
                        valid_params = true;
                    }

                    if (key == "ra") {
                      char *e;
                      errno = 0;

                      // ra = atof(value.c_str());
                      ra = std::strtod(value.c_str(), &e);
                      if (*e != '\0' || // error, we didn't consume the entire
                                        // string
                          errno != 0)   // error, overflow or underflow
                        ra = NAN;
                    }

                    if (key == "dec") {
                      char *e;
                      errno = 0;

                      // dec = atof(value.c_str());
                      dec = std::strtod(value.c_str(), &e);
                      if (*e != '\0' || // error, we didn't consume the entire
                                        // string
                          errno != 0)   // error, overflow or underflow
                        dec = NAN;
                    }

                    if (key == "radius") {
                      char *e;
                      errno = 0;

                      // radius = atof(value.c_str());
                      radius = std::strtod(value.c_str(), &e);
                      if (*e != '\0' || // error, we didn't consume the entire
                                        // string
                          errno != 0)   // error, overflow or underflow
                        radius = NAN;
                    }

                    if (key == "where") {
                      CURL *curl = curl_easy_init();

                      char *str = curl_easy_unescape(curl, value.c_str(),
                                                     value.length(), NULL);
                      where = std::string(str);
                      curl_free(str);

                      curl_easy_cleanup(curl);
                    }
                  }
                }

                std::cout << "ra:" << ra << ", dec:" << dec
                          << ", radius:" << radius;

                if (where != "")
                  std::cout << " where " << where;
                std::cout << std::endl;

                if (valid_params) {
                  print_search_criteria(search);

                  std::size_t id = std::hash<std::string>{}(uri);

                  std::stringstream sstream;
                  sstream << std::hex << id;
                  std::string uuid = sstream.str();

                  return execute_gaia(res, search, where, uuid);
                } else {
                  const std::string not_found =
                      "ERROR: please specify valid search parameters";
                  res->end(not_found.data(), not_found.length());
                  return;
                }
              } else {
                const std::string not_found =
                    "ERROR: URL parameters not found.";
                res->end(not_found.data(), not_found.length());
                return;
              }
            }

            return serve_file(res, "htdocs" + uri);
          });

          // This makes use of the SO_REUSEPORT of the Linux kernel
          // Other solutions include listening to one port per thread
          // with or without some kind of proxy inbetween
          if (!h.listen(server_port, nullptr, uS::ListenOptions::REUSE_PORT)) {
            std::cout << "Failed to listen\n";
          }

          std::cout << "Launching a uWS::HTTP/WS thread\n";

          h.run();
        });
      });

  std::for_each(threads.begin(), threads.end(),
                [](std::thread *t) { t->join(); });
}
