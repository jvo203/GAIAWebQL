#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_SUB 0

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

static const double halfpi = 1.570796326794896619231321691639751442099;
#define HTMfunc_Pr 3.1415926535897932385E0 / 180.0

/*#include <math.h>

const double R[3][3] = {{-0.0548755604162154, -0.8734370902348850, -0.4838350155487132}, {+0.4941094278755837, -0.4448296299600112, +0.7469822444972189}, {-0.8676661490190047, -0.1980763734312015, +0.4559837761750669}};
const double _theta = 0.003253017785496385;
const double H[3][3] = {{cos(_theta), 0.0, sin(_theta)}, {0.0, 1.0, 0.0}, {-sin(_theta), 0.0, cos(_theta)}};
const double dGC = 8300.0;*/

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define SERVER_PORT 8081
#define SERVER_STRING "GAIAWebQL v" STR(VERSION_MAJOR) "." STR(VERSION_MINOR) "." STR(VERSION_SUB)
#define VERSION_STRING "SV2019-06-27.0"

#include <sys/types.h>
#include <pwd.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>

//base64 encoding with SSL
#include <openssl/sha.h>
#include <openssl/hmac.h>
#include <openssl/evp.h>
#include <openssl/bio.h>
#include <openssl/buffer.h>

#include "PJMCoords.h"

char *base64(const unsigned char *input, int length);

#include <iostream>
#include <fstream>
#include <thread>
#include <unordered_map>
#include <ctime>
#include <chrono>
using std::chrono::steady_clock;

#include <curl/curl.h>
#include <uuid/uuid.h>
#include <libpq-fe.h>
#include "healpix_base.h"
#include "json.h"

//CERN ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TH2.h>
#include <TH3.h>
#include <TH2D.h>
#include <TImage.h>
#include <TCanvas.h>
#include <TGaxis.h>
//#include <TThread.h>

void ReverseYAxis(TH1 *h);
void ReverseYData(TH2 *h);
void SetView2D(TH2 *h);
void SetView3D(TH3 *h);

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lockfree/stack.hpp>
#include <boost/atomic.hpp>

struct search_criteria
{
    double X_min, X_max;
    double Y_min, Y_max;
    double Z_min, Z_max;
    double R_min, R_max;
    double Phi_min, Phi_max;
    double M_G_min, M_G_max;
    double parallax_error_min;

    search_criteria()
    {
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

void print_search_criteria(struct search_criteria *search)
{
    printf("X: [%f - %f]\n", search->X_min, search->X_max);
    printf("Y: [%f - %f]\n", search->Y_min, search->Y_max);
    printf("Z: [%f - %f]\n", search->Z_min, search->Z_max);
    printf("R: [%f - %f]\n", search->R_min, search->R_max);
    printf("Φ: [%f - %f]\n", search->Phi_min, search->Phi_max);
    printf("M_G: [%f - %f]\n", search->M_G_min, search->M_G_max);
    printf("parallax_error: ≥%f\n", search->parallax_error_min);
}

struct plot_data
{
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

#include <omp.h>
#include <uWS/uWS.h>

typedef uWS::WebSocket<uWS::SERVER> uWS_SERVER_CONNECTION;
std::unordered_map<std::string, uWS_SERVER_CONNECTION *> progress_list;
std::mutex progress_mtx;

int server_port = SERVER_PORT;

struct db_entry
{
    std::string schema_name;
    std::string table_name;
    std::string owner;
    std::string host;
    int port;
};

std::unordered_map<int, struct db_entry> db_index;

struct db_search_job
{
    //steady_clock::time_point timestamp = steady_clock::now();
    time_t timestamp = std::time(nullptr);
    std::vector<int> completed;
    std::mutex completed_mtx;
    bool abort = false;
};

std::unordered_map<std::string, std::shared_ptr<struct db_search_job>> requests;
std::mutex requests_mtx;

struct gaia_plots
{
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

struct gaia_hist
{
    TH2 *HR;
    TH2 *XY;
    TH2 *RZ;
    TH3 *XYVR;
    TH3 *XYVPhi;
    TH3 *XYVZ;
    TH3 *RZVR;
    TH3 *RZVPhi;
    TH3 *RZVZ;
};

std::mutex root_mtx;

inline const char *
check_null(const char *str)
{
    if (str != nullptr)
        return str;
    else
        return "\"\"";
};

void signalHandler(int signum)
{
    std::cout << "Interrupt signal (" << signum << ") received.\n";

    curl_global_cleanup();

    // cleanup and close up stuff here
    // terminate program

    std::cout << "remaining requests: " << requests.size() << "\tremaining results: " << results.size() << std::endl;

    std::cout << "GAIAWebQL shutdown completed." << std::endl;

    exit(signum);
}

//resource not found
void http_not_found(uWS::HttpResponse *res)
{
    const std::string not_found = "HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n";
    res->write(not_found.data(), not_found.length());
}

//server error
void http_internal_server_error(uWS::HttpResponse *res)
{
    const std::string server_error = "HTTP/1.1 500 Internal Server Error\r\nContent-Length: 0\r\n\r\n";
    res->write(server_error.data(), server_error.length());
}

//request accepted but not ready yet
void http_accepted(uWS::HttpResponse *res)
{
    const std::string accepted = "HTTP/1.1 202 Accepted\r\nContent-Length: 0\r\n\r\n";
    res->write(accepted.data(), accepted.length());
}

//functionality not implemented/not available
void http_not_implemented(uWS::HttpResponse *res)
{
    const std::string not_implemented = "HTTP/1.1 501 Not Implemented\r\nContent-Length: 0\r\n\r\n";
    res->write(not_implemented.data(), not_implemented.length());
}

void write_status(uWS::HttpResponse *res, int code, std::string message)
{
    std::string status = "HTTP/1.1 " + std::to_string(code) + " " + message + "\r\n";
    res->write(status.data(), status.length());
}

void write_content_length(uWS::HttpResponse *res, size_t length)
{
    std::string content_length = "Content-Length: " + std::to_string(length) + "\r\n";
    res->write(content_length.data(), content_length.length());
}

void write_content_type(uWS::HttpResponse *res, std::string mime)
{
    std::string content_type = "Content-Type: " + mime + "\r\n";
    res->write(content_type.data(), content_type.length());
}

void write_key_value(uWS::HttpResponse *res, std::string key, std::string value)
{
    std::string content_type = key + ": " + value + "\r\n";
    res->write(content_type.data(), content_type.length());
}

void serve_file(uWS::HttpResponse *res, std::string uri)
{
    std::string resource = "htdocs" + uri;

    //strip '?' from the requested file name
    size_t pos = resource.find("?");

    if (pos != std::string::npos)
        resource = resource.substr(0, pos);

    std::cout << "serving " << resource << std::endl;

    //mmap a disk resource
    int fd = -1;
    void *buffer = NULL;

    struct stat64 st;
    stat64(resource.c_str(), &st);
    long size = st.st_size;

    fd = open(resource.c_str(), O_RDONLY);

    if (fd != -1)
    {
        buffer = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);

        if (buffer != NULL)
        {
            write_status(res, 200, "OK");
            write_content_length(res, size);

            //detect mime-types
            size_t pos = resource.find_last_of(".");

            if (pos != std::string::npos)
            {
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
        }
        else
        {
            perror("error mapping a file");
            http_not_found(res);
        }

        if (munmap(buffer, size) == -1)
            perror("un-mapping error");

        close(fd);
    }
    else
        http_not_found(res);
}

void load_db_index(std::string filename)
{
    std::ifstream index_file(filename);

    for (std::string line; getline(index_file, line);)
    {
        std::cout << line << std::endl;

        std::vector<std::string> columns;
        boost::split(columns, line, [](char c) { return c == '|'; });

        //extract the database connection info based on the healpix index
        if (columns.size() == 6)
        {
            int hpx;
            struct db_entry entry;

            entry.schema_name = columns[0];
            entry.table_name = columns[1];
            entry.owner = columns[3];
            entry.host = columns[4];
            entry.port = std::stoi(columns[5]);

            sscanf(entry.table_name.c_str(), "gaia_source_%d", &hpx);

            db_index.insert(std::make_pair(hpx, entry));
        }
    }

    std::cout << "PostgreSQL HEALPix index contains " << db_index.size() << " entries." << std::endl;
}

bool search_gaia_db(int hpx, struct db_entry &entry, std::string uuid, struct search_criteria *search, std::string where, OmniCoords &coords, boost::lockfree::stack<struct plot_data> &plot_stack, struct gaia_hist &hist)
{
    std::stringstream msg;
    bool abort_search = false;

    msg << uuid << ":\t" << entry.schema_name << "/" << entry.table_name << "/" << entry.owner << "/" << entry.host << ":" << entry.port << "\t";

    std::string conn_str = "dbname=gaiadr2 host=" + entry.host + " port=" + std::to_string(entry.port) + " user=" + entry.owner + " password=jvo!";

    PGconn *gaia_db = PQconnectdb(conn_str.c_str());
    uint64 count = 0;

    if (PQstatus(gaia_db) != CONNECTION_OK)
    {
        fprintf(stderr, "PostgreSQL connection failed: %s\n", PQerrorMessage(gaia_db));

        PQfinish(gaia_db);
        gaia_db = NULL;
    }
    else
    {
        printf("PostgreSQL connection successful.\n");

        //perform a cone search
        std::string sql = "select ra,dec,phot_g_mean_mag,bp_rp,parallax,pmra,pmdec,radial_velocity from " + entry.schema_name + "." + entry.table_name + " where parallax > 0"; // and parallax_over_error > 10;"; //" limit 1;";

        //add optional search conditions
        if (where != "")
            sql.append(" and " + where);

        //parallax error
        if (!std::isnan(search->parallax_error_min))
            sql.append(" and parallax_over_error > " + std::to_string(search->parallax_error_min));

        sql.append(";"); //finish the sql

        //std::cout << sql << std::endl;

        if (PQsendQuery(gaia_db, sql.c_str()))
        {
            if (PQsetSingleRowMode(gaia_db))
            {
                PGresult *res = NULL;

                while ((res = PQgetResult(gaia_db)) != NULL)
                {
                    if (PQresultStatus(res) == PGRES_SINGLE_TUPLE)
                    {
                        count++;

                        std::stringstream res_str;
                        res_str << count << ":\t";

                        int nRows = PQntuples(res);
                        int nFields = PQnfields(res);

                        /*for (int i = 0; i < nFields; i++)
                            res_str << PQfname(res, i) << "\t";
                        res_str << std::endl;*/

                        for (int i = 0; i < nRows; i++)
                        {
                            if (nFields >= 8)
                            {
                                char *e;
                                errno = 0;

                                double ra = std::strtod(PQgetvalue(res, i, 0), &e);
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
                                    radial_velocity = NAN;

                                if (!std::isnan(ra) && !std::isnan(dec) && !std::isnan(phot_g_mean_mag) && !std::isnan(bp_rp) && !std::isnan(parallax) && !std::isnan(pmra) && !std::isnan(pmdec) && !std::isnan(radial_velocity))
                                {
                                    double M_G = phot_g_mean_mag + 5.0 + 5.0 * log10(parallax / 1000.0);

                                    //res_str << "ra: " << ra << "\t dec: " << dec << "phot_g_mean_mag = " << phot_g_mean_mag << "\tbp_rp = " << bp_rp << "\tparallax = " << parallax << "\tM_G = " << M_G << std::endl;

                                    double alpha = ra * HTMfunc_Pr;  //[rad]
                                    double delta = dec * HTMfunc_Pr; //[rad]
                                    double d = 1000.0 / parallax;    //distance [parsec, pc]

                                    /*double rICRS[3] = {d * cos(alpha) * cos(delta), d * sin(alpha) * cos(delta), d * sin(delta)};
                                    double rGC[3] = {0.0, 0.0, 0.0};
                                    double r[3] = {-dGC, 0.0, 0.0};

                                    for (int i = 0; i < 3; i++)
                                        for (int j = 0; j < 3; j++)
                                            r[i] += R[i][j] * rICRS[j];

                                    for (int i = 0; i < 3; i++)
                                        for (int j = 0; j < 3; j++)
                                            rGC[i] += H[i][j] * r[j];

                                    double x = rGC[0];
                                    double y = rGC[1];
                                    double z = rGC[2];*/

                                    vec6 sHEQ, sGCA, sGCY;
                                    sHEQ[0] = d / 1000.0;      //[kpc]
                                    sHEQ[1] = ra;              //[deg]
                                    sHEQ[2] = dec;             //[deg]
                                    sHEQ[3] = radial_velocity; //[km/s]
                                    sHEQ[4] = pmra;
                                    sHEQ[5] = pmdec;

                                    coords.take_HEQ_units(sHEQ);
                                    sGCA = coords.give_GCA_units();

                                    coords.take_HEQ_units(sHEQ);
                                    sGCY = coords.give_GCY_units();

                                    double X = sGCA[0]; //[kpc]
                                    double Y = sGCA[1]; //[kpc]
                                    double Z = sGCA[2]; //[kpc]

                                    double R = sGCY[0];    //[kpc]
                                    double Phi = sGCY[2];  //[rad]
                                    double VR = sGCY[3];   //[km/s]
                                    double VZ = sGCY[4];   //[km/s]
                                    double VPhi = sGCY[5]; //[km/s]

                                    bool data_ok = true;

                                    //validate data against the search critera
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

                                    //the data is OK, add values to histograms
                                    if (data_ok)
                                    {
                                        struct plot_data data = {bp_rp, M_G, X, Y, Z, R, VR, VPhi, VZ};
                                        while (!plot_stack.push(data))
                                            ;
                                    }
                                }
                            }
                        }
                    };

                    PQclear(res);
                };
            }
            else
                std::cout << "error setting PQsetSingleRowMode.\n";
        }
        else
            std::cout << "PQsendQuery error.\n";
    }

    if (gaia_db != NULL)
        PQfinish(gaia_db);

    try
    {
        auto entry = requests.at(uuid);
        std::lock_guard lock(entry->completed_mtx);
        entry->completed.push_back(hpx);
        size_t len = entry->completed.size();
        abort_search = entry->abort;

        try
        {
            std::lock_guard lock(progress_mtx);
            auto ws = progress_list.at(uuid);

            //send completed via websockets
            std::ostringstream json;
            json << "{ \"type\" : \"progress\",  \"completed\" : " << len << ", \"elapsed\" : " << (std::time(nullptr) - entry->timestamp) << " }";
            //((steady_clock::now() - entry->timestamp).count()) * steady_clock::period::num / static_cast<double>(steady_clock::period::den)
            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
        }
        catch (const std::out_of_range &err)
        {
            printf("no websocket connection found for a job request %s\n", uuid.c_str());
        }
    }
    catch (const std::out_of_range &err)
    {
        printf("no entry found for a job request %s\n", uuid.c_str());
        abort_search = true;
    }

    msg << "processed " << count << " records." << std::endl;
    std::cout << msg.str();

    return abort_search;
}

bool gaia_db_cone_search(int hpx, struct db_entry &entry, std::string uuid, double longitude, double latitude, double radius, std::string where, OmniCoords &coords, boost::lockfree::stack<struct plot_data> &plot_stack, struct gaia_hist &hist)
{
    std::stringstream msg;
    bool abort_search = false;

    msg << uuid << ":\t" << entry.schema_name << "/" << entry.table_name << "/" << entry.owner << "/" << entry.host << ":" << entry.port << "\t";

    std::string conn_str = "dbname=gaiadr2 host=" + entry.host + " port=" + std::to_string(entry.port) + " user=" + entry.owner + " password=jvo!";

    PGconn *gaia_db = PQconnectdb(conn_str.c_str());
    uint64 count = 0;

    if (PQstatus(gaia_db) != CONNECTION_OK)
    {
        fprintf(stderr, "PostgreSQL connection failed: %s\n", PQerrorMessage(gaia_db));

        PQfinish(gaia_db);
        gaia_db = NULL;
    }
    else
    {
        printf("PostgreSQL connection successful.\n");

        //perform a cone search
        std::string sql = "select ra,dec,phot_g_mean_mag,bp_rp,parallax,pmra,pmdec,radial_velocity,_point_radec from " + entry.schema_name + "." + entry.table_name + " where (_point_radec <-> spoint(" + std::to_string(longitude) + " , " + std::to_string(latitude) + ")) < " + std::to_string(radius) + " and parallax > 0"; // and parallax_over_error > 10;"; //" limit 1;";

        //add optional search conditions
        if (where != "")
            sql.append(" and " + where);

        sql.append(";"); //finish the sql

        //std::cout << sql << std::endl;

        if (PQsendQuery(gaia_db, sql.c_str()))
        {
            if (PQsetSingleRowMode(gaia_db))
            {
                PGresult *res = NULL;

                while ((res = PQgetResult(gaia_db)) != NULL)
                {
                    if (PQresultStatus(res) == PGRES_SINGLE_TUPLE)
                    {
                        count++;

                        std::stringstream res_str;
                        res_str << count << ":\t";

                        int nRows = PQntuples(res);
                        int nFields = PQnfields(res);

                        /*for (int i = 0; i < nFields; i++)
                            res_str << PQfname(res, i) << "\t";
                        res_str << std::endl;*/

                        for (int i = 0; i < nRows; i++)
                        {
                            if (nFields >= 8)
                            {
                                char *e;
                                errno = 0;

                                double ra = std::strtod(PQgetvalue(res, i, 0), &e);
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
                                    radial_velocity = NAN;

                                if (!std::isnan(ra) && !std::isnan(dec) && !std::isnan(phot_g_mean_mag) && !std::isnan(bp_rp) && !std::isnan(parallax) && !std::isnan(pmra) && !std::isnan(pmdec) && !std::isnan(radial_velocity))
                                {
                                    double M_G = phot_g_mean_mag + 5.0 + 5.0 * log10(parallax / 1000.0);

                                    //res_str << "ra: " << ra << "\t dec: " << dec << "phot_g_mean_mag = " << phot_g_mean_mag << "\tbp_rp = " << bp_rp << "\tparallax = " << parallax << "\tM_G = " << M_G << std::endl;

                                    double alpha = ra * HTMfunc_Pr;  //[rad]
                                    double delta = dec * HTMfunc_Pr; //[rad]
                                    double d = 1000.0 / parallax;    //distance [parsec, pc]

                                    /*double rICRS[3] = {d * cos(alpha) * cos(delta), d * sin(alpha) * cos(delta), d * sin(delta)};
                                    double rGC[3] = {0.0, 0.0, 0.0};
                                    double r[3] = {-dGC, 0.0, 0.0};

                                    for (int i = 0; i < 3; i++)
                                        for (int j = 0; j < 3; j++)
                                            r[i] += R[i][j] * rICRS[j];

                                    for (int i = 0; i < 3; i++)
                                        for (int j = 0; j < 3; j++)
                                            rGC[i] += H[i][j] * r[j];

                                    double x = rGC[0];
                                    double y = rGC[1];
                                    double z = rGC[2];*/

                                    vec6 sHEQ, sGCA, sGCY;
                                    sHEQ[0] = d / 1000.0;      //[kpc]
                                    sHEQ[1] = ra;              //[deg]
                                    sHEQ[2] = dec;             //[deg]
                                    sHEQ[3] = radial_velocity; //[km/s]
                                    sHEQ[4] = pmra;
                                    sHEQ[5] = pmdec;

                                    coords.take_HEQ_units(sHEQ);
                                    sGCA = coords.give_GCA_units();

                                    coords.take_HEQ_units(sHEQ);
                                    sGCY = coords.give_GCY_units();

                                    double X = sGCA[0]; //[kpc]
                                    double Y = sGCA[1]; //[kpc]
                                    double Z = sGCA[2]; //[kpc]

                                    double R = sGCY[0];    //[kpc]
                                    double VR = sGCY[3];   //[km/s]
                                    double VZ = sGCY[4];   //[km/s]
                                    double VPhi = sGCY[5]; //[km/s]

                                    //the data is OK, add values to histograms
                                    struct plot_data data = {bp_rp, M_G, X, Y, Z, R, VR, VPhi, VZ};
                                    while (!plot_stack.push(data))
                                        ;

                                    /*std::lock_guard lock(hist.plots_mtx);
                                    hist.HR->Fill(bp_rp, M_G);
                                    hist.XY->Fill(X, Y);
                                    hist.RZ->Fill(R, Z);
                                    hist.XYVR->Fill(X, Y, VR);
                                    hist.XYVPhi->Fill(X, Y, VPhi);
                                    hist.XYVZ->Fill(X, Y, VZ);
                                    hist.RZVR->Fill(R, Z, VR);
                                    hist.RZVPhi->Fill(R, Z, VPhi);
                                    hist.RZVZ->Fill(R, Z, VZ);*/
                                }
                            }

                            /*for (int j = 0; j < nFields; j++)
                                res_str << PQgetvalue(res, i, j) << "\t";
                            res_str << std::endl
                                    << std::endl;*/
                        }

                        //std::cout << res_str.str();
                    };

                    PQclear(res);
                };
            }
            else
                std::cout << "error setting PQsetSingleRowMode.\n";
        }
        else
            std::cout << "PQsendQuery error.\n";
    }

    if (gaia_db != NULL)
        PQfinish(gaia_db);

    try
    {
        auto entry = requests.at(uuid);
        std::lock_guard lock(entry->completed_mtx);
        entry->completed.push_back(hpx);
        size_t len = entry->completed.size();
        abort_search = entry->abort;

        try
        {
            std::lock_guard lock(progress_mtx);
            auto ws = progress_list.at(uuid);

            //send completed via websockets
            std::ostringstream json;
            json << "{ \"type\" : \"progress\",  \"completed\" : " << len << ", \"elapsed\" : " << (std::time(nullptr) - entry->timestamp) << " }";
            //((steady_clock::now() - entry->timestamp).count()) * steady_clock::period::num / static_cast<double>(steady_clock::period::den)
            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
        }
        catch (const std::out_of_range &err)
        {
            printf("no websocket connection found for a job request %s\n", uuid.c_str());
        }
    }
    catch (const std::out_of_range &err)
    {
        printf("no entry found for a job request %s\n", uuid.c_str());
        abort_search = true;
    }

    msg << "processed " << count << " records." << std::endl;
    std::cout << msg.str();

    return abort_search;
}

void execute_gaia(uWS::HttpResponse *res, struct search_criteria *search, std::string where)
{
    std::string uuid = []() -> std::string {
        uuid_t binuuid;
        uuid_generate(binuuid);
        char *uuid = (char *)malloc(37);

        if (uuid != NULL)
        {
            uuid_unparse(binuuid, uuid);
            std::string uuid_str(uuid);
            free(uuid);
            return uuid_str;
        }
        else
            return std::string("NULL");
    }();

    {
        std::lock_guard req_lock(requests_mtx);
        requests.insert(std::make_pair(uuid, std::make_shared<struct db_search_job>()));
    }

    {
        std::lock_guard res_lock(results_mtx);
        results.insert(std::make_pair(uuid, std::make_shared<struct gaia_plots>()));
    }

    std::thread([=]() {
        printf("starting a GAIA database search thread, job id:%s\n", uuid.c_str());

        int max_threads = omp_get_max_threads();
        //std::vector<std::shared_ptr<struct gaia_hist>> thread_hist; //(max_threads);
        std::vector<OmniCoords> thread_coords(max_threads);

        struct gaia_hist global_hist;
        char name[255];

        sprintf(name, "%s/HR", uuid.c_str());
        global_hist.HR = new TH2D(name, "Hertzsprung-Russell diagram", 600, -1.0, 5.0, 600, -5.0, 15.0);
        /*global_hist.hr_hist = new TH2D(uuid.c_str(), "Hertzsprung-Russell diagram", 600, -1.0, 5.0, 600, -1.0, 1.0);
        global_hist.hr_hist->SetCanExtend(TH1::kAllAxes);*/
        //global_hist.hr_hist->SetBit(TH1::kAutoBinPTwo);

        sprintf(name, "%s/XY", uuid.c_str());
        global_hist.XY = new TH2D(name, "X-Y", 600, -0.1, 0.1, 600, -0.1, 0.1);
        global_hist.XY->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZ", uuid.c_str());
        global_hist.RZ = new TH2D(name, "R-Z", 600, -0.1, 0.1, 600, -0.1, 0.1);
        global_hist.RZ->SetCanExtend(TH1::kAllAxes);

        //3D histograms
        sprintf(name, "%s/XYVR", uuid.c_str());
        global_hist.XYVR = new TH3D(name, "X-Y-V_{R}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVR->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/XYVPhi", uuid.c_str());
        global_hist.XYVPhi = new TH3D(name, "X-Y-V_{\\Phi}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVPhi->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/XYVZ", uuid.c_str());
        global_hist.XYVZ = new TH3D(name, "X-Y-V_{Z}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVZ->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVR", uuid.c_str());
        global_hist.RZVR = new TH3D(name, "R-Z-V_{R}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVR->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVPhi", uuid.c_str());
        global_hist.RZVPhi = new TH3D(name, "R-Z-V_{\\Phi}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVPhi->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVZ", uuid.c_str());
        global_hist.RZVZ = new TH3D(name, "R-Z-V_{Z}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVZ->SetCanExtend(TH1::kAllAxes);

        for (int i = 0; i < max_threads; i++)
        {
            thread_coords[i].change_sol_pos(8.3, 0.027);

            char name[255];
            sprintf(name, "%s/%d", uuid.c_str(), (i + 1));

            //struct gaia_hist hist;
            //hist.hr_hist = new TH2D(name, "M_{G} vs. BP-RP", 600, -1.0, 5.0, 600, -5.0, 15.0);
            /*hist->SetCanExtend(TH1::kAllAxes);
            hist->SetBit(TH1::kAutoBinPTwo);*/
            //thread_hist[i] = hist;
            //thread_hist.push_back(std::make_shared<struct gaia_hist>(std::move(hist)));
        }

        boost::lockfree::stack<struct plot_data> plot_stack(100000);
        boost::atomic<bool> search_done(false);

        std::thread plot_thread([&]() {
            std::cout << "starting a histogram thread function.\n";

            unsigned long counter = 0;
            struct plot_data data;

            auto plotter = [&](struct plot_data data) {
                global_hist.HR->Fill(data.bp_rp, data.M_G);
                global_hist.XY->Fill(data.X, data.Y);
                global_hist.RZ->Fill(data.R, data.Z);
                global_hist.XYVR->Fill(data.X, data.Y, data.VR);
                global_hist.XYVPhi->Fill(data.X, data.Y, data.VPhi);
                global_hist.XYVZ->Fill(data.X, data.Y, data.VZ);
                global_hist.RZVR->Fill(data.R, data.Z, data.VR);
                global_hist.RZVPhi->Fill(data.R, data.Z, data.VPhi);
                global_hist.RZVZ->Fill(data.R, data.Z, data.VZ);
            };

            //older Boost on py1 does not have consume_all!!!

            while (!search_done)
            {
                //counter += plot_stack.consume_all(plotter);
                while (plot_stack.pop(data))
                {
                    plotter(data);
                    counter++;
                }
            }

            //counter += plot_stack.consume_all(plotter);
            while (plot_stack.pop(data))
            {
                plotter(data);
                counter++;
            }

            std::cout << "a histogram thread ended after processing " << counter << " values.\n";
        });

        boost::atomic<bool> search_aborted(false);

#pragma omp parallel shared(global_hist)
        {
#pragma omp single
            {
                for (auto &it : db_index)
                    if (!search_aborted)
                    {
#pragma omp task
                        {
                            int hpx = it.first;
                            auto entry = it.second;

                            //for (size_t i = 0; i < pixels->size(); i++)
                            int tid = omp_get_thread_num();
                            bool abort_search = false;

                            printf("parametric search through hpx %d\n", hpx);

                            abort_search = search_gaia_db(hpx, entry, uuid, search, where, thread_coords[tid], plot_stack, global_hist);

                            if (abort_search)
                            {
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

        if (!search_aborted)
        {

            for (int i = 0; i < max_threads; i++)
            {
                //global_hist.hr_hist->Add(thread_hist[i]->hr_hist);
                //delete thread_hist[i]->hr_hist;
            }

            //the H-R diagram
            {
                ReverseYData(global_hist.HR);
                global_hist.HR->SetStats(false);
                global_hist.HR->GetXaxis()->SetTitle("(BP-RP) [mag]");
                global_hist.HR->GetYaxis()->SetTitle("M_{G} [mag]");
                //global_hist->GetZaxis()->SetTitle("star density");

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.HR->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);
                ReverseYAxis(global_hist.HR);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "HR diagram PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (1) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->hr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y plot
            {
                global_hist.XY->SetStats(false);
                global_hist.XY->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XY->GetYaxis()->SetTitle("Y [kpc]");
                SetView2D(global_hist.XY);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XY->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (2) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xy = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z plot
            {
                global_hist.RZ->SetStats(false);
                global_hist.RZ->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZ->GetYaxis()->SetTitle("Z [kpc]");
                SetView2D(global_hist.RZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZ->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (3) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VR plot
            {
                global_hist.XYVR->SetStats(false);
                global_hist.XYVR->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVR->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVR->GetZaxis()->SetTitle("V_{R} [km/s]");
                SetView3D(global_hist.XYVR);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVR->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VR plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (4) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VPhi plot
            {
                global_hist.XYVPhi->SetStats(false);
                global_hist.XYVPhi->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVPhi->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
                SetView3D(global_hist.XYVPhi);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVPhi->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VPhi plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (5) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvphi = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VZ plot
            {
                global_hist.XYVZ->SetStats(false);
                global_hist.XYVZ->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVZ->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
                SetView3D(global_hist.XYVZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVZ->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VZ plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (6) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VR plot
            {
                global_hist.RZVR->SetStats(false);
                global_hist.RZVR->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVR->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVR->GetZaxis()->SetTitle("V_{R} [km/s]");
                SetView3D(global_hist.RZVR);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVR->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VR plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (7) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VPhi plot
            {
                global_hist.RZVPhi->SetStats(false);
                global_hist.RZVPhi->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVPhi->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
                SetView3D(global_hist.RZVPhi);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVPhi->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VPhi plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (8) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvphi = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VZ plot
            {
                global_hist.RZVZ->SetStats(false);
                global_hist.RZVZ->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVZ->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
                SetView3D(global_hist.RZVZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVZ->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VZ plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (9) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }
        }
        else
        {
            //purge the results
            std::lock_guard res_lock(results_mtx);
            results.erase(uuid);
        }

        delete global_hist.HR;
        delete global_hist.XY;
        delete global_hist.RZ;
        delete global_hist.XYVR;
        delete global_hist.XYVPhi;
        delete global_hist.XYVZ;
        delete global_hist.RZVR;
        delete global_hist.RZVPhi;
        delete global_hist.RZVZ;

        std::lock_guard lock(requests_mtx);
        requests.erase(uuid);
    })
        .detach();

    /*std::string html = "GAIA DR2 WebQL::<ra:" + std::to_string(ra) + ", dec:" + std::to_string(dec) + ", radius:" + std::to_string(radius) + "> --> job id:" + uuid + ", hpx:" + std::to_string(idx) + ", searching pixels: ";

    for (int i = 0; i < pixels->size(); i++)
        html += std::to_string((*pixels)[i]) + " ";*/

    std::string html = "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n";
    html.append("<link href=\"https://fonts.googleapis.com/css?family=Inconsolata\" rel=\"stylesheet\"/>\n");
    html.append("<link href=\"https://fonts.googleapis.com/css?family=Material+Icons\" rel=\"stylesheet\"/>\n");
    html.append("<script src=\"reconnecting-websocket.js?" VERSION_STRING "\" defer></script>\n");
    html.append("<script src=\"//cdnjs.cloudflare.com/ajax/libs/numeral.js/2.0.6/numeral.min.js\"></script>\n");
    html.append("<script src=\"https://cdn.jsdelivr.net/gh/jvo203/fits_web_ql/htdocs/fitswebql/ra_dec_conversion.js\"></script>\n");

    //bootstrap
    html.append("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, user-scalable=no, minimum-scale=1, maximum-scale=1\">\n");
    html.append("<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">\n");
    html.append("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js\"></script>\n");
    html.append("<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>\n");

    //GAIAWebQL main JavaScript + CSS
    html.append("<script src=\"gaiawebql.js?" VERSION_STRING "\" defer></script>\n");
    //html.append("<link rel=\"stylesheet\" href=\"gaiawebql.css?" VERSION_STRING "\"/>\n");

#ifdef PRODUCTION
    html.append(R"(<script>var WS_SOCKET = 'wss://';</script>)");
#else
    html.append(R"(<script>var WS_SOCKET = 'ws://';</script>)");
#endif

    //HTML content
    html.append("<title>GAIA DR2 WebQL</title></head><body>\n");
    html.append("<div id='session-data' style='width: 0; height: 0;' data-search-count='" + std::to_string(db_index.size()) + "' ");
    html.append("data-uuid='" + uuid + "' ");
    html.append("data-where='" + where + "' ");
    html.append("data-server-version='" + std::string(VERSION_STRING) + "' data-server-string='" + SERVER_STRING + "'></div>\n");

    //add the main HTML container element
    html.append("<div id="
                "main"
                " class="
                "container"
                ">\n");
    //html.append("<h1>GAIA DR2 WebQL</h1>");
    html.append("<h3 id=\"coords\"></h3>");
    html.append("<h3 id=\"processing\">processing request; #HEALPix search pixels: " + std::to_string(db_index.size()) + "</h3>");
    html.append("<h3 id="
                "completed"
                "></h3>");
    html.append("<div class=\"progress\">");
    html.append("<div id=\"progress-bar\" class=\"progress-bar progress-bar-info progress-bar-striped\" role=\"progressbar\" aria-valuenow=0 aria-valuemin=0 aria-valuemax=" + std::to_string(db_index.size()) + " style=\"width:0%\">0/0</div></div>");
    html.append("</div");

    //html.append(R"(<script>main();</script>)");

    html.append("</body></html>");

    size_t size = html.length();
    write_status(res, 200, "OK");
    write_content_length(res, size);
    write_content_type(res, "text/html");
    res->write("\r\n", 2);
    res->write((const char *)html.data(), size);
    res->write("\r\n\r\n", 4);
}

void execute_gaia_cone_search(uWS::HttpResponse *res, double ra, double dec, double radius, std::string where)
{
    int order = 4;
    T_Healpix_Base<int> hpx(order, RING);

    double theta, phi, rad;

    theta = halfpi - dec * HTMfunc_Pr;
    phi = ra * HTMfunc_Pr;
    rad = radius * HTMfunc_Pr;
    pointing ang(theta, phi);

    int idx = hpx.ang2pix(ang);
    auto pixels = std::make_shared<std::vector<int>>(std::move(hpx.query_disc_inclusive(ang, rad).toVector()));

    std::string uuid = []() -> std::string {
        uuid_t binuuid;
        uuid_generate(binuuid);
        char *uuid = (char *)malloc(37);

        if (uuid != NULL)
        {
            uuid_unparse(binuuid, uuid);
            std::string uuid_str(uuid);
            free(uuid);
            return uuid_str;
        }
        else
            return std::string("NULL");
    }();

    {
        std::lock_guard req_lock(requests_mtx);
        requests.insert(std::make_pair(uuid, std::make_shared<struct db_search_job>()));
    }

    {
        std::lock_guard res_lock(results_mtx);
        results.insert(std::make_pair(uuid, std::make_shared<struct gaia_plots>()));
    }

    std::thread([=]() {
        printf("starting a GAIA database search thread, job id:%s\n", uuid.c_str());

        double longitude = ra * HTMfunc_Pr;
        double latitude = dec * HTMfunc_Pr;

        int max_threads = omp_get_max_threads();
        //std::vector<std::shared_ptr<struct gaia_hist>> thread_hist; //(max_threads);
        std::vector<OmniCoords> thread_coords(max_threads);

        struct gaia_hist global_hist;
        char name[255];

        sprintf(name, "%s/HR", uuid.c_str());
        global_hist.HR = new TH2D(name, "Hertzsprung-Russell diagram", 600, -1.0, 5.0, 600, -5.0, 15.0);
        /*global_hist.hr_hist = new TH2D(uuid.c_str(), "Hertzsprung-Russell diagram", 600, -1.0, 5.0, 600, -1.0, 1.0);
        global_hist.hr_hist->SetCanExtend(TH1::kAllAxes);*/
        //global_hist.hr_hist->SetBit(TH1::kAutoBinPTwo);

        sprintf(name, "%s/XY", uuid.c_str());
        global_hist.XY = new TH2D(name, "X-Y", 600, -0.1, 0.1, 600, -0.1, 0.1);
        global_hist.XY->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZ", uuid.c_str());
        global_hist.RZ = new TH2D(name, "R-Z", 600, -0.1, 0.1, 600, -0.1, 0.1);
        global_hist.RZ->SetCanExtend(TH1::kAllAxes);

        //3D histograms
        sprintf(name, "%s/XYVR", uuid.c_str());
        global_hist.XYVR = new TH3D(name, "X-Y-V_{R}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVR->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/XYVPhi", uuid.c_str());
        global_hist.XYVPhi = new TH3D(name, "X-Y-V_{\\Phi}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVPhi->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/XYVZ", uuid.c_str());
        global_hist.XYVZ = new TH3D(name, "X-Y-V_{Z}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.XYVZ->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVR", uuid.c_str());
        global_hist.RZVR = new TH3D(name, "R-Z-V_{R}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVR->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVPhi", uuid.c_str());
        global_hist.RZVPhi = new TH3D(name, "R-Z-V_{\\Phi}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVPhi->SetCanExtend(TH1::kAllAxes);

        sprintf(name, "%s/RZVZ", uuid.c_str());
        global_hist.RZVZ = new TH3D(name, "R-Z-V_{Z}", 100, -0.1, 0.1, 100, -0.1, 0.1, 100, -0.1, 0.1);
        global_hist.RZVZ->SetCanExtend(TH1::kAllAxes);

        for (int i = 0; i < max_threads; i++)
        {
            thread_coords[i].change_sol_pos(8.3, 0.027);

            char name[255];
            sprintf(name, "%s/%d", uuid.c_str(), (i + 1));

            //struct gaia_hist hist;
            //hist.hr_hist = new TH2D(name, "M_{G} vs. BP-RP", 600, -1.0, 5.0, 600, -5.0, 15.0);
            /*hist->SetCanExtend(TH1::kAllAxes);
            hist->SetBit(TH1::kAutoBinPTwo);*/
            //thread_hist[i] = hist;
            //thread_hist.push_back(std::make_shared<struct gaia_hist>(std::move(hist)));
        }

        boost::lockfree::stack<struct plot_data> plot_stack(100000);
        boost::atomic<bool> search_done(false);

        std::thread plot_thread([&]() {
            std::cout << "starting a histogram thread function.\n";

            unsigned long counter = 0;
            struct plot_data data;

            auto plotter = [&](struct plot_data data) {
                global_hist.HR->Fill(data.bp_rp, data.M_G);
                global_hist.XY->Fill(data.X, data.Y);
                global_hist.RZ->Fill(data.R, data.Z);
                global_hist.XYVR->Fill(data.X, data.Y, data.VR);
                global_hist.XYVPhi->Fill(data.X, data.Y, data.VPhi);
                global_hist.XYVZ->Fill(data.X, data.Y, data.VZ);
                global_hist.RZVR->Fill(data.R, data.Z, data.VR);
                global_hist.RZVPhi->Fill(data.R, data.Z, data.VPhi);
                global_hist.RZVZ->Fill(data.R, data.Z, data.VZ);
            };

            //older Boost on py1 does not have consume_all!!!

            while (!search_done)
            {
                //counter += plot_stack.consume_all(plotter);
                while (plot_stack.pop(data))
                {
                    plotter(data);
                    counter++;
                }
            }

            //counter += plot_stack.consume_all(plotter);
            while (plot_stack.pop(data))
            {
                plotter(data);
                counter++;
            }

            std::cout << "a histogram thread ended after processing " << counter << " values.\n";
        });

        boost::atomic<bool> search_aborted(false);

#pragma omp parallel shared(global_hist)
#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < pixels->size(); i++)
        {
            int tid = omp_get_thread_num();
            bool abort_search = false;
            int hpx = (*pixels)[i];

            printf("cone search for hpx %d\n", hpx);

            try
            {
                auto entry = db_index.at(hpx);
                abort_search = gaia_db_cone_search(hpx, entry, uuid, longitude, latitude, rad, where, thread_coords[tid], plot_stack, global_hist);
            }
            catch (const std::out_of_range &err)
            {
                printf("no entry found for hpx=%d\n", hpx);
            }

            if (abort_search)
            {
                search_aborted = true;
                printf("cancelling a parallel OpenMP search loop.\n");

#pragma omp cancel for
            }
        }

        printf("OpenMP parallel for loop done.\n");

        search_done = true;
        plot_thread.join();

        if (!search_aborted)
        {

            for (int i = 0; i < max_threads; i++)
            {
                //global_hist.hr_hist->Add(thread_hist[i]->hr_hist);
                //delete thread_hist[i]->hr_hist;
            }

            //the H-R diagram
            {
                ReverseYData(global_hist.HR);
                global_hist.HR->SetStats(false);
                global_hist.HR->GetXaxis()->SetTitle("(BP-RP) [mag]");
                global_hist.HR->GetYaxis()->SetTitle("M_{G} [mag]");
                //global_hist->GetZaxis()->SetTitle("star density");

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.HR->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);
                ReverseYAxis(global_hist.HR);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "HR diagram PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (1) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->hr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y plot
            {
                global_hist.XY->SetStats(false);
                global_hist.XY->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XY->GetYaxis()->SetTitle("Y [kpc]");
                SetView2D(global_hist.XY);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XY->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (2) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xy = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z plot
            {
                global_hist.RZ->SetStats(false);
                global_hist.RZ->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZ->GetYaxis()->SetTitle("Z [kpc]");
                SetView2D(global_hist.RZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZ->Draw("COLZ"); //COLZ or CONTZ
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (3) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VR plot
            {
                global_hist.XYVR->SetStats(false);
                global_hist.XYVR->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVR->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVR->GetZaxis()->SetTitle("V_{R} [km/s]");
                SetView3D(global_hist.XYVR);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVR->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VR plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (4) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VPhi plot
            {
                global_hist.XYVPhi->SetStats(false);
                global_hist.XYVPhi->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVPhi->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
                SetView3D(global_hist.XYVPhi);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVPhi->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VPhi plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (5) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvphi = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the X-Y-VZ plot
            {
                global_hist.XYVZ->SetStats(false);
                global_hist.XYVZ->GetXaxis()->SetTitle("X [kpc]");
                global_hist.XYVZ->GetYaxis()->SetTitle("Y [kpc]");
                global_hist.XYVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
                SetView3D(global_hist.XYVZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.XYVZ->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "X-Y-VZ plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (6) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->xyvz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VR plot
            {
                global_hist.RZVR->SetStats(false);
                global_hist.RZVR->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVR->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVR->GetZaxis()->SetTitle("V_{R} [km/s]");
                SetView3D(global_hist.RZVR);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVR->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VR plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (7) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvr = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VPhi plot
            {
                global_hist.RZVPhi->SetStats(false);
                global_hist.RZVPhi->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVPhi->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVPhi->GetZaxis()->SetTitle("V_{\\Phi} [km/s]");
                SetView3D(global_hist.RZVPhi);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVPhi->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VPhi plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (8) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvphi = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }

            //the R-Z-VZ plot
            {
                global_hist.RZVZ->SetStats(false);
                global_hist.RZVZ->GetXaxis()->SetTitle("R [kpc]");
                global_hist.RZVZ->GetYaxis()->SetTitle("Z [kpc]");
                global_hist.RZVZ->GetZaxis()->SetTitle("V_{Z} [km/s]");
                SetView3D(global_hist.RZVZ);

                std::lock_guard lock(root_mtx);
                //gStyle->SetImageScaling(3.);//the HTML canvas image is too big
                TCanvas *c = new TCanvas("", "", 600, 600);
                c->SetBatch(true);
                c->SetGrid(true);
                global_hist.RZVZ->Draw("ISO");
                c->SetRightMargin(0.13);

                TImage *img = TImage::Create();
                img->FromPad(c);

                char *image_buffer;
                int image_size;

                img->GetImageBuffer(&image_buffer, &image_size);
                std::cout << "R-Z-VZ plot PNG graphic size: " << image_size << std::endl;

                //base64 conversion
                char *output = base64((const unsigned char *)image_buffer, image_size);

                if (output != NULL)
                {
                    char *encoded = json_encode_string(output);

                    if (encoded != NULL)
                    {
                        std::lock_guard res_lock(results_mtx);

                        //send json via websockets
                        try
                        {
                            //send completed search jobs via websockets
                            std::ostringstream json;

                            json << "{";
                            json << "\"type\" : \"plot\",";
                            /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                            json << "\"thread\" : " << (9) << ",";
                            //json << "\"total\" : " << (1) << ",";
                            json << "\"density_plot\" : " << encoded;
                            json << "}";

                            std::lock_guard lock(progress_mtx);
                            auto ws = progress_list.at(uuid);

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                            //remove the results
                            results.erase(uuid);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no websocket connection found for a job request %s\n", uuid.c_str());

                            try
                            {
                                auto result = results.at(uuid);
                                result->rzvz = std::string(encoded);
                            }
                            catch (const std::out_of_range &err)
                            {
                                printf("no entry found in results for a job request %s\n", uuid.c_str());
                            }
                        }

                        free(encoded);
                    }

                    free(output);
                }

                free(image_buffer);

                delete img;
                delete c;
            }
        }
        else
        {
            //purge the results
            std::lock_guard res_lock(results_mtx);
            results.erase(uuid);
        }

        delete global_hist.HR;
        delete global_hist.XY;
        delete global_hist.RZ;
        delete global_hist.XYVR;
        delete global_hist.XYVPhi;
        delete global_hist.XYVZ;
        delete global_hist.RZVR;
        delete global_hist.RZVPhi;
        delete global_hist.RZVZ;

        std::lock_guard lock(requests_mtx);
        requests.erase(uuid);
    })
        .detach();

    /*std::string html = "GAIA DR2 WebQL::<ra:" + std::to_string(ra) + ", dec:" + std::to_string(dec) + ", radius:" + std::to_string(radius) + "> --> job id:" + uuid + ", hpx:" + std::to_string(idx) + ", searching pixels: ";

    for (int i = 0; i < pixels->size(); i++)
        html += std::to_string((*pixels)[i]) + " ";*/

    std::string html = "<!DOCTYPE html>\n<html>\n<head>\n<meta charset=\"utf-8\">\n";
    html.append("<link href=\"https://fonts.googleapis.com/css?family=Inconsolata\" rel=\"stylesheet\"/>\n");
    html.append("<link href=\"https://fonts.googleapis.com/css?family=Material+Icons\" rel=\"stylesheet\"/>\n");
    html.append("<script src=\"reconnecting-websocket.js?" VERSION_STRING "\" defer></script>\n");
    html.append("<script src=\"//cdnjs.cloudflare.com/ajax/libs/numeral.js/2.0.6/numeral.min.js\"></script>\n");
    html.append("<script src=\"https://cdn.jsdelivr.net/gh/jvo203/fits_web_ql/htdocs/fitswebql/ra_dec_conversion.js\"></script>\n");

    //bootstrap
    html.append("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, user-scalable=no, minimum-scale=1, maximum-scale=1\">\n");
    html.append("<link rel=\"stylesheet\" href=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css\">\n");
    html.append("<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/3.1.1/jquery.min.js\"></script>\n");
    html.append("<script src=\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js\"></script>\n");

    //GAIAWebQL main JavaScript + CSS
    html.append("<script src=\"gaiawebql.js?" VERSION_STRING "\" defer></script>\n");
    //html.append("<link rel=\"stylesheet\" href=\"gaiawebql.css?" VERSION_STRING "\"/>\n");

#ifdef PRODUCTION
    html.append(R"(<script>var WS_SOCKET = 'wss://';</script>)");
#else
    html.append(R"(<script>var WS_SOCKET = 'ws://';</script>)");
#endif

    //HTML content
    html.append("<title>GAIA DR2 WebQL</title></head><body>\n");
    html.append("<div id='session-data' style='width: 0; height: 0;' data-search-count='" + std::to_string(pixels->size()) + "' ");
    html.append("data-uuid='" + uuid + "' ");
    html.append("data-ra='" + std::to_string(ra) + "' ");
    html.append("data-dec='" + std::to_string(dec) + "' ");
    html.append("data-radius='" + std::to_string(radius) + "' ");
    html.append("data-where='" + where + "' ");
    html.append("data-server-version='" + std::string(VERSION_STRING) + "' data-server-string='" + SERVER_STRING + "'></div>\n");

    //add the main HTML container element
    html.append("<div id="
                "main"
                " class="
                "container"
                ">\n");
    //html.append("<h1>GAIA DR2 WebQL</h1>");
    html.append("<h3 id=\"coords\"></h3>");
    html.append("<h3 id=\"processing\">processing request; #HEALPix search pixels: " + std::to_string(pixels->size()) + "</h3>");
    html.append("<h3 id="
                "completed"
                "></h3>");
    html.append("<div class=\"progress\">");
    html.append("<div id=\"progress-bar\" class=\"progress-bar progress-bar-info progress-bar-striped\" role=\"progressbar\" aria-valuenow=0 aria-valuemin=0 aria-valuemax=" + std::to_string(pixels->size()) + " style=\"width:0%\">0/0</div></div>");
    html.append("</div");

    //html.append(R"(<script>main();</script>)");

    html.append("</body></html>");

    size_t size = html.length();
    write_status(res, 200, "OK");
    write_content_length(res, size);
    write_content_type(res, "text/html");
    res->write("\r\n", 2);
    res->write((const char *)html.data(), size);
    res->write("\r\n\r\n", 4);
}

char *base64(const unsigned char *input, int length)
{
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

void ReverseYAxis(TH1 *h)
{
    // Remove the current axis
    h->GetYaxis()->SetLabelOffset(999);
    h->GetYaxis()->SetTickLength(0);
    // Redraw the new axis
    gPad->Update();
    TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
                                 gPad->GetUymax(),
                                 gPad->GetUxmin() - 0.001,
                                 gPad->GetUymin(),
                                 h->GetYaxis()->GetXmin(),
                                 h->GetYaxis()->GetXmax(),
                                 510, "+");
    newaxis->SetLabelOffset(-0.03);
    newaxis->Draw();
}

void ReverseYData(TH2 *h)
{
    Int_t nx = h->GetNbinsX();
    Int_t ny = h->GetNbinsY();

    for (Int_t i = 0; i < nx; i++)
    {
        for (Int_t j = 0; j < ny / 2; j++)
        {
            Int_t a = h->GetBin(i, j);
            Int_t b = h->GetBin(i, ny - 1 - j);

            auto tmp = h->GetBinContent(a);
            auto tmp2 = h->GetBinContent(b);

            h->SetBinContent(a, tmp2);
            h->SetBinContent(b, tmp);
        }
    }

    h->ComputeIntegral();
}

void SetView2D(TH2 *h)
{
    auto mean_x = h->GetMean(1);
    auto std_x = h->GetStdDev(1);
    auto mean_y = h->GetMean(2);
    auto std_y = h->GetStdDev(2);

    h->GetXaxis()->SetRangeUser(mean_x - 10.0 * std_x, mean_x + 10.0 * std_x);
    h->GetYaxis()->SetRangeUser(mean_y - 10.0 * std_y, mean_y + 10.0 * std_y);
}

void SetView3D(TH3 *h)
{
    auto mean_x = h->GetMean(1);
    auto std_x = h->GetStdDev(1);
    auto mean_y = h->GetMean(2);
    auto std_y = h->GetStdDev(2);
    auto mean_z = h->GetMean(3);
    auto std_z = h->GetStdDev(3);

    h->GetXaxis()->SetRangeUser(mean_x - 10.0 * std_x, mean_x + 10.0 * std_x);
    h->GetYaxis()->SetRangeUser(mean_y - 10.0 * std_y, mean_y + 10.0 * std_y);
    h->GetZaxis()->SetRangeUser(mean_z - 10.0 * std_z, mean_z + 10.0 * std_z);
}

int main(int argc, char *argv[])
{
    curl_global_init(CURL_GLOBAL_ALL);
    ROOT::EnableThreadSafety();

    //load the db healpix index file
    load_db_index("gaiadr2-table.dat");

    // register signal SIGINT and signal handler
    signal(SIGINT, signalHandler);

    //parse local command-line options
    if (argc > 2)
    {
        for (int i = 1; i < argc - 1; i++)
        {
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
    std::transform(threads.begin(), threads.end(), threads.begin(), [](std::thread *t) {
        return new std::thread([]() {
            uWS::Hub h;

            h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
                std::string uri = req.getUrl().toString();

                std::cout << "HTTP request for " << uri << std::endl;

                //root
                if (uri == "/")
                    return serve_file(res, "/index.html");

                //GAIAWebQL entry
                if (uri.find("GAIAWebQL.html") != std::string::npos)
                {
                    //get a position of '?'
                    size_t pos = uri.find("?");

                    if (pos != std::string::npos)
                    {
                        double ra = NAN;
                        double dec = NAN;
                        double radius = NAN;
                        std::string where;

                        struct search_criteria search;
                        bool valid_params = false;

                        //using std::string for now as std::string_view is broken
                        //in the Intel C++ compiler v19 Update 1
                        //(works OK with v19 Update 2)
                        //LLVM CLANG works OK with std::string_view

                        std::string query = uri.substr(pos + 1, std::string::npos);
                        std::cout << "query: (" << query << ")" << std::endl;

                        std::vector<std::string> params;
                        boost::split(params, query, [](char c) { return c == '&'; });

                        for (auto const &s : params)
                        {
                            //find '='
                            size_t pos = s.find("=");

                            if (pos != std::string::npos)
                            {
                                std::string key = s.substr(0, pos);
                                std::string value = s.substr(pos + 1, std::string::npos);

                                if (key == "xmin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.X_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.X_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "xmax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.X_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.X_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "ymin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Y_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Y_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "ymax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Y_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Y_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "zmin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Z_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Z_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "zmax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Z_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Z_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "rmin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.R_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.R_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "rmax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.R_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.R_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "phimin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Phi_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Phi_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "phimax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.Phi_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.Phi_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "mgmin")
                                {
                                    char *e;
                                    errno = 0;

                                    search.M_G_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.M_G_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "mgmax")
                                {
                                    char *e;
                                    errno = 0;

                                    search.M_G_max = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.M_G_max = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "parallax_over_error")
                                {
                                    char *e;
                                    errno = 0;

                                    search.parallax_error_min = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        search.parallax_error_min = NAN;
                                    else
                                        valid_params = true;
                                }

                                if (key == "ra")
                                {
                                    char *e;
                                    errno = 0;

                                    //ra = atof(value.c_str());
                                    ra = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        ra = NAN;
                                }

                                if (key == "dec")
                                {
                                    char *e;
                                    errno = 0;

                                    //dec = atof(value.c_str());
                                    dec = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        dec = NAN;
                                }

                                if (key == "radius")
                                {
                                    char *e;
                                    errno = 0;

                                    //radius = atof(value.c_str());
                                    radius = std::strtod(value.c_str(), &e);
                                    if (*e != '\0' || // error, we didn't consume the entire string
                                        errno != 0)   // error, overflow or underflow
                                        radius = NAN;
                                }

                                if (key == "where")
                                {
                                    CURL *curl = curl_easy_init();

                                    char *str = curl_easy_unescape(curl, value.c_str(), value.length(), NULL);
                                    where = std::string(str);
                                    curl_free(str);

                                    curl_easy_cleanup(curl);
                                }
                            }
                        }

                        std::cout << "ra:" << ra << ", dec:" << dec << ", radius:" << radius;

                        if (where != "")
                            std::cout << " where " << where;
                        std::cout << std::endl;

                        if (valid_params)
                        {
                            print_search_criteria(&search);
                            return execute_gaia(res, &search, where);
                        }

                        if (std::isnan(ra) || std::isnan(dec) || std::isnan(radius))
                        {
                            const std::string not_found = "ERROR: please specify valid parameters ra/dec/radius.";
                            res->end(not_found.data(), not_found.length());
                            return;
                        }
                        else
                            return execute_gaia_cone_search(res, ra, dec, radius, where);
                    }
                    else
                    {
                        const std::string not_found = "ERROR: URL parameters not found.";
                        res->end(not_found.data(), not_found.length());
                        return;
                    }
                }

                return serve_file(res, uri);
            });

            h.onConnection([](uWS::WebSocket<uWS::SERVER> *ws, uWS::HttpRequest req) {
                std::string url = req.getUrl().toString();
                std::cout << "[µWS] onConnection URL(" << url << ")" << std::endl;

                ws->setUserData(NULL);

                std::size_t slash = url.find_last_of("/");
                if (slash != std::string::npos)
                {
                    if (url.size() > slash)
                    {
                        std::string job_id = url.substr(slash + 1, std::string::npos);
                        std::cout << "[µWS] job id: " << job_id << std::endl;

                        ws->setUserData(strdup(job_id.c_str()));

                        try
                        {
                            std::lock_guard res_lock(results_mtx);
                            auto result = results.at(job_id);
                            bool erase = false;

                            if (result->hr != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (1) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->hr.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->xy != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (2) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->xy.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->rz != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (3) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->rz.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->xyvr != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (4) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->xyvr.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->xyvphi != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (5) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->xyvphi.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->xyvz != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (6) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->xyvz.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->rzvr != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (7) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->rzvr.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->rzvphi != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (8) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->rzvphi.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (result->rzvz != "")
                            {
                                //send the result via websockets
                                std::ostringstream json;

                                json << "{";
                                json << "\"type\" : \"plot\",";
                                /*json << "\"thread\" : " << (max_threads + 1) << ",";
                        json << "\"total\" : " << (max_threads + 1) << ",";*/
                                json << "\"thread\" : " << (9) << ",";
                                //json << "\"total\" : " << (1) << ",";
                                json << "\"density_plot\" : " << result->rzvz.c_str();
                                json << "}";

                                ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                                erase = true;
                            }

                            if (erase)
                                results.erase(job_id);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no entry found in results for a job request %s\n", job_id.c_str());
                            ws->terminate(); //force a connection close in a thread-safe manner (the status code will be 1006)
                        }

                        try
                        {
                            std::lock_guard req_lock(requests_mtx);
                            auto entry = requests.at(job_id);

                            {
                                std::lock_guard lock(progress_mtx);
                                progress_list.insert(std::make_pair(job_id, ws));

                                // show content:
                                /*for (auto &x : progress_list)
                                    std::cout << x.first << ": " << x.second << std::endl;*/
                            }

                            //send all completed jobs
                            std::lock_guard lock(entry->completed_mtx);
                            size_t len = entry->completed.size();

                            std::ostringstream json;
                            json << "{ \"completed\" : " << len << ", \"elapsed\" : " << (std::time(nullptr) - entry->timestamp) << " }";

                            ws->send(json.str().c_str(), json.tellp(), uWS::OpCode::TEXT);
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no entry found for a job request %s\n", job_id.c_str());
                            //ws->terminate(); //force a connection close in a thread-safe manner (the status code will be 1006)
                        }
                    }
                }
            });

            h.onDisconnection([](uWS::WebSocket<uWS::SERVER> *ws, int code, char *message, size_t length) {
                std::string msg = std::string(message, length);
                std::cout << "[µWS] onDisconnection " << msg << " code: " << code << std::endl;

                void *job_id = ws->getUserData();

                if (job_id != NULL)
                {
                    //remove itself from the progress_list
                    std::lock_guard lock(progress_mtx);
                    progress_list.erase(std::string((char *)job_id));

                    //client leaving, abort any search under way
                    if (code == 1001)
                    {
                        try
                        {
                            std::lock_guard req_lock(requests_mtx);
                            auto entry = requests.at(std::string((char *)job_id));

                            std::lock_guard lock(entry->completed_mtx);
                            entry->abort = true;
                        }
                        catch (const std::out_of_range &err)
                        {
                            printf("no entry found for a job request %s\n", (char *)job_id);
                        }
                    }

                    //remove the results too
                    /*std::lock_guard res_lock(results_mtx);
                    results.erase(std::string((char *)job_id));*/

                    // show content:
                    /*for (auto &x : progress_list)
                        std::cout << x.first << ": " << x.second << std::endl;*/

                    free(job_id);
                }
            });

            h.onMessage([](uWS::WebSocket<uWS::SERVER> *ws, char *message, size_t length, uWS::OpCode opCode) {
                //ws->send(message, length, opCode);

                std::string msg = std::string(message, length);

                std::cout << "[µWS] " << msg << std::endl;
            });

            // This makes use of the SO_REUSEPORT of the Linux kernel
            // Other solutions include listening to one port per thread
            // with or without some kind of proxy inbetween
            if (!h.listen(server_port, nullptr, uS::ListenOptions::REUSE_PORT))
            {
                std::cout << "Failed to listen\n";
            }

            std::cout << "Launching a uWS::HTTP/WS thread\n";

            h.run();
        });
    });

    std::for_each(threads.begin(), threads.end(), [](std::thread *t) {
        t->join();
    });
}
