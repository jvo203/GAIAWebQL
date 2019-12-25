#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>

#include "SeedHist2D.hpp"

SeedH2::SeedH2(bool _invert) {
  title = "_x <=> _y";
  invert = _invert;
  init_done = false;
  data.reserve(BURN_IN);
  x_min = -FLT_MAX;
  x_max = FLT_MAX;
  y_min = -FLT_MAX;
  y_max = FLT_MAX;
  printf("created a SeedH2::%s\n", title.c_str());
}
SeedH2::~SeedH2() {
  printf("%s x_min: %f x_max: %f; y_min: %f, y_max: %f\n", title.c_str(), x_min,
         x_max, y_min, y_max);
};

void SeedH2::update(float _x, float _y) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y));
  else {
    // update a Boost.Histogram
    _hist(_x, _y);
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH2::flush() {
  if (data.size() == 0)
    return;

  double meanX = 0.0;
  double meanY = 0.0;

  // gather the statistics
  for (auto &x : data) {
    meanX += std::get<0>(x);
    meanY += std::get<1>(x);
  }

  meanX /= double(data.size());
  meanY /= double(data.size());

  double stdX = 0.0;
  double stdY = 0.0;

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    stdX += (_x - meanX) * (_x - meanX);
    stdY += (_y - meanY) * (_y - meanY);
  }

  stdX = sqrt(stdX / double(data.size()));
  stdY = sqrt(stdY / double(data.size()));

  x_min = meanX - 3.0 * stdX;
  x_max = meanX + 3.0 * stdX;

  y_min = meanY - 3.0 * stdY;
  y_max = meanY + 3.0 * stdY;

  printf("[%s] x_min: %f x_max: %f y_min: %f y_max: %f\n", title.c_str(), x_min,
         x_max, y_min, y_max);

  // make a custom histogram
  /*double dx = (x_max - x_min) / double(NO_BINS);
  double dy = (y_max - y_min) / double(NO_BINS);

  for (int i = 0; i < NO_BINS + 1; i++) {
    x_axis[i] = x_min + double(i) * dx;
    y_axis[i] = y_min + double(i) * dy;
  }*/

  // allocate a new Boost.Histogram
  _hist = bh::make_histogram(
      bh::axis::regular<double, bh::use_default, bh::use_default,
                        bh::axis::option::growth_t>(NBINS, x_min, x_max, "_x"),
      bh::axis::regular<double, bh::use_default, bh::use_default,
                        bh::axis::option::growth_t>(NBINS, y_min, y_max, "_y"));

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    _hist(_x, _y);
  }

  data.clear();
  init_done = true;
}

void SeedH2::save(std::string uuid, std::string type) {
  std::string filename = "DATA/" + uuid + "/" + type + ".json";

  std::cout << "saving " << title << " into " << filename << std::endl;

  // mkdir DATA/<uuid>.tmp
  std::string tmp = "DATA/" + uuid + ".tmp";

  if (mkdir(tmp.c_str(), 0777) != 0) {
    // return only in case of errors other than "directory already exists"
    if (errno != EEXIST) {
      perror(title.c_str());
      return;
    }
  }

  x_min = _hist.axis(0).value(0);
  x_max = _hist.axis(0).value(NBINS - 1);

  if (!invert) {
    y_min = _hist.axis(1).value(0);
    y_max = _hist.axis(1).value(NBINS - 1);
  } else {
    y_max = _hist.axis(1).value(0);
    y_min = _hist.axis(1).value(NBINS - 1);
  }

  // write the bin data
  {
    filename = tmp + "/" + type + ".json";
    std::ofstream json(filename);

    json << "{\n";
    json << "\t\"NBINS\" : " << NBINS << ",\n";
    json << "\t\"xmin\" : " << x_min << ",\n";
    json << "\t\"xmax\" : " << x_max << ",\n";
    json << "\t\"ymin\" : " << y_min << ",\n";
    json << "\t\"ymax\" : " << y_max << ",\n";

    json << "\t\"bins\" : [";

    for (int j = 0; j < NBINS; j++) {
      json << "[";
      for (int i = 0; i < NBINS; i++) {
        json << _hist.at(i, invert
                                ? (NBINS - 1 - j)
                                : j); // invert the Y axis for the H-R diagram

        if (i != NBINS - 1)
          json << ",";
      }
      json << "]";

      if (j != NBINS - 1)
        json << ",";
    }

    json << "]}\n";
    json.close();
  }
};