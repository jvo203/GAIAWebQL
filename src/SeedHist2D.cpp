#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "SeedHist2D.hpp"

void SeedH2::update(float _x, float _y) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y));
  else {
    // a custom histogram
    this->fill(_x, _y);

    // append to the histogram
    this->hist->Fill(_x, _y);
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH2::fill(float _x, float _y) {
  if (_x < x_min || _x >= x_max || _y < y_min || _y >= y_max) {
    // re-bin the histogram
  }

  double x = double(_x - x_min) / double(x_max - x_min);
  int idx = std::clamp((int)(x * NO_BINS), 0, NO_BINS - 1);

  double y = double(_y - y_min) / double(y_max - y_min);
  int idy = std::clamp((int)(y * NO_BINS), 0, NO_BINS - 1);

  bin_data[idx][idy]++;
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

  for (int i = 0; i < NO_BINS; i++)
    for (int j = 0; j < NO_BINS; j++)
      bin_data[i][j] = 0;

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    fill(_x, _y);
  }

  // allocate a new Boost.Histogram
  /*_hist = make_histogram(axis::regular<float>(600, x_min, x_max,
     "_x"), axis::regular<float>(600, y_min, y_max, "_y"));*/

  // allocate a new ROOT histogram
  boost::uuids::random_generator gen;
  boost::uuids::uuid id = gen();
  const std::string name = boost::uuids::to_string(id);

  this->hist = new TH2D((const char *)name.c_str(), this->title.c_str(), 600,
                        x_min, x_max, 600, y_min, y_max);
  this->hist->SetCanExtend(TH1::kAllAxes);

  if (this->hist == NULL) {
    printf("[error] TH2D histogram object cannot be created\n");
    return;
  }

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    hist->Fill(_x, _y);
  }

  data.clear();

  init_done = true;
}