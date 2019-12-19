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
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH2::fill(float _x, float _y) {
  if (_x < x_min) {
    double d = x_min - _x;
    double x_min_new = x_min - SCALE * d;

    rebin_x(x_min_new, x_max);
  }

  if (_x > x_max) {
    double d = _x - x_max;
    double x_max_new = x_max + SCALE * d;

    rebin_x(x_min, x_max_new);
  }

  if (_y < y_min) {
    double d = y_min - _y;
    double y_min_new = y_min - SCALE * d;

    rebin_y(y_min_new, y_max);
  }

  if (_y > y_max) {
    double d = _y - y_max;
    double y_max_new = y_max + SCALE * d;

    rebin_y(y_min, y_max_new);
  }

  double x = double(_x - x_min) / double(x_max - x_min);
  int idx = std::clamp((int)(x * NBINS), 0, NBINS - 1);

  double y = double(_y - y_min) / double(y_max - y_min);
  int idy = std::clamp((int)(y * NBINS), 0, NBINS - 1);

  bin_data[idy][idx]++;
}

void SeedH2::rebin_x(double x_min_new, double x_max_new) {
  uint64_t bins[NBINS];

  printf("rebin_x::%s x_min: %f ==> %f, x_max: %f ==> %f\n", title.c_str(),
         x_min, x_min_new, x_max, x_max_new);

  double width = (x_min - x_max) / double(NBINS);

  for (int j = 0; j < NBINS; j++) {
    memset(bins, 0, NBINS * sizeof(uint64_t));

    for (int i = 0; i < NBINS; i++) {
      // get the bin centre
      uint64_t count = bin_data[j][i];
      double centre = x_min + i * width + 0.5 * width;

      // update a new bin
      double x = double(centre - x_min_new) / double(x_max_new - x_min_new);
      int idx = std::clamp((int)(x * NBINS), 0, NBINS - 1);
      bins[idx] += count;
    }

    // copy new bins to bin_data
    memcpy(bin_data[j], bins, NBINS * sizeof(uint64_t));
  }

  // update the range
  x_min = x_min_new;
  x_max = x_max_new;
}

void SeedH2::rebin_y(double y_min_new, double y_max_new) {
  uint64_t bins[NBINS];

  printf("rebin_y::%s y_min: %f ==> %f, y_max: %f ==> %f\n", title.c_str(),
         y_min, y_min_new, y_max, y_max_new);

  double height = (y_min - y_max) / double(NBINS);

  for (int j = 0; j < NBINS; j++) {
    memset(bins, 0, NBINS * sizeof(uint64_t));

    for (int i = 0; i < NBINS; i++) {
      // get the bin centre
      uint64_t count = bin_data[i][j];
      double centre = y_min + i * height + 0.5 * height;

      // update a new bin
      double y = double(centre - y_min_new) / double(y_max_new - y_min_new);
      int idy = std::clamp((int)(y * NBINS), 0, NBINS - 1);
      bins[idy] += count;
    }

    // copy new bins to bin_data
    for (int i = 0; i < NBINS; i++)
      bin_data[i][j] = bins[i];
  }

  // update the range
  y_min = y_min_new;
  y_max = y_max_new;
};

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

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    fill(_x, _y);
  }

  // allocate a new Boost.Histogram
  /*_hist = make_histogram(axis::regular<float>(600, x_min, x_max,
     "_x"), axis::regular<float>(600, y_min, y_max, "_y"));*/

  data.clear();

  init_done = true;
}