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

  bin_data[idx][idy]++;
}

void SeedH2::rebin_x(double x_min_new, double x_max_new) {
  uint64_t bins[NBINS];

  memset(bins, 0, NBINS * sizeof(uint64_t));
}

void SeedH2::rebin_y(double y_min_new, double y_max_new){};

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