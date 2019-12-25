#pragma once

#include <cfloat>
#include <string.h>
#include <tuple>
#include <vector>

#include <boost/histogram.hpp>
namespace bh = boost::histogram;

#define BURN_IN 1000000
#define NBINS 600
#define SCALE 1.67

using histogram_t =
    decltype(bh::make_histogram(bh::axis::regular<>(), bh::axis::regular<>()));

class SeedH2 {
public:
  SeedH2();
  ~SeedH2();

public:
  void set_title(std::string _title) { title = _title; };
  void update(float _x, float _y);
  void flush();
  void save(std::string uuid, std::string type);

private:
  void fill(float _x, float _y);
  void rebin_x(double x_min_new, double x_max_new);
  void rebin_y(double y_min_new, double y_max_new);

private:
  std::string title;
  double x_min, x_max;
  double y_min, y_max;
  std::vector<std::tuple<float, float>> data;
  bool init_done;

  histogram_t hist_;

  // a custom histogram
  /*float x_axis[NO_BINS + 1];
  float y_axis[NO_BINS + 1];*/
  // uint64_t bin_data[NBINS][NBINS];
  uint64_t **bin_data;
};