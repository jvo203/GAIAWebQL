#pragma once

#include <cfloat>
#include <string.h>
#include <tuple>
#include <vector>

// CERN ROOT
#include <TH2.h>

/*#include <boost/histogram.hpp>
namespace bh = boost::histogram;*/

#define BURN_IN 1000000
#define NBINS 600

/*using histogram_t = decltype(bh::make_histogram(
    bh::axis::regular<double, bh::use_default, bh::use_default,
                      bh::axis::option::growth_t>(),
    bh::axis::regular<double, bh::use_default, bh::use_default,
                      bh::axis::option::growth_t>()));*/

class SeedH2 {
public:
  SeedH2(bool _invert = false);
  ~SeedH2();

public:
  void set_title(std::string _title) { title = _title; };
  void update(float _x, float _y);
  void flush();
  void save(std::string uuid, std::string type);

private:
  std::string title;
  double x_min, x_max;
  double y_min, y_max;
  std::vector<std::tuple<float, float>> data;
  bool init_done;
  bool invert;
  // histogram_t _hist;
  TH2 *_hist;
};