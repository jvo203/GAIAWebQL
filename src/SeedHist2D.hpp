#pragma once

#include <tuple>
#include <vector>

#include <boost/histogram.hpp>
using namespace boost::histogram;

// CERN ROOT
#include <TH2.h>

#define BURN_IN 1000000
#define NO_BINS 600

class SeedH2 {
public:
  SeedH2() {
    title = "_x <=> _y";
    init_done = false;
    hist = NULL;
    data.reserve(BURN_IN);
  }
  ~SeedH2() {
    if (hist != NULL)
      delete hist;
  };

public:
  void set_title(std::string _title) { title = _title; };
  void update(float _x, float _y);
  void flush();

public:
  TH2 *hist;
  // histogram _hist;

private:
  std::string title;
  std::vector<std::tuple<float, float>> data;
  bool init_done;

  // a custom histogram
  float x_axis[NO_BINS];
  float y_axis[NO_BINS];
  float bin_data[NO_BINS][NO_BINS];
};