#pragma once

#include <tuple>
#include <vector>

//#include <boost/histogram.hpp>

// CERN ROOT
#include <TH2.h>

#define BURN_IN 100000

class SeedH2 {
public:
  SeedH2() {
    init_done = false;
    hist = NULL;
    data.reserve(BURN_IN);
  }
  ~SeedH2() {
    if (hist != NULL)
      delete hist;
  };

public:
  void update(float _x, float _y);
  void flush();

private:
  std::vector<std::tuple<float, float>> data;
  TH2 *hist;
  bool init_done;
};