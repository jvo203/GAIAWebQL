#pragma once

#include <tuple>
#include <vector>

//#include <boost/histogram.hpp>

// CERN ROOT
#include <TH2.h>

#define BURN_IN 10000

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

public:
  TH2 *hist;

private:
  std::vector<std::tuple<float, float>> data;
  bool init_done;
};