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
    hist = NULL;
    data.reserve(BURN_IN);
  }
  ~SeedH2(){};

public:
  void update(float _x, float _y);

private:
  std::vector<std::tuple<float, float>> data;
  TH2 *hist;
};