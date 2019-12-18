#pragma once

#include <cfloat>
#include <tuple>
#include <vector>

/*#include <boost/histogram.hpp>
using namespace boost::histogram;*/

#define BURN_IN 1000000
#define NBINS 100
// cannot deal with 600*600, need to allocate the bin data dynamically

class SeedH2 {
public:
  SeedH2() {
    title = "_x <=> _y";
    init_done = false;
    data.reserve(BURN_IN);
    x_min = -FLT_MAX;
    x_max = FLT_MAX;
    y_min = -FLT_MAX;
    y_max = FLT_MAX;
    printf("created a SeedH2::%s\n", title.c_str());
  }
  ~SeedH2() {
    // print a custom histogram
    for (int i = 0; i < NBINS; i++)
      for (int j = 0; j < NBINS; j++)
        printf("%zu", bin_data[i][j]);
    printf("\n");
  };

public:
  void set_title(std::string _title) { title = _title; };
  void update(float _x, float _y);
  void flush();

private:
  void fill(float _x, float _y);

private:
  std::string title;
  double x_min, x_max;
  double y_min, y_max;
  std::vector<std::tuple<float, float>> data;
  bool init_done;

  // a custom histogram
  /*float x_axis[NO_BINS + 1];
  float y_axis[NO_BINS + 1];*/
  uint64_t bin_data[NBINS][NBINS];
};