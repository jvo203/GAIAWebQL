#pragma once

#include <cfloat>
#include <string.h>
#include <tuple>
#include <vector>

/*#include <boost/histogram.hpp>
using namespace boost::histogram;*/

#define BURN_IN 1000000
#define NBINS 600
#define SCALE 1.67

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

    // allocate the bin memory
    bin_data = (uint64_t **)malloc(sizeof(uint64_t *) * NBINS);

    if (bin_data != NULL)
      for (int i = 0; i < NBINS; i++) {
        bin_data[i] = (uint64_t *)malloc(sizeof(uint64_t) * NBINS);

        if (bin_data[i] != NULL)
          memset(bin_data[i], 0, NBINS * sizeof(uint64_t));
        else
          printf("error allocating bin_data[%d]\n", i);
      }

    printf("created a SeedH2::%s\n", title.c_str());
  }
  ~SeedH2() {
    // print a custom histogram
    for (int i = 0; i < NBINS; i++)
      for (int j = 0; j < NBINS; j++)
        printf("%zu ", bin_data[i][j]);
    printf("\n");

    // deallocate the bins
    if (bin_data != NULL)
      for (int i = 0; i < NBINS; i++)
        if (bin_data[i] != NULL)
          free(bin_data[i]);
  };

public:
  void set_title(std::string _title) { title = _title; };
  void update(float _x, float _y);
  void flush();

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

  // a custom histogram
  /*float x_axis[NO_BINS + 1];
  float y_axis[NO_BINS + 1];*/
  // uint64_t bin_data[NBINS][NBINS];
  uint64_t **bin_data;
};