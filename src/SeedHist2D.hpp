#pragma once

#include <cfloat>
#include <string.h>
#include <tuple>
#include <vector>

// CERN ROOT
#include <TH2.h>

#define BURN_IN 1000000
#define NBINS 600

class SeedH2 {
public:
  SeedH2(bool _invert = false);
  ~SeedH2();

public:
  void set_title(std::string _name, std::string _title, std::string _x_axis,
                 std::string _y_axis) {
    name = _name;
    title = _title;
    x_title = _x_axis;
    y_title = _y_axis;
  };
  void update(float _x, float _y);
  void flush();
  void save(std::string uuid, std::string type);
  void export_root(std::string uuid, std::string type);

private:
  void ReverseYAxis(TH1 *h);
  void ReverseYData(TH2 *h);

private:
  std::string name, title, x_title, y_title;
  double x_min, x_max;
  double y_min, y_max;
  std::vector<std::tuple<float, float>> data;
  bool init_done;
  bool invert;
  TH2 *_hist;
};