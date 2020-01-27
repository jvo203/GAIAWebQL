#pragma once

#include <cfloat>
#include <string.h>
#include <tuple>
#include <vector>

// CERN ROOT
#include <TProfile2D.h>

#define BURN_IN 1000000
#define NBINS3 100

class SeedH3 {
public:
  SeedH3(bool _invert = false);
  ~SeedH3();

public:
  void set_title(std::string _name, std::string _title, std::string _x_axis,
                 std::string _y_axis, std::string _z_axis) {
    name = _name;
    title = _title;
    x_title = _x_axis;
    y_title = _y_axis;
    z_title = _z_axis;
  };
  void update(float _x, float _y, float _z);
  void flush();
  void export_root(std::string uuid, std::string docs_root, std::string type);

private:
  std::string name, title, x_title, y_title, z_title;
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  std::vector<std::tuple<float, float, float>> data;
  bool init_done;
  bool invert;
  TProfile2D *_hist;
};