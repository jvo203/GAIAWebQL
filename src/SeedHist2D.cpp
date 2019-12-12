#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "SeedHist2D.hpp"

void SeedH2::update(float _x, float _y) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y));
  else {
    // append to the histogram
    this->hist->Fill(_x, _y);
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH2::flush() {
  if (data.size() == 0)
    return;

  double meanX = 0.0;
  double meanY = 0.0;

  // gather the statistics
  for (auto &x : data) {
    meanX += std::get<0>(x);
    meanY += std::get<1>(x);
  }

  meanX /= double(data.size());
  meanY /= double(data.size());

  double stdX = 0.0;
  double stdY = 0.0;

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    stdX += (_x - meanX) * (_x - meanX);
    stdY += (_y - meanY) * (_y - meanY);
  }

  stdX = sqrt(stdX / double(data.size()));
  stdY = sqrt(stdY / double(data.size()));

  double x_min = meanX - 3.0 * stdX;
  double x_max = meanX + 3.0 * stdX;

  double y_min = meanY - 3.0 * stdY;
  double y_max = meanY + 3.0 * stdY;

  // allocate a new ROOT histogram
  boost::uuids::random_generator gen;
  boost::uuids::uuid id = gen();
  const std::string name = boost::uuids::to_string(id);

  this->hist = new TH2D((const char *)name.c_str(), "_x <=> _y", 600, x_min,
                        x_max, 600, y_min, y_max);

  if (this->hist == NULL) {
    printf("[error] TH2D histogram object cannot be created\n");
    return;
  }

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    hist->Fill(_x, _y);
  }

  data.clear();

  init_done = true;
}