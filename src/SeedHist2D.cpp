#include "SeedHist2D.hpp"

void SeedH2::update(float _x, float _y) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y));
  else {
    // append to the histogram
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH2::flush() { init_done = true; }