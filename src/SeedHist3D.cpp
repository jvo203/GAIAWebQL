#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <sys/stat.h>

#include "SeedHist3D.hpp"

#include <TFile.h>
#include <TH2.h>
#include <TThread.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

SeedH3::SeedH3(bool _invert) {
  title = "_x <=> _y <=> _z";
  _hist = NULL;
  init_done = false;
  invert = _invert;
  data.reserve(BURN_IN);
  x_min = -FLT_MAX;
  x_max = FLT_MAX;
  y_min = -FLT_MAX;
  y_max = FLT_MAX;
  z_min = -FLT_MAX;
  z_max = FLT_MAX;
  printf("created a SeedH3::%s\n", title.c_str());
}

SeedH3::~SeedH3() {
  printf("%s x_min: %f x_max: %f; y_min: %f, y_max: %f, z_min: %f, z_max: %f\n",
         title.c_str(), x_min, x_max, y_min, y_max, z_min, z_max);

  if (_hist != NULL)
    delete _hist;
};

void SeedH3::update(float _x, float _y, float _z) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y, _z));
  else {
    // update a histogram
    _hist->Fill(_x, _y, _z);
  }

  if (data.size() == BURN_IN)
    this->flush();
}

void SeedH3::flush() {
  if (data.size() == 0)
    return;

  double meanX = 0.0;
  double meanY = 0.0;
  double meanZ = 0.0;

  // gather the statistics
  for (auto &x : data) {
    meanX += std::get<0>(x);
    meanY += std::get<1>(x);
    meanZ += std::get<2>(x);
  }

  meanX /= double(data.size());
  meanY /= double(data.size());
  meanZ /= double(data.size());

  double stdX = 0.0;
  double stdY = 0.0;
  double stdZ = 0.0;

  for (auto &x : data) {
    double _x, _y, _z;
    std::tie(_x, _y, _z) = x;

    stdX += (_x - meanX) * (_x - meanX);
    stdY += (_y - meanY) * (_y - meanY);
    stdZ += (_z - meanZ) * (_z - meanZ);
  }

  stdX = sqrt(stdX / double(data.size()));
  stdY = sqrt(stdY / double(data.size()));
  stdZ = sqrt(stdZ / double(data.size()));

  x_min = meanX - 3.0 * stdX;
  x_max = meanX + 3.0 * stdX;

  y_min = meanY - 3.0 * stdY;
  y_max = meanY + 3.0 * stdY;

  z_min = meanZ - 3.0 * stdZ;
  z_max = meanZ + 3.0 * stdZ;

  printf("[%s] x_min: %f x_max: %f y_min: %f y_max: %f z_min: %f z_max: %f\n",
         title.c_str(), x_min, x_max, y_min, y_max, z_min, z_max);

  _hist = new TH3D(name.c_str(), title.c_str(), NBINS3, x_min, x_max, NBINS3,
                   y_min, y_max, NBINS3, z_min, z_max);
  _hist->SetCanExtend(TH1::kAllAxes);
  _hist->SetStats(false);
  _hist->GetXaxis()->SetTitle(x_title.c_str());
  _hist->GetYaxis()->SetTitle(y_title.c_str());
  _hist->GetZaxis()->SetTitle(z_title.c_str());

  for (auto &x : data) {
    double _x, _y, _z;
    std::tie(_x, _y, _z) = x;

    _hist->Fill(_x, _y, _z);
  }

  data.clear();
  init_done = true;
}

void SeedH3::export_root(std::string uuid, std::string docs_root,
                         std::string type) {
  std::string filename =
      docs_root + "/gaiawebql/DATA/" + uuid + "/" + type + ".root";

  std::cout << "saving " << title << " into " << filename << std::endl;

  // mkdir DATA/<uuid>.tmp
  std::string tmp = docs_root + "/gaiawebql/DATA/" + uuid + ".tmp";

  if (mkdir(tmp.c_str(), 0777) != 0) {
    // return only in case of errors other than "directory already exists"
    if (errno != EEXIST) {
      perror(title.c_str());
      return;
    }
  }

  TH2D *hist = (TH2D *)_hist->Project3D("xy");

  filename = tmp + "/" + type + ".root";

  TThread::Lock();

  TFile outputFile(filename.c_str(), "RECREATE");
  outputFile.SetCompressionLevel(
      ROOT::RCompressionSetting::ELevel::kDefaultZLIB);
  hist->Write();
  outputFile.Close();

  TThread::UnLock();

  delete hist;
}