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

#include <TCanvas.h>
#include <TFile.h>
#include <TThread.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

SeedH3::SeedH3(bool _invert) {
  title = "_x <=> _y <=> _z";
  _hist = NULL;
  init_done = false;
  has_data = false;
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
  has_data = true;

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

  _hist = new TProfile2D(name.c_str(), title.c_str(), NBINS3, x_min, x_max,
                         NBINS3, y_min, y_max, z_min, z_max);

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
  if (!has_data)
    return;

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

  filename = tmp + "/" + type + "_mean.root";

  TThread::Lock();

  TFile outputFile(filename.c_str(), "RECREATE");
  outputFile.SetCompressionLevel(
      ROOT::RCompressionSetting::ELevel::kDefaultZLIB);

  TCanvas *c = new TCanvas((name + "_mean").c_str(), title.c_str(), 1000, 600);
  c->SetGrid(true);
  _hist->Draw("CONTZ"); // COLZ or CONTZ
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);

  c->Write();
  outputFile.Close();
  delete c;

  // export the errors too
  TH2D *error = _hist->ProjectionXY();
  error->SetStats(false);
  error->GetXaxis()->SetTitle(x_title.c_str());
  error->GetYaxis()->SetTitle(y_title.c_str());
  error->GetZaxis()->SetTitle(
      (std::string("#sigma_{") + z_value + std::string("} ") + z_unit).c_str());

  // TO DO: fix the main title too
  std::string error_title = title;
  size_t pos = title.find_last_of("-");

  if (pos != std::string::npos) {
    std::string prefix = title.substr(0, pos + 1);
    error_title =
        prefix + std::string("#sigma_{") + z_value + std::string("} ");
    error->SetTitle(error_title.c_str());
  }

  FillRMS(error);

  filename = tmp + "/" + type + "_error.root";

  TFile errorFile(filename.c_str(), "RECREATE");
  errorFile.SetCompressionLevel(
      ROOT::RCompressionSetting::ELevel::kDefaultZLIB);

  c = new TCanvas((name + "_error").c_str(), error_title.c_str(), 1000, 600);
  c->SetGrid(true);
  error->Draw("CONTZ"); // COLZ or CONTZ
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.15);

  c->Write();
  errorFile.Close();

  delete c;
  delete error;

  TThread::UnLock();
}

void SeedH3::FillRMS(TH2D *error) {
  Int_t nx = _hist->GetNbinsX();
  Int_t ny = _hist->GetNbinsY();

  for (Int_t i = 0; i < nx; i++) {
    for (Int_t j = 0; j < ny; j++) {
      Int_t bin = _hist->GetBin(i, j);
      auto value = _hist->GetBinError(bin);
      error->SetBinContent(bin, value);
    }
  }
}