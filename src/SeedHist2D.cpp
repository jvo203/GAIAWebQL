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

#include "SeedHist2D.hpp"

#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TThread.h>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

SeedH2::SeedH2(bool _invert) {
  title = "_x <=> _y";
  _hist = NULL;
  init_done = false;
  invert = _invert;
  data.reserve(BURN_IN);
  x_min = -FLT_MAX;
  x_max = FLT_MAX;
  y_min = -FLT_MAX;
  y_max = FLT_MAX;
  printf("created a SeedH2::%s\n", title.c_str());
}
SeedH2::~SeedH2() {
  printf("%s x_min: %f x_max: %f; y_min: %f, y_max: %f\n", title.c_str(), x_min,
         x_max, y_min, y_max);

  if (_hist != NULL)
    delete _hist;
};

void SeedH2::update(float _x, float _y) {
  if (!init_done)
    data.push_back(std::make_tuple(_x, _y));
  else {
    // update a histogram
    _hist->Fill(_x, _y);
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

  x_min = meanX - 3.0 * stdX;
  x_max = meanX + 3.0 * stdX;

  y_min = meanY - 3.0 * stdY;
  y_max = meanY + 3.0 * stdY;

  printf("[%s] x_min: %f x_max: %f y_min: %f y_max: %f\n", title.c_str(), x_min,
         x_max, y_min, y_max);

  /*boost::uuids::uuid uuid = boost::uuids::random_generator()();
  std::string name = boost::lexical_cast<std::string>(uuid);*/

  _hist = new TH2D(name.c_str(), title.c_str(), NBINS, x_min, x_max, NBINS,
                   y_min, y_max);
  _hist->SetCanExtend(TH1::kAllAxes);
  _hist->SetStats(false);
  _hist->GetXaxis()->SetTitle(x_title.c_str());
  _hist->GetYaxis()->SetTitle(y_title.c_str());

  for (auto &x : data) {
    double _x, _y;
    std::tie(_x, _y) = x;

    _hist->Fill(_x, _y);
  }

  data.clear();
  init_done = true;
}

void SeedH2::ReverseYAxis(TH1 *h) {
  // Remove the current axis
  h->GetYaxis()->SetLabelOffset(999);
  h->GetYaxis()->SetTickLength(0);
  // Redraw the new axis
  gPad->Update();
  TGaxis *newaxis =
      new TGaxis(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmin() - 0.001,
                 gPad->GetUymin(), h->GetYaxis()->GetXmin(),
                 h->GetYaxis()->GetXmax(), 510, "+");
  newaxis->SetLabelOffset(-0.03);
  newaxis->Draw();
}

void SeedH2::ReverseYData(TH2 *h) {
  Int_t nx = h->GetNbinsX();
  Int_t ny = h->GetNbinsY();

  for (Int_t i = 0; i < nx; i++) {
    for (Int_t j = 0; j < ny / 2; j++) {
      Int_t a = h->GetBin(i, j);
      Int_t b = h->GetBin(i, ny - 1 - j);

      auto tmp = h->GetBinContent(a);
      auto tmp2 = h->GetBinContent(b);

      h->SetBinContent(a, tmp2);
      h->SetBinContent(b, tmp);
    }
  }

  h->ComputeIntegral();
}

void SeedH2::export_root(std::string uuid, std::string docs_root,
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

  filename = tmp + "/" + type + ".root";

  TThread::Lock();

  TFile outputFile(filename.c_str(), "RECREATE");
  outputFile.SetCompressionLevel(
      ROOT::RCompressionSetting::ELevel::kDefaultZLIB);
  _hist->Write();
  outputFile.Close();

  TThread::UnLock();
}

void SeedH2::save(std::string uuid, std::string docs_root, std::string type) {
  std::string filename =
      docs_root + "/gaiawebql/DATA/" + uuid + "/" + type + ".json";

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

  TAxis *axis;

  axis = _hist->GetXaxis();
  x_min = axis->GetBinCenter(0);
  x_max = axis->GetBinCenter(NBINS - 1);

  axis = _hist->GetYaxis();
  y_min = axis->GetBinCenter(0);
  y_max = axis->GetBinCenter(NBINS - 1);

  // write the bin data
  {
    filename = tmp + "/" + type + ".json";
    std::ofstream json(filename);

    json << "{\n";
    json << "\t\"NBINS\" : " << NBINS << ",\n";
    json << "\t\"xmin\" : " << x_min << ",\n";
    json << "\t\"xmax\" : " << x_max << ",\n";
    json << "\t\"ymin\" : " << y_min << ",\n";
    json << "\t\"ymax\" : " << y_max << ",\n";

    json << "\t\"bins\" : [";

    for (int j = 0; j < NBINS; j++) {
      json << "[";
      for (int i = 0; i < NBINS; i++) {
        Int_t a = _hist->GetBin(i, j);
        auto value = _hist->GetBinContent(a);
        json << (value > 0.0) ? std::to_string(value) : "null";

        if (i != NBINS - 1)
          json << ",";
      }
      json << "]";

      if (j != NBINS - 1)
        json << ",";
    }

    json << "]}\n";
    json.close();

    // gz-compress with Boost
    std::ifstream inStream(filename, std::ios_base::in);
    std::ofstream outStream(filename + ".gz", std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    in.push(boost::iostreams::gzip_compressor());
    in.push(inStream);
    boost::iostreams::copy(in, outStream);
  }
};