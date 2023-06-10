#include "QnDataFrame.hpp"

vector <vector<string>> Q1Q1=
{
  {"F1_PLAIN", "F2_PLAIN"},
  {"F1_PLAIN", "F3_PLAIN"},
  {"F2_PLAIN", "F3_PLAIN"},

  {"Tp_PLAIN", "F1_PLAIN"},
  {"Tp_PLAIN", "F2_PLAIN"},
  {"Tp_PLAIN", "F3_PLAIN"},

  {"Tpi_PLAIN", "F1_PLAIN"},
  {"Tpi_PLAIN", "F2_PLAIN"},
  {"Tpi_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> u1Q1=
{
  {"proton_PLAIN", "F1_PLAIN"},
  {"proton_PLAIN", "F2_PLAIN"},
  {"proton_PLAIN", "F3_PLAIN"},
  {"proton_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u1=
{
  {"rnd_proton_PLAIN"},
  {"rnd_sub_PLAIN"},
};

void mcpico_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"bimp", 14, 0, 14});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:Q1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
  }
  for (auto &corr:u1Q1)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr : u1 )
  {
    std::array<std::string, 1> qn{ corr.at(0) };
    string corrName=corr.at(0);
    auto corr_function = []( const Qn::QVector &a ){ return a.x(1); };
    auto sumW = [](const Qn::QVector &a) { return a.sumweights(); };
    corrBuilder.AddCorrelationWithInternalReader(corrName, corr_function, sumW, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
}
