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

vector <vector<string>> Q2Q1Q1{
  {"Tpi2_PLAIN", "F1_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> u1Q1=
{
  {"proton_PLAIN", "F1_PLAIN"},
  {"proton_PLAIN", "F2_PLAIN"},
  {"proton_PLAIN", "F3_PLAIN"},
  {"proton_PLAIN", "psi_rp_PLAIN"},

  {"pi_neg_PLAIN", "F1_PLAIN"},
  {"pi_neg_PLAIN", "F2_PLAIN"},
  {"pi_neg_PLAIN", "F3_PLAIN"},
  {"pi_neg_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q2=
{
  {"proton_PLAIN", "psi_rp_PLAIN"},
  {"pi_neg_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q1Q1=
{
  {"proton_PLAIN", "F1_PLAIN", "F1_PLAIN"},
  {"proton_PLAIN", "F2_PLAIN", "F2_PLAIN"},
  {"proton_PLAIN", "F3_PLAIN", "F3_PLAIN"},

  {"proton_PLAIN", "F1_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> u1Q2Q1=
{
  {"proton_PLAIN", "Tpi2_PLAIN", "F1_PLAIN"},
  {"proton_PLAIN", "Tpi2_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> u1=
{
  {"rnd_proton_PLAIN"},
  {"rnd_sub_PLAIN"},
};

void mcpico_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"b_norm", 20, 0, 2});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  ROOT::RDataFrame d( "tree", inputFiles.c_str() );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wUnity3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return 1; };

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
  for (auto &corr:u2Q1Q1)
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, qn, qn);
  }

  for (auto &corr:u2Q2)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, qn, qn);
  }

  for (auto &corr:u1Q2Q1)
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x2x1", P3::xxx(1, 2, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y2y1", P3::xyy(1, 2, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x2y1", P3::yxy(1, 2, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y2x1", P3::yyx(1, 2, 1), wSumWu3part, wy, qn, qn);
  }
  for (auto &corr:Q2Q1Q1)
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wUnity3part, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wUnity3part, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wUnity3part, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wUnity3part, wn, qn, qn);
  }

  {
    auto corr = u1.at(0);
    std::array<std::string, 1> qn{corr.at(0)};
    string corrName = corr.at(0) + ".x1";
    string corrName2 = corr.at(0) + ".x2";
    auto corr_function = [](const Qn::QVector &a) { return a.x(1); };
    auto corr_function2 = [](const Qn::QVector &a) { return a.x(2); };
    auto sumW = [](const Qn::QVector &a) { return a.sumweights(); };
    corrBuilder.AddCorrelationWithInternalReader(corrName, corr_function, sumW, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName2, corr_function2, sumW, wy, qn, qn);
  }

  {
    auto corr = u1.at(1);
    std::array<std::string, 1> qn{ corr.at(0) };
    string corrName=corr.at(0)+".x1";
    auto corr_function = []( const Qn::QVector &a ){ return a.x(1); };
    auto unity = [](const Qn::QVector &a) { return 1; };
    corrBuilder.AddCorrelationWithInternalReader(corrName, corr_function, unity, wn, qn, qn);
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
