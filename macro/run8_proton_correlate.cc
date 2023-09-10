#include "QnDataFrame.hpp"

vector <vector<string>> Q1Q1_plain =
{
  {"F1_PLAIN", "F2_PLAIN"},
  {"F1_PLAIN", "F3_PLAIN"},
  {"F2_PLAIN", "F3_PLAIN"},

  {"Tneg_PLAIN", "F1_PLAIN"},
  {"Tneg_PLAIN", "F2_PLAIN"},
  {"Tneg_PLAIN", "F3_PLAIN"},

  {"Tpos_PLAIN", "F1_PLAIN"},
  {"Tpos_PLAIN", "F2_PLAIN"},
  {"Tpos_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> Q1Q1_recentered =
{
  {"F1_RECENTERED", "F2_RECENTERED"},
  {"F1_RECENTERED", "F3_RECENTERED"},
  {"F2_RECENTERED", "F3_RECENTERED"},

  {"Tneg_RESCALED", "F1_RECENTERED"},
  {"Tneg_RESCALED", "F2_RECENTERED"},
  {"Tneg_RESCALED", "F3_RECENTERED"},

  {"Tpos_RESCALED", "F1_RECENTERED"},
  {"Tpos_RESCALED", "F2_RECENTERED"},
  {"Tpos_RESCALED", "F3_RECENTERED"},
};

vector <vector<string>> Q1Q1_rescaled =
{
  {"F1_RESCALED", "F2_RESCALED"},
  {"F1_RESCALED", "F3_RESCALED"},
  {"F2_RESCALED", "F3_RESCALED"},

  {"Tneg_RESCALED", "F1_RESCALED"},
  {"Tneg_RESCALED", "F2_RESCALED"},
  {"Tneg_RESCALED", "F3_RESCALED"},

  {"Tpos_RESCALED", "F1_RESCALED"},
  {"Tpos_RESCALED", "F2_RESCALED"},
  {"Tpos_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u1Q1_plain=
{
  {"proton_RESCALED", "F1_PLAIN"},
  {"proton_RESCALED", "F2_PLAIN"},
  {"proton_RESCALED", "F3_PLAIN"},
};

vector <vector<string>> u1Q1_recentered=
{
  {"proton_RESCALED", "F1_RECENTERED"},
  {"proton_RESCALED", "F2_RECENTERED"},
  {"proton_RESCALED", "F3_RECENTERED"},

  {"proton700_RESCALED", "F1_RECENTERED"},
  {"proton700_RESCALED", "F2_RECENTERED"},
  {"proton700_RESCALED", "F3_RECENTERED"},

  {"proton400_RESCALED", "F1_RECENTERED"},
  {"proton400_RESCALED", "F2_RECENTERED"},
  {"proton400_RESCALED", "F3_RECENTERED"},
};

vector <vector<string>> u1Q1_rescaled=
{
  {"proton_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED"},

  {"proton700_RESCALED", "F1_RESCALED"},
  {"proton700_RESCALED", "F2_RESCALED"},
  {"proton700_RESCALED", "F3_RESCALED"},

  {"proton400_RESCALED", "F1_RESCALED"},
  {"proton400_RESCALED", "F2_RESCALED"},
  {"proton400_RESCALED", "F3_RESCALED"},
};


void run8_proton_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 8, 0, 40});
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

  for (auto &corr: Q1Q1_rescaled)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
  }

  for ( auto &corr: u1Q1_rescaled )
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
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
