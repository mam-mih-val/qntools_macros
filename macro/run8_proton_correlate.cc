#include "QnDataFrame.hpp"

// ********************************
// ----------- Q1_plain ----------- 
// ********************************
vector <vector<string>> u1_plain_Q1_plain=
{
  {"proton_PLAIN", "F1_PLAIN"},
  {"proton_PLAIN", "F2_PLAIN"},
  {"proton_PLAIN", "F3_PLAIN"},
};

vector <vector<string>> u1_recentered_Q1_plain=
{
  {"proton_RECENTERED", "F1_PLAIN"},
  {"proton_RECENTERED", "F2_PLAIN"},
  {"proton_RECENTERED", "F3_PLAIN"},
};

vector <vector<string>> u1_rescaled_Q1_plain=
{
  {"proton_RESCALED", "F1_PLAIN"},
  {"proton_RESCALED", "F2_PLAIN"},
  {"proton_RESCALED", "F3_PLAIN"},
};

vector <vector<string>> Q1_plain_Q1_plain =
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

// *************************************
// ----------- Q1_recentered ----------- 
// *************************************

vector <vector<string>> u1_plain_Q1_recentered=
{
  {"proton_PLAIN", "F1_RECENTERED"},
  {"proton_PLAIN", "F2_RECENTERED"},
  {"proton_PLAIN", "F3_RECENTERED"},
};

vector <vector<string>> u1_recentered_Q1_recentered=
{
  {"proton_RECENTERED", "F1_RECENTERED"},
  {"proton_RECENTERED", "F2_RECENTERED"},
  {"proton_RECENTERED", "F3_RECENTERED"},
};

vector <vector<string>> u1_rescaled_Q1_recentered=
{
  {"proton_RESCALED", "F1_RECENTERED"},
  {"proton_RESCALED", "F2_RECENTERED"},
  {"proton_RESCALED", "F3_RECENTERED"},
};

vector <vector<string>> Q1_recentered_Q1_recentered =
{
  {"F1_RECENTERED", "F2_RECENTERED"},
  {"F1_RECENTERED", "F3_RECENTERED"},
  {"F2_RECENTERED", "F3_RECENTERED"},

  {"Tneg_RECENTERED", "F1_RECENTERED"},
  {"Tneg_RECENTERED", "F2_RECENTERED"},
  {"Tneg_RECENTERED", "F3_RECENTERED"},

  {"Tpos_RECENTERED", "F1_RECENTERED"},
  {"Tpos_RECENTERED", "F2_RECENTERED"},
  {"Tpos_RECENTERED", "F3_RECENTERED"},
};

// *************************************
// ----------- Q1_rescaled ----------- 
// *************************************

vector <vector<string>> u1_plain_Q1_rescaled=
{
  {"proton_PLAIN", "F1_RESCALED"},
  {"proton_PLAIN", "F2_RESCALED"},
  {"proton_PLAIN", "F3_RESCALED"},
};

vector <vector<string>> u1_recentered_Q1_rescaled=
{
  {"proton_RECENTERED", "F1_RESCALED"},
  {"proton_RECENTERED", "F2_RESCALED"},
  {"proton_RECENTERED", "F3_RESCALED"},
};

vector <vector<string>> u1_rescaled_Q1_rescaled=
{
  {"proton_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED"},

  {"proton_r_RESCALED", "F1_RESCALED"},
  {"proton_r_RESCALED", "F2_RESCALED"},
  {"proton_r_RESCALED", "F3_RESCALED"},

  {"proton_l_RESCALED", "F1_RESCALED"},
  {"proton_l_RESCALED", "F2_RESCALED"},
  {"proton_l_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> Q1_rescaled_Q1_rescaled =
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



vector <vector<string>> u2Q1Q1_rescaled=
{
  {"proton_RESCALED", "F1_RESCALED", "F3_RESCALED"},

  // {"pi_pos_RESCALED", "F1_RESCALED"},
  // {"pi_pos_RESCALED", "F2_RESCALED"},
  // {"pi_pos_RESCALED", "F3_RESCALED"},

  // {"pi_neg_RESCALED", "F1_RESCALED"},
  // {"pi_neg_RESCALED", "F2_RESCALED"},
  // {"pi_neg_RESCALED", "F3_RESCALED"},

  // {"deuteron_RESCALED", "F1_RESCALED"},
  // {"deuteron_RESCALED", "F2_RESCALED"},
  // {"deuteron_RESCALED", "F3_RESCALED"},
};


void run8_proton_correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 4, 0, 40});
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

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};

  // *******************************************
  // -------------- Q1 RESCALED --------------
  // *******************************************

  for ( auto &corr: u1_rescaled_Q1_rescaled )
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr: Q1_rescaled_Q1_rescaled)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
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
