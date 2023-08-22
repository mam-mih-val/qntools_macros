#include "QnDataFrame.hpp"

vector <vector<string>> Q1Q1{
  {"F1_RESCALED", "F2_RESCALED"},
  {"F1_RESCALED", "F3_RESCALED"},
  {"F2_RESCALED", "F3_RESCALED"},

  {"Tp_RESCALED", "F1_RESCALED"},
  {"Tp_RESCALED", "F2_RESCALED"},
  {"Tp_RESCALED", "F3_RESCALED"},

  {"Tpi_RESCALED", "F1_RESCALED"},
  {"Tpi_RESCALED", "F2_RESCALED"},
  {"Tpi_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u1Q1{
  {"proton_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED"},

  {"pi_pos_RESCALED", "F1_RESCALED"},
  {"pi_pos_RESCALED", "F2_RESCALED"},
  {"pi_pos_RESCALED", "F3_RESCALED"},

  {"pi_neg_RESCALED", "F1_RESCALED"},
  {"pi_neg_RESCALED", "F2_RESCALED"},
  {"pi_neg_RESCALED", "F3_RESCALED"},

  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"tru_proton_PLAIN", "psi_rp_PLAIN"},

  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},

  {"pi_neg_RESCALED", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
};

vector <vector<string>> u2Q1Q1{
  {"proton_RESCALED", "F1_RESCALED", "F3_RESCALED"},
  {"pi_pos_RESCALED", "F1_RESCALED", "F3_RESCALED"},
  {"pi_neg_RESCALED", "F1_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u2Q2{
  {"proton_RESCALED", "psi_rp_PLAIN"},
  {"pi_pos_RESCALED", "psi_rp_PLAIN"},
  {"pi_neg_RESCALED", "psi_rp_PLAIN"},

  {"tru_proton_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_pos_PLAIN", "psi_rp_PLAIN"},
  {"tru_pi_neg_PLAIN", "psi_rp_PLAIN"},
};


void mpd_fixed_target_correlate(string inputFiles="qn.root", string outputFile="corr.root")
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
  for (auto &corr:Q1Q1){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, qn, qn);
  }

  for (auto &corr:u1Q1){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr:u2Q2){
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, qn, qn);
  }

  for (auto &corr : u2Q1Q1 ){
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, qn, qn);
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
