#include "QnDataFrame.hpp"

vector <vector<string>> u1_rescaled_Q1_rescaled=
{
  {"deuteron_RESCALED", "F1_RESCALED"},
  {"deuteron_RESCALED", "F2_RESCALED"},
  {"deuteron_RESCALED", "F3_RESCALED"},
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
  {"proton_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "F3_RESCALED"},

  {"proton_RESCALED", "F1_RESCALED", "Tneg_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "Tneg_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "Tneg_RESCALED"},
  
  {"proton_RESCALED", "F1_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F3_RESCALED"},
  {"proton_RESCALED", "F1_RESCALED", "F3_RESCALED"},
};

vector <vector<string>> u3Q1Q1Q1_rescaled=
{
  {"proton_RESCALED", "F1_RESCALED", "F2_RESCALED", "F3_RESCALED"},
  {"proton_RESCALED", "F1_RESCALED", "F3_RESCALED", "Tneg_RESCALED"},
  {"proton_RESCALED", "F1_RESCALED", "F2_RESCALED", "Tneg_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F3_RESCALED", "Tneg_RESCALED"},

  {"proton_RESCALED", "F1_RESCALED", "F1_RESCALED", "F1_RESCALED"},
  {"proton_RESCALED", "F2_RESCALED", "F2_RESCALED", "F2_RESCALED"},
  {"proton_RESCALED", "F3_RESCALED", "F3_RESCALED", "F3_RESCALED"},

};

namespace P4{

  inline auto xxxx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.x(h_b) * c.x(h_c) * d.x(h_d);
    };
  }

  inline auto xxyy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.x(h_b) * c.y(h_c) * d.y(h_d);
    };
  }

  inline auto xyxy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.y(h_b) * c.x(h_c) * d.y(h_d);
    };
  }


  inline auto xyyx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.x(h_a) * b.y(h_b) * c.y(h_c) * d.x(h_d);
    };
  }

  inline auto yxxy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.x(h_b) * c.x(h_c) * d.y(h_d);
    };
  }

  inline auto yxyx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.x(h_b) * c.y(h_c) * d.x(h_d);
    };
  }

  inline auto yyxx(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.y(h_b) * c.x(h_c) * d.x(h_d);
    };
  }

  inline auto yyyy(unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d) {
    return [h_a, h_b, h_c, h_d](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) {
      return a.y(h_a) * b.y(h_b) * c.y(h_c) * d.y(h_d);
    };
  }

}
void run8_deuteron_correlate(string inputFiles="qn.root", string outputFile="corr.root")
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
  auto wSumWu4part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) { return a.sumweights(); };

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
