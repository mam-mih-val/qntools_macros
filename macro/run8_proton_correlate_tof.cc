#include "QnDataFrame.hpp"

vector <vector<string>> u1_rescaled_Q1_rescaled=
{
  {"proton_RESCALED", "F1_RECENTERED"},
  {"proton_RESCALED", "F2_RECENTERED"},
  {"proton_RESCALED", "F3_RECENTERED"},

  {"proton_RECENTERED", "F1_RECENTERED"},
  {"proton_RECENTERED", "F2_RECENTERED"},
  {"proton_RECENTERED", "F3_RECENTERED"},

  {"proton_PLAIN", "F1_RECENTERED"},
  {"proton_PLAIN", "F2_RECENTERED"},
  {"proton_PLAIN", "F3_RECENTERED"},

  // {"proton_400_RESCALED", "F1_RECENTERED"},
  // {"proton_400_RESCALED", "F2_RECENTERED"},
  // {"proton_400_RESCALED", "F3_RECENTERED"},
  
  // {"proton_700_RESCALED", "F1_RECENTERED"},
  // {"proton_700_RESCALED", "F2_RECENTERED"},
  // {"proton_700_RESCALED", "F3_RECENTERED"},
};

vector <vector<string>> Q1_rescaled_Q1_rescaled =
{
  {"F1_RECENTERED", "F2_RECENTERED"},
  {"F1_RECENTERED", "F3_RECENTERED"},
  {"F2_RECENTERED", "F3_RECENTERED"},

  // {"Tneg_RESCALED", "F1_RECENTERED"},
  // {"Tneg_RESCALED", "F2_RECENTERED"},
  // {"Tneg_RESCALED", "F3_RECENTERED"},

  // {"Tpos_RESCALED", "F1_RECENTERED"},
  // {"Tpos_RESCALED", "F2_RECENTERED"},
  // {"Tpos_RESCALED", "F3_RECENTERED"},
};

// vector <vector<string>> u2Q1Q1_rescaled=
// {
//   {"proton_RESCALED", "F1_RECENTERED", "F1_RECENTERED"},
//   {"proton_RESCALED", "F2_RECENTERED", "F2_RECENTERED"},
//   {"proton_RESCALED", "F3_RECENTERED", "F3_RECENTERED"},

//   {"proton_400_RESCALED", "F1_RECENTERED", "F1_RECENTERED"},
//   {"proton_400_RESCALED", "F2_RECENTERED", "F2_RECENTERED"},
//   {"proton_400_RESCALED", "F3_RECENTERED", "F3_RECENTERED"},

//   {"proton_700_RESCALED", "F1_RECENTERED", "F1_RECENTERED"},
//   {"proton_700_RESCALED", "F2_RECENTERED", "F2_RECENTERED"},
//   {"proton_700_RESCALED", "F3_RECENTERED", "F3_RECENTERED"},
// };

vector <vector<string>> u2Q1Q1_rescaled=
{
  {"proton_RESCALED", "F1_RECENTERED", "F2_RECENTERED"},
  {"proton_RESCALED", "F2_RECENTERED", "F3_RECENTERED"},
  {"proton_RESCALED", "F1_RECENTERED", "F3_RECENTERED"},

  {"proton_RECENTERED", "F1_RECENTERED", "F2_RECENTERED"},
  {"proton_RECENTERED", "F2_RECENTERED", "F3_RECENTERED"},
  {"proton_RECENTERED", "F1_RECENTERED", "F3_RECENTERED"},

  {"proton_PLAIN", "F1_RECENTERED", "F2_RECENTERED"},
  {"proton_PLAIN", "F2_RECENTERED", "F3_RECENTERED"},
  {"proton_PLAIN", "F1_RECENTERED", "F3_RECENTERED"},

  {"proton_RESCALED", "F1_RECENTERED", "F1_RECENTERED"},
  {"proton_RESCALED", "F2_RECENTERED", "F2_RECENTERED"},
  {"proton_RESCALED", "F3_RECENTERED", "F3_RECENTERED"},

  {"proton_RECENTERED", "F1_RECENTERED", "F1_RECENTERED"},
  {"proton_RECENTERED", "F2_RECENTERED", "F2_RECENTERED"},
  {"proton_RECENTERED", "F3_RECENTERED", "F3_RECENTERED"},

  {"proton_PLAIN", "F1_RECENTERED", "F1_RECENTERED"},
  {"proton_PLAIN", "F2_RECENTERED", "F2_RECENTERED"},
  {"proton_PLAIN", "F3_RECENTERED", "F3_RECENTERED"},

  // {"proton_400_RESCALED", "F1_RECENTERED", "F2_RECENTERED"},
  // {"proton_400_RESCALED", "F2_RECENTERED", "F3_RECENTERED"},
  // {"proton_400_RESCALED", "F1_RECENTERED", "F3_RECENTERED"},

  // {"proton_700_RESCALED", "F1_RECENTERED", "F2_RECENTERED"},
  // {"proton_700_RESCALED", "F2_RECENTERED", "F3_RECENTERED"},
  // {"proton_700_RESCALED", "F1_RECENTERED", "F3_RECENTERED"},
};


vector <vector<string>> u3Q1Q1Q1_rescaled=
{
  {"proton_RESCALED", "F1_RECENTERED", "F2_RECENTERED", "F3_RECENTERED"},
  {"proton_RECENTERED", "F1_RECENTERED", "F2_RECENTERED", "F3_RECENTERED"},
  {"proton_PLAIN", "F1_RECENTERED", "F2_RECENTERED", "F3_RECENTERED"},

  {"proton_RESCALED", "F1_RECENTERED", "F1_RECENTERED", "F1_RECENTERED"},
  {"proton_RECENTERED", "F1_RECENTERED", "F1_RECENTERED", "F1_RECENTERED"},
  {"proton_PLAIN", "F1_RECENTERED", "F1_RECENTERED", "F1_RECENTERED"},

  {"proton_RESCALED", "F2_RECENTERED", "F2_RECENTERED", "F2_RECENTERED"},
  {"proton_RECENTERED", "F2_RECENTERED", "F2_RECENTERED", "F2_RECENTERED"},
  {"proton_PLAIN", "F2_RECENTERED", "F2_RECENTERED", "F2_RECENTERED"},

  {"proton_RESCALED", "F3_RECENTERED", "F3_RECENTERED", "F3_RECENTERED"},
  {"proton_RECENTERED", "F3_RECENTERED", "F3_RECENTERED", "F3_RECENTERED"},
  {"proton_PLAIN", "F3_RECENTERED", "F3_RECENTERED", "F3_RECENTERED"},

  // {"proton_RESCALED", "F1_RECENTERED", "F3_RECENTERED", "Tneg_RESCALED"},
  // {"proton_RESCALED", "F1_RECENTERED", "F2_RECENTERED", "Tneg_RESCALED"},
  // {"proton_RESCALED", "F2_RECENTERED", "F3_RECENTERED", "Tneg_RESCALED"},

  // {"proton_RESCALED", "F1_RECENTERED", "F1_RECENTERED", "F1_RECENTERED"},
  // {"proton_RESCALED", "F2_RECENTERED", "F2_RECENTERED", "F2_RECENTERED"},
  // {"proton_RESCALED", "F3_RECENTERED", "F3_RECENTERED", "F3_RECENTERED"},
};

namespace P4{

  inline auto x1y1y1y1( unsigned int h_a, unsigned int h_b, unsigned int h_c ) {
    return [ h_a, h_b, h_c ](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c ) {
      return a.x(h_a)*a.y(h_a) * b.y(h_b) * c.y(h_c) ;
    };
  }

  inline auto x1x1y1y1( unsigned int h_a, unsigned int h_b, unsigned int h_c ) {
    return [ h_a, h_b, h_c ](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c ) {
      return a.x(h_a)*a.x(h_a) * b.y(h_b) * c.y(h_c) ;
    };
  }

  inline auto y1y1y1y1( unsigned int h_a, unsigned int h_b, unsigned int h_c ) {
    return [ h_a, h_b, h_c ](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c ) {
      return a.y(h_a)*a.y(h_a) * b.y(h_b) * c.y(h_c) ;
    };
  }
  
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
void run8_proton_correlate_tof(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 4, 0, 40});
  auto axes_correlation = Qn::MakeAxes(centAxis);
  std::string treename = "tree";
  auto* chain = new TChain( treename.c_str() );
  chain->AddFile( inputFiles.c_str() );
  if( chain->GetEntries() <= 0 )
    return;
  ROOT::RDataFrame d( *chain );
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
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1", P2::xx(2, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1", P2::yy(2, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1", P2::xy(2, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1", P2::yx(2, 1), wSumWu, wy, qn, qn);
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

  for ( auto &corr: u2Q1Q1_rescaled )
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    // diagonal/non-zero
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, qn, qn);
    // // non-diagonal/zero 
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1x1", P3::yxx(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1y1", P3::yyy(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1y1", P3::xxy(2, 1, 1), wSumWu3part, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1x1", P3::xyx(2, 1, 1), wSumWu3part, wy, qn, qn);

    //double harmonic
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1y1y1", P4::x1y1y1y1(1, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1y1y1", P4::x1x1y1y1(1, 1, 1), wSumWu3part, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1y1y1", P4::x1x1y1y1(1, 1, 1), wSumWu3part, wy, qn, qn);
  }

  // for ( auto &corr: u3Q1Q1Q1_rescaled )
  // {
  //   std::array<std::string, 4> qn{corr.at(0), corr.at(1), corr.at(2), corr.at(3)};
  //   string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2)+"."+corr.at(3);
  //   corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1x1x1", P4::xxxx(3, 1, 1, 1), wSumWu4part, wy, qn, qn);
  //   corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1y1y1", P4::yyyy(3, 1, 1, 1), wSumWu4part, wy, qn, qn);

  //   corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1x1x1", P4::xxxx(1, 1, 1, 1), wSumWu4part, wy, qn, qn);
  //   corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1y1y1", P4::yyyy(1, 1, 1, 1), wSumWu4part, wy, qn, qn);
  // }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open(outputFile.c_str(), "RECREATE"), []( TFile* f ){ f->Close(); } };
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  // corrFile->Close();
}
