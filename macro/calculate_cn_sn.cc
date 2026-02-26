#include "QnDataFrame.hpp"

std::string u1_vector{ "proton_PLAIN" };

std::string f1_vector{ "F1_PLAIN" };
std::string f2_vector{ "F2_PLAIN" };
std::string f3_vector{ "F3_PLAIN" };

std::string tp_vector{ "Tpos_PLAIN" };
std::string tn_vector{ "Tneg_PLAIN" };

std::vector < std::array<std::string, 1> > arr_u1 {
};

std::vector < std::array<std::string, 1> > arr_Q1 {
  std::array<std::string, 1>{u1_vector},
  std::array<std::string, 1>{f1_vector},
  std::array<std::string, 1>{f2_vector},  
  std::array<std::string, 1>{f3_vector},
  std::array<std::string, 1>{tp_vector},
  std::array<std::string, 1>{tn_vector},
};

namespace P1 {
  inline auto x( unsigned int h_a ) {
    return [ h_a ](const Qn::QVector &a ) {
      return a.x(h_a);
    };
  }

  inline auto y( unsigned int h_a ) {
    return [ h_a ](const Qn::QVector &a ) {
      return a.y(h_a);
    };
  }

  inline auto xx( unsigned int h_a, unsigned int h_b ) {
    return [ h_a, h_b ](const Qn::QVector &a ) {
      return a.x(h_a) * a.x(h_b);
    };
  }

  inline auto yy( unsigned int h_a, unsigned int h_b ) {
    return [ h_a, h_b ](const Qn::QVector &a ) {
      return a.y(h_a) * a.y(h_b);
    };
  }

  inline auto xy( unsigned int h_a, unsigned int h_b ) {
    return [ h_a, h_b ](const Qn::QVector &a ) {
      return a.x(h_a) * a.y(h_b);
    };
  }

  inline auto yx( unsigned int h_a, unsigned int h_b ) {
    return [ h_a, h_b ](const Qn::QVector &a ) {
      return a.y(h_a) * a.x(h_b);
    };
  }

}

void calculate_cn_sn(string inputFiles="qn.root", string outputFile="CnSn.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 6, 0, 60});
  Qn::AxisD runIdAxis({ "runId", 17, 6600, 8300 });
  auto axes_correlation = Qn::MakeAxes(centAxis, runIdAxis);
  std::string treename = "tree";
  auto* chain = new TChain( treename.c_str() );
  chain->AddFile( inputFiles.c_str() );
  if( chain->GetEntries() <= 0 )
    return;
  ROOT::RDataFrame d( *chain );
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity1part = [](const Qn::QVector &a) { return 1; };
  auto wUnity2part = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu1part = [](const Qn::QVector &a) { return a.sumweights(); };
  auto wSumWu2part = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};

  // *******************************************
  // -------------- Q1 RESCALED --------------
  // *******************************************

  for ( auto &corr: arr_u1 ){
    string corrName=corr.at(0);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1", P1::x(1), wSumWu1part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1", P1::y(1), wSumWu1part, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2", P1::x(2), wSumWu1part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2", P1::y(2), wSumWu1part, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3", P1::x(3), wSumWu1part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3", P1::y(3), wSumWu1part, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4", P1::x(4), wSumWu1part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4", P1::y(4), wSumWu1part, wy, corr, corr);

  }

  for ( auto &corr: arr_Q1 ){
    string corrName=corr.at(0);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1", P1::x(1), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1", P1::y(1), wUnity1part, wn, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2", P1::x(2), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2", P1::y(2), wUnity1part, wn, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3", P1::x(3), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3", P1::y(3), wUnity1part, wn, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x4", P1::x(4), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y4", P1::y(4), wUnity1part, wn, corr, corr);
  }

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
