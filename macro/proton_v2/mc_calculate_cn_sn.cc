#include "QnDataFrame.hpp"

std::string u1_vector{ "proton_PLAIN" };
std::string tru_p_vector{ "tru_proton_PLAIN" };

std::string f1_vector{ "F1_PLAIN" };
std::string f2_vector{ "F2_PLAIN" };
std::string f3_vector{ "F3_PLAIN" };
std::string f4_vector{ "F4_PLAIN" };

std::string tp_vector{ "Tpos_PLAIN" };
std::string tn_vector{ "Tneg_PLAIN" };

std::string psi_vector{ "psi_rp_PLAIN" };

std::vector < std::array<std::string, 1> > arr_u1 { 
  std::array<std::string, 1>{u1_vector},
  std::array<std::string, 1>{tru_p_vector},
};

std::vector < std::array<std::string, 1> > arr_Q1 {
  std::array<std::string, 1>{f1_vector},
  std::array<std::string, 1>{f2_vector},  
  std::array<std::string, 1>{f3_vector},
  std::array<std::string, 1>{f4_vector},
  std::array<std::string, 1>{tp_vector},
  std::array<std::string, 1>{tn_vector},
  std::array<std::string, 1>{psi_vector},
};

namespace P1 {
  inline auto mag( unsigned int h_a ) {
    return [ h_a ](const Qn::QVector &a ) {
      return sqrt( a.x(h_a)*a.x(h_a) + a.y(h_a)*a.y(h_a) );
    };
  }
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

  inline auto x2( unsigned int h_a ) {
    return [ h_a ](const Qn::QVector &a ) {
      return a.x(h_a)*a.x(h_a) - a.y(h_a)*a.y(h_a);
    };
  }

  inline auto y2( unsigned int h_a ) {
    return [ h_a ](const Qn::QVector &a ) {
      return 2 * a.x(h_a) * a.y(h_a);
    };
  }
}

void mc_calculate_cn_sn(string inputFiles="qn.root", string outputFile="CnSn.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 6, 0, 60});
  Qn::AxisD runIdAxis({ "runId", 240, 7100, 8300 });
  Qn::AxisD vtxXAxis({ "vtxX", 5, -1.0, 1.0 });
  Qn::AxisD vtxYAxis({ "vtxY", 5, -1.0, 1.0 });
  
  auto axes_correlation = Qn::MakeAxes(centAxis);
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

    for( auto i=size_t{1}; i<=8; ++i ){
      auto full_name = std::string{ corrName }.append("x").append( std::to_string(i) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::x(i), wSumWu1part, wy, corr, corr);
      full_name = std::string{ corrName }.append("y").append( std::to_string(i) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::y(i), wSumWu1part, wy, corr, corr);
      full_name = std::string{ corrName }.append("x").append( std::to_string(i) ).append( std::to_string(i*2) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::x2(i), wSumWu1part, wy, corr, corr);
      full_name = std::string{ corrName }.append("y").append( std::to_string(i) ).append( std::to_string(i*2) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::y2(i), wSumWu1part, wy, corr, corr);
    }
  }

  for ( auto &corr: arr_Q1 ){
    string corrName=corr.at(0);
    
    for( auto i=size_t{1}; i<=8; ++i ){
      auto full_name = std::string{ corrName }.append("x").append( std::to_string(i) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::x(i), wUnity1part, wn, corr, corr);
      full_name = std::string{ corrName }.append("y").append( std::to_string(i) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::y(i), wUnity1part, wn, corr, corr);
      full_name = std::string{ corrName }.append("x").append( std::to_string(i) ).append( std::to_string(i*2) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::x2(i), wUnity1part, wn, corr, corr);
      full_name = std::string{ corrName }.append("y").append( std::to_string(i) ).append( std::to_string(i*2) );
      corrBuilder.AddCorrelationWithInternalReader(full_name, P1::y2(i), wUnity1part, wn, corr, corr);
    }
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
