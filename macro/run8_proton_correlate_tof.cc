#include "QnDataFrame.hpp"

std::string u1_vector{ "proton_DECOMPOSED" };

std::string f1_vector{ "F1_DECOMPOSED" };
std::string f2_vector{ "F2_DECOMPOSED" };
std::string f3_vector{ "F3_DECOMPOSED" };

std::string tp_vector{ "Tpos_DECOMPOSED" };
std::string tn_vector{ "Tneg_DECOMPOSED" };

std::vector < std::array<std::string, 1> > arr_u1 {
  std::array<std::string, 1>{u1_vector}
};

std::vector < std::array<string, 1> > arr_Q1 {
  std::array<std::string, 1>{f1_vector},
  std::array<std::string, 1>{f2_vector},  
  std::array<std::string, 1>{f3_vector},
};

std::vector < std::array<string, 2> > arr_u1Q1 {
  std::array<std::string, 2>{u1_vector, f1_vector},
  std::array<std::string, 2>{u1_vector, f2_vector},
  std::array<std::string, 2>{u1_vector, f3_vector},
};

std::vector < std::array<string, 2> > arr_Q1Q1 {
  std::array<std::string, 2>{f1_vector, f2_vector},
  std::array<std::string, 2>{f1_vector, f3_vector},
  std::array<std::string, 2>{f2_vector, f3_vector},
  
  std::array<std::string, 2>{f1_vector, tp_vector},
  std::array<std::string, 2>{f2_vector, tp_vector},
  std::array<std::string, 2>{f3_vector, tp_vector},

  std::array<std::string, 2>{f1_vector, tn_vector},
  std::array<std::string, 2>{f2_vector, tn_vector},
  std::array<std::string, 2>{f3_vector, tn_vector},
};

std::vector < std::array<string, 3> > arr_u2Q1Q1 {
  std::array<std::string, 3>{u1_vector, f1_vector, f2_vector},
  std::array<std::string, 3>{u1_vector, f1_vector, f3_vector},
  std::array<std::string, 3>{u1_vector, f2_vector, f3_vector},
};

std::vector < std::array<string, 4> > arr_u3Q1Q1Q1 {
  std::array<std::string, 4>{u1_vector, f1_vector, f2_vector, f3_vector},
};

std::vector < std::array<string, 3> > arr_Q1Q1Q1 {
  std::array<std::string, 3>{f1_vector, f2_vector, f3_vector},
};

namespace P4{
  inline auto yyyy( unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d ) {
    return [ h_a, h_b, h_c, h_d ]( const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d ) {
      return a.y(h_a) * b.y(h_b) * c.y(h_c) * d.y(h_d);
    };
  }

  inline auto xyyy( unsigned int h_a, unsigned int h_b, unsigned int h_c, unsigned int h_d ) {
    return [ h_a, h_b, h_c, h_d ]( const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d ) {
      return a.x(h_a) * b.y(h_b) * c.y(h_c) * d.y(h_d);
    };
  }
}

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

}

void run8_proton_correlate_tof(string inputFiles="qn.root", string outputFile="corr.root")
{
  int nSamples = 100;
  Qn::AxisD centAxis({"centrality", 6, 0, 60});
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
  auto wUnity1part = [](const Qn::QVector &a) { return 1; };
  auto wUnity3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return 1; };
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu1part = [](const Qn::QVector &a) { return a.sumweights(); };
  auto wSumWu3part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wSumWu4part = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c, const Qn::QVector &d) { return a.sumweights(); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};

  // *******************************************
  // -------------- Q1 RESCALED --------------
  // *******************************************

  for ( auto &corr: arr_u1 )
  {
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

  for ( auto &corr: arr_Q1 )
  {
    string corrName = corr.at(0);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1", P1::x(1), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1", P1::y(1), wUnity1part, wn, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2", P1::x(2), wUnity1part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2", P1::y(2), wUnity1part, wn, corr, corr);
  }

  for ( auto &corr: arr_u1Q1 )
  {
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wSumWu, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wSumWu, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1", P2::xx(2, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1", P2::yy(2, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1", P2::xy(2, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1", P2::yx(2, 1), wSumWu, wy, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1", P2::xx(3, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1", P2::yy(3, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1", P2::xy(3, 1), wSumWu, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1", P2::yx(3, 1), wSumWu, wy, corr, corr);
  }

  for (auto &corr: arr_Q1Q1)
  {
    string corrName=corr.at(0)+"."+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1", P2::xx(1, 1), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1", P2::yy(1, 1), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1", P2::xy(1, 1), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1", P2::yx(1, 1), wUnity, wn, corr, corr);

    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x2", P2::xx(2, 2), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y2", P2::yy(2, 2), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y2", P2::xy(2, 2), wUnity, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x2", P2::yx(2, 2), wUnity, wn, corr, corr);
  }

  for ( auto &corr: arr_u2Q1Q1 )
  {
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    // diagonal/non-zero
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1x1", P3::xxx(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1y1", P3::xyy(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1y1", P3::yxy(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1x1", P3::yyx(2, 1, 1), wSumWu3part, wy, corr, corr);
    // non-diagonal/zero 
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2x1x1", P3::yxx(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y2y1y1", P3::yyy(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2x1y1", P3::xxy(2, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x2y1x1", P3::xyx(2, 1, 1), wSumWu3part, wy, corr, corr);
    // third harmonic
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1x1", P3::xxx(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1y1", P3::xyy(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1y1", P3::yxy(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1x1", P3::yyx(3, 1, 1), wSumWu3part, wy, corr, corr);
    // non-diagonal/zero 
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3x1x1", P3::yxx(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1y1", P3::yyy(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3x1y1", P3::xxy(3, 1, 1), wSumWu3part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1x1", P3::xyx(3, 1, 1), wSumWu3part, wy, corr, corr);
  }

  for ( auto &corr: arr_Q1Q1Q1 )
  {
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2);
    // diagonal/non-zero
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1x1", P3::xxx(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1y1", P3::xyy(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1y1", P3::yxy(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1x1", P3::yyx(1, 1, 1), wUnity3part, wn, corr, corr);
    // non-diagonal/zero 
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1x1x1", P3::yxx(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1y1", P3::yyy(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1x1y1", P3::xxy(1, 1, 1), wUnity3part, wn, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1x1", P3::xyx(1, 1, 1), wUnity3part, wn, corr, corr);
  }

  for ( auto &corr: arr_u3Q1Q1Q1 )
  {
    string corrName=corr.at(0)+"."+corr.at(1)+"."+corr.at(2)+"."+corr.at(3);
    // diagonal/non-zero
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y1y1y1y1", P4::yyyy(1, 1, 1, 1), wSumWu4part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x1y1y1y1", P4::xyyy(1, 1, 1, 1), wSumWu4part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".x3y1y1y1", P4::yyyy(3, 1, 1, 1), wSumWu4part, wy, corr, corr);
    corrBuilder.AddCorrelationWithInternalReader(corrName+".y3y1y1y1", P4::xyyy(3, 1, 1, 1), wSumWu4part, wy, corr, corr);
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
