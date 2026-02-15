#include "matrix.h"
#include "DataContainer.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"
#include "QVector.hpp"
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

std::tuple< std::vector<Qn::DataContainerStatCalculate>, std::vector<Qn::DataContainerStatCalculate> > ReadCnSn( std::string str_vec_name, TFile* calib_file ){
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_s = std::vector<Qn::DataContainerStatCalculate>{};
  vec_c.reserve(4);
  vec_s.reserve(4);
  for( auto i=size_t{0}; i<4; ++i ){
    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centrality"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.emplace_back( *tmp );

    corr_name = str_vec_name+".y"+std::to_string(i+1)+"centrality"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.emplace_back( *tmp );
  }

  return {vec_c, vec_s};
}

void run8_proton_decompose(std::string in_file_name, std::string in_calib_file){

  const std::string f1_name {"F1_PLAIN"};
  const std::string f2_name {"F2_PLAIN"};
  const std::string f3_name {"F3_PLAIN"};
  const std::string tp_name {"Tpos_PLAIN"};
  const std::string tn_name {"Tneg_PLAIN"};
  const std::string proton_name {"proton_PLAIN"};

  auto calib_file = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( in_calib_file.c_str(), "READ" ), [](auto f){f ->Close(); } };
  assert(calib_file);

  auto [vec_c_f1, vec_s_f1] = ReadCnSn(f1_name, calib_file.get());
  auto [vec_c_f2, vec_s_f2] = ReadCnSn(f2_name, calib_file.get());
  auto [vec_c_f3, vec_s_f3] = ReadCnSn(f3_name, calib_file.get());
  auto [vec_c_tp, vec_s_tp] = ReadCnSn(tp_name, calib_file.get());
  auto [vec_c_tn, vec_s_tn] = ReadCnSn(tn_name, calib_file.get());
  auto [vec_c_p, vec_s_p] = ReadCnSn(proton_name, calib_file.get());

  const auto correction_generator = []( const std::vector<Qn::DataContainerStatCalculate>& vec_c, const std::vector<Qn::DataContainerStatCalculate>& vec_s ){
    return [&vec_c, &vec_s]( Qn::DataContainerQVector qvec, Double_t centrality ) -> Qn::DataContainerQVector {
      if( centrality < 1 || centrality > 60 )
        return qvec;
      auto new_qvec = qvec;
      auto c_axis = vec_c[0].GetAxis( "centrality" ); 
      auto c_bin = c_axis.FindBin( centrality );
      auto bin_lo = c_axis.GetLowerBinEdge( c_bin );
      auto bin_hi = c_axis.GetUpperBinEdge( c_bin );
      auto new_c_axis = Qn::AxisD{ "centrality", 1, bin_lo, bin_hi };

      auto vec_c_c = vec_c;
      auto vec_s_c = vec_s;

      for( auto& el : vec_c_c ) {
        if( el.GetAxes().size() > 1 )
          el = el.Select( new_c_axis );
        else
          el = el.Rebin( new_c_axis );
      }
      for( auto& el : vec_s_c ) {
        if( el.GetAxes().size() > 1 )
          el = el.Select( new_c_axis );
        else
          el = el.Rebin( new_c_axis );
      }

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        auto c1 = vec_c_c[0].At(i).Mean();
        auto c2 = vec_c_c[1].At(i).Mean();
        auto c3 = vec_c_c[2].At(i).Mean();
        auto c4 = vec_c_c[3].At(i).Mean();

        auto s1 = vec_s_c[0].At(i).Mean();
        auto s2 = vec_s_c[1].At(i).Mean();
        auto s3 = vec_s_c[2].At(i).Mean();
        auto s4 = vec_s_c[3].At(i).Mean();

        auto x1_old = qvec.At(i).x(1) - c1;
        auto y1_old = qvec.At(i).y(1) - s1;
        auto x2_old = qvec.At(i).x(2) - c2;
        auto y2_old = qvec.At(i).y(2) - s2;

        auto M = Matrix<double, 2>{};
        M[0] = std::array<double, 2>{ 1+c2,     s2 };
        M[1] = std::array<double, 2>{ s2,     1-c2 };
        // M[0] = std::array<double, 4>{ 1+c2,     s2,   c1+c3,  s1+s3 };
        // M[1] = std::array<double, 4>{ s2,     1-c2,   s3-s1,  c1-c3 };
        // M[2] = std::array<double, 4>{ c1+c3, s3-s1,   1+c4,      s4 };
        // M[3] = std::array<double, 4>{ s1+s3, c1-c3,     s4,    1-c4 };

        auto Minv = Inverse(M);
        
        auto x1_new = Minv[0][0]*x1_old + Minv[1][0]*y1_old;
        auto y1_new = Minv[0][1]*x1_old + Minv[1][1]*y1_old;

        // auto x1_new = Minv[0][0]*x1_old + Minv[1][0]*y1_old + Minv[2][0]*x2_old + Minv[3][0]*y2_old;
        // auto y1_new = Minv[0][1]*x1_old + Minv[1][1]*y1_old + Minv[2][1]*x2_old + Minv[3][1]*y2_old;
        // auto x2_new = Minv[0][2]*x1_old + Minv[1][2]*y1_old + Minv[2][2]*x2_old + Minv[3][2]*y2_old;
        // auto y2_new = Minv[0][3]*x1_old + Minv[1][3]*y1_old + Minv[2][3]*x2_old + Minv[3][3]*y2_old;
        
        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_old, y2_old );
      }

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Define("F1_DECOMPOSED", correction_generator(vec_c_f1, vec_s_f1), { f1_name, "centrality" } )
    .Define("F2_DECOMPOSED", correction_generator(vec_c_f2, vec_s_f2), { f2_name, "centrality" } )
    .Define("F3_DECOMPOSED", correction_generator(vec_c_f3, vec_s_f3), { f3_name, "centrality" } )

    .Define("Tpos_DECOMPOSED", correction_generator(vec_c_tp, vec_s_tp), { tp_name, "centrality" } )
    .Define("Tneg_DECOMPOSED", correction_generator(vec_c_tn, vec_s_tn), { tn_name, "centrality" } )
    
    .Define("proton_DECOMPOSED", correction_generator(vec_c_p, vec_s_p), { proton_name, "centrality" } )
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, tp{}, tn{}, p{};
  double cent{};
  tree->Branch( "centrality", &cent );

  tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );

  dd.Foreach( [tree, &cent, &f1, &f2, &f3, &tp, &tn, &p]( 
    double centrality,
    Qn::DataContainerQVector f1_ev, 
    Qn::DataContainerQVector f2_ev, 
    Qn::DataContainerQVector f3_ev, 
    Qn::DataContainerQVector tp_ev, 
    Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev ) mutable {
    cent = centrality;
    f1 = f1_ev;
    f2 = f2_ev;
    f3 = f3_ev;
    tp = tp_ev;
    tn = tn_ev;
    p = p_ev;
    tree->Fill();
  }, 
  std::vector<std::string>{ 
    "centrality",
    "F1_DECOMPOSED",
    "F2_DECOMPOSED",
    "F3_DECOMPOSED",
    "Tpos_DECOMPOSED",
    "Tneg_DECOMPOSED",
    "proton_DECOMPOSED",
  } );
  file_out->cd();
  tree->Write();

}