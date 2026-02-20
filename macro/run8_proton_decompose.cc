#include "DataContainer.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"
#include "QVector.hpp"
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

Qn::DataContainer<Eigen::Matrix4d, Qn::AxisD> MakeCorrectionMatrix(const std::vector<Qn::DataContainerStatCalculate>& vec_c, const std::vector<Qn::DataContainerStatCalculate>& vec_s, const std::vector<Qn::DataContainerStatCalculate>& vec_cov ){
  auto axes = vec_c[0].GetAxes();
  Qn::DataContainer<Eigen::Matrix4d, Qn::AxisD> corr_matrix{axes};
  for( auto i = size_t{0}; i<vec_c[0].size(); ++i ){
    auto c1 = vec_c[0].At(i).Mean();
    auto c2 = vec_c[1].At(i).Mean();
    auto c3 = vec_c[2].At(i).Mean();
    auto c4 = vec_c[3].At(i).Mean();

    auto s1 = vec_s[0].At(i).Mean();
    auto s2 = vec_s[1].At(i).Mean();
    auto s3 = vec_s[2].At(i).Mean();
    auto s4 = vec_s[3].At(i).Mean();

    auto M11 = vec_cov[0].At(i).Mean() - c1*c1;
    auto M12 = vec_cov[1].At(i).Mean() - c1*s1;
    auto M13 = vec_cov[2].At(i).Mean() - c1*c2;
    auto M14 = vec_cov[3].At(i).Mean() - c1*s2;
    
    auto M22 = vec_cov[4].At(i).Mean() - s1*s1;
    auto M23 = vec_cov[5].At(i).Mean() - s1*c2;
    auto M24 = vec_cov[6].At(i).Mean() - s1*s2;

    auto M33 = vec_cov[7].At(i).Mean() - c2*c2;
    auto M34 = vec_cov[8].At(i).Mean() - c2*s2;
    
    auto M44 = vec_cov[9].At(i).Mean() - s2*s2;

    auto M = Eigen::Matrix4d{
      { M11, M12, M13, M14 },
      { M12, M22, M23, M24 },
      { M13, M23, M33, M34 },
      { M14, M24, M34, M44 },
    };

    auto solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d>{ M };
    auto Minv = solver.operatorInverseSqrt() * pow(2.0, 0.5);
    corr_matrix.At(i) = Minv;
  }
  return corr_matrix;
}


std::tuple< std::vector<Qn::DataContainerStatCalculate>, std::vector<Qn::DataContainerStatCalculate>, std::vector<Qn::DataContainerStatCalculate> > ReadCnSn( std::string str_vec_name, TFile* calib_file ){
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_s = std::vector<Qn::DataContainerStatCalculate>{};
  vec_c.reserve(4);
  vec_s.reserve(4);
  for( auto i=size_t{0}; i<4; ++i ){
    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centralityrunId"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.emplace_back( *tmp );
    corr_name = str_vec_name+".y"+std::to_string(i+1)+"centralityrunId"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.emplace_back( *tmp );
  }

  auto vec_cov = std::vector<Qn::DataContainerStatCalculate>{};
  auto components = std::vector<std::string>{ "x1", "y1", "x2", "y2" };
  for( auto i=size_t{0}; i<components.size(); ++i ){
    for( auto j=size_t{i}; j<components.size(); ++j ){
      auto corr_name = str_vec_name + "." + components.at(i) + components.at(j) + "centralityrunId"s;
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );
    }
  }
  return {vec_c, vec_s, vec_cov};
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

  auto [vec_c_f1, vec_s_f1, vec_cov_f1] = ReadCnSn(f1_name, calib_file.get());
  auto [vec_c_f2, vec_s_f2, vec_cov_f2] = ReadCnSn(f2_name, calib_file.get());
  auto [vec_c_f3, vec_s_f3, vec_cov_f3] = ReadCnSn(f3_name, calib_file.get());
  auto [vec_c_tp, vec_s_tp, vec_cov_tp] = ReadCnSn(tp_name, calib_file.get());
  auto [vec_c_tn, vec_s_tn, vec_cov_tn] = ReadCnSn(tn_name, calib_file.get());
  auto [vec_c_p, vec_s_p, vec_cov_p] = ReadCnSn(proton_name, calib_file.get());

  auto f1_corr = MakeCorrectionMatrix(vec_c_f1, vec_s_f1, vec_cov_f1);
  auto f2_corr = MakeCorrectionMatrix(vec_c_f2, vec_s_f2, vec_cov_f2);
  auto f3_corr = MakeCorrectionMatrix(vec_c_f3, vec_s_f3, vec_cov_f3);
  
  auto tp_corr = MakeCorrectionMatrix(vec_c_tp, vec_s_tp, vec_cov_tp);
  auto tn_corr = MakeCorrectionMatrix(vec_c_tn, vec_s_tn, vec_cov_tn);
  
  auto p_corr = MakeCorrectionMatrix(vec_c_p, vec_s_p, vec_cov_p);

  const auto correction_generator = []( const std::vector<Qn::DataContainerStatCalculate>& vec_c, const std::vector<Qn::DataContainerStatCalculate>& vec_s, const Qn::DataContainer<Eigen::Matrix4d, Qn::AxisD>& vec_cov ){
    return [&vec_c, &vec_s, &vec_cov]( Qn::DataContainerQVector qvec, Double_t centrality, Double_t run_id ) -> Qn::DataContainerQVector {
      auto new_qvec = qvec;
      auto c_axis = vec_c[0].GetAxis( "centrality" ); 
      auto c_bin = c_axis.FindBin( centrality );
      auto c_bin_lo = c_axis.GetLowerBinEdge( c_bin );
      auto c_bin_hi = c_axis.GetUpperBinEdge( c_bin );
      auto new_c_axis = Qn::AxisD{ "centrality", 1, c_bin_lo, c_bin_hi };

      auto r_axis = vec_c[0].GetAxis( "runId" ); 
      auto r_bin = r_axis.FindBin( run_id );
      auto r_bin_lo = r_axis.GetLowerBinEdge( r_bin );
      auto r_bin_hi = r_axis.GetUpperBinEdge( r_bin );
      auto new_r_axis = Qn::AxisD{ "runId", 1, r_bin_lo, r_bin_hi };

      auto vec_c_c = vec_c;
      auto vec_s_c = vec_s;

      for( auto& el : vec_c_c ) {
        el.Select( new_r_axis );
        if( el.GetAxes().size() > 1 ){
          el = el.Select( new_c_axis );
        } else {
          el = el.Rebin( new_c_axis );
        }
      }
      for( auto& el : vec_s_c ) {
        el.Select( new_r_axis );
        if( el.GetAxes().size() > 1 ){
          el = el.Select( new_c_axis );
        } else {
          el = el.Rebin( new_c_axis );
        }
      }
      
      auto vec_cov_c = vec_cov.Select( new_r_axis );
      if( vec_cov_c.GetAxes().size() > 1 ){
          vec_cov_c = vec_cov_c.Select( new_c_axis );
        } else {
          vec_cov_c = vec_cov_c.Rebin( new_c_axis, [](auto a, auto b){ return (a + b) * 0.5; } );
        }

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        auto c1 = vec_c_c[0].At(i).Mean();
        auto c2 = vec_c_c[1].At(i).Mean();
        auto s1 = vec_s_c[0].At(i).Mean();
        auto s2 = vec_s_c[1].At(i).Mean();

        auto x1_old = qvec.At(i).x(1) - c1;
        auto x2_old = qvec.At(i).x(2) - c2;
        auto y1_old = qvec.At(i).y(1) - s1;
        auto y2_old = qvec.At(i).y(2) - s2;
    
        auto Minv = vec_cov_c.At(i);
        auto Xold =  Eigen::Vector4d{ x1_old, y1_old, x2_old, y2_old };
        auto Xnew = Minv * Xold;

        auto x1_new = static_cast<double>(Xnew(0));
        auto y1_new = static_cast<double>(Xnew(1));
        auto x2_new = static_cast<double>(Xnew(2));
        auto y2_new = static_cast<double>(Xnew(3));
        
        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_new, y2_new );
      }

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Define("F1_DECOMPOSED", correction_generator(vec_c_f1, vec_s_f1, f1_corr), { f1_name, "centrality", "runId" } )
    .Define("F2_DECOMPOSED", correction_generator(vec_c_f2, vec_s_f2, f2_corr), { f2_name, "centrality", "runId" } )
    .Define("F3_DECOMPOSED", correction_generator(vec_c_f3, vec_s_f3, f3_corr), { f3_name, "centrality", "runId" } )

    .Define("Tpos_DECOMPOSED", correction_generator(vec_c_tp, vec_s_tp, tp_corr), { tp_name, "centrality", "runId" } )
    .Define("Tneg_DECOMPOSED", correction_generator(vec_c_tn, vec_s_tn, tn_corr), { tn_name, "centrality", "runId" } )
    
    .Define("proton_DECOMPOSED", correction_generator(vec_c_p, vec_s_p, p_corr), { proton_name, "centrality", "runId" } )
    .Filter( "1 < centrality && centrality < 60" )
    .Filter( "6600 < runId && runId < 8300" )
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, tp{}, tn{}, p{};
  double cent{};
  double runid{};
  tree->Branch( "centrality", &cent );
  tree->Branch( "runId", &runid );

  tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );

  dd.Foreach( [tree, &cent, &runid, &f1, &f2, &f3, &tp, &tn, &p]( 
    double centrality,
    double run_id,
    Qn::DataContainerQVector f1_ev, 
    Qn::DataContainerQVector f2_ev, 
    Qn::DataContainerQVector f3_ev, 
    Qn::DataContainerQVector tp_ev, 
    Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev ) mutable {
    cent = centrality;
    runid = runid;
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
    "runId",
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