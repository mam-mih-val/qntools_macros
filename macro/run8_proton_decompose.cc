#include "DataContainer.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"
#include "QVector.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

template<typename T>
using vector1d = std::vector<T>;

template<typename T>
using vector2d = std::vector<std::vector<T>>;

template<typename T>
using vector3d = std::vector<std::vector<std::vector<T>>>;

using DataContainerMatrix = Qn::DataContainer<Eigen::Matrix4d, Qn::AxisD>;

DataContainerMatrix MakeCorrectionMatrix(const vector1d<Qn::DataContainerStatCalculate>& vec_c, const vector1d<Qn::DataContainerStatCalculate>& vec_s, const vector1d<Qn::DataContainerStatCalculate>& vec_cov ){
  std::cout << __func__ << std::endl;
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
    auto Minv = solver.operatorInverseSqrt() * pow(2.0, - 0.5);
    corr_matrix.At(i) = Minv;
  }
  return corr_matrix;
}


std::tuple< vector1d<Qn::DataContainerStatCalculate>, vector1d<Qn::DataContainerStatCalculate>, vector1d<Qn::DataContainerStatCalculate> > ReadCnSn( std::string str_vec_name, TFile* calib_file ){
  std::cout << __func__ << std::endl;
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_s = std::vector<Qn::DataContainerStatCalculate>{};
  vec_c.reserve(4);
  vec_s.reserve(4);
  for( auto i=size_t{0}; i<4; ++i ){
    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centralityrunId"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.emplace_back( *tmp );
    std::cout << "Extracting " << corr_name << "\n";
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
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );
    }
  }
  return {vec_c, vec_s, vec_cov};
}

template<typename T>
vector2d < Qn::DataContainer<T, Qn::AxisD> >  ExtractEventAxes( const Qn::DataContainer<T, Qn::AxisD>& container, const std::vector< Qn::AxisD >& axes ){
  std::cout << __func__ << std::endl;
  vector2d < Qn::DataContainer<T, Qn::AxisD> > result;
  result.reserve( axes.at(0).size() );
  for( auto bin_1 = size_t{0}; bin_1 < axes.at(0).size(); ++bin_1 ){
    auto lo_1 = axes.at(0).GetLowerBinEdge( bin_1 );
    auto hi_1 = axes.at(0).GetUpperBinEdge( bin_1 );
    auto name_1 = axes.at(0).Name();
    auto new_axes_1 = Qn::AxisD{ name_1, 1, lo_1, hi_1 };
    auto container_a1 = container.Select( new_axes_1 );
    result.emplace_back();
    result.back().reserve( axes.at(1).size() );
    for( auto bin_2 = size_t{0}; bin_2 < axes.at(1).size(); ++bin_2 ){
      auto lo_2 = axes.at(1).GetLowerBinEdge( bin_2 );
      auto hi_2 = axes.at(1).GetUpperBinEdge( bin_2 );
      auto name_2 = axes.at(1).Name();
      auto new_axes_2 = Qn::AxisD{ name_2, 1, lo_2, hi_2 };
      auto container_a2 = container_a1.GetAxes().size() > 1 ? container_a1.Select( new_axes_2 ) : container_a1.Rebin( new_axes_2, [](const auto& a, const auto& b){ return (a+b)/2; } );
      result.back().emplace_back(container_a2);
    }
  }
  return result;
}

vector3d< Qn::DataContainerStatCalculate > ExtractPack( const vector1d<Qn::DataContainerStatCalculate>& vec_containers, const vector1d<Qn::AxisD>& axes ){
  std::cout << __func__ << std::endl;
  auto result = vector3d< Qn::DataContainerStatCalculate >{};
  result.reserve( vec_containers.size() );
  for( const auto& container : vec_containers ){
    result.emplace_back( ExtractEventAxes(container, axes) );
  }
  return result;
}

void run8_proton_decompose(std::string in_file_name, std::string in_calib_file){

  auto event_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "runId", 17, 6600, 8300 },
  };

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

  auto c_f1 = ExtractPack(vec_c_f1, event_axes);
  auto s_f1 = ExtractPack(vec_s_f1, event_axes);
  auto cov_f1 = ExtractEventAxes(f1_corr, event_axes);

  auto c_f2 = ExtractPack(vec_c_f2, event_axes);
  auto s_f2 = ExtractPack(vec_s_f2, event_axes);
  auto cov_f2 = ExtractEventAxes(f2_corr, event_axes);

  auto c_f3 = ExtractPack(vec_c_f3, event_axes);
  auto s_f3 = ExtractPack(vec_s_f3, event_axes);
  auto cov_f3 = ExtractEventAxes(f3_corr, event_axes);

  auto c_tp = ExtractPack(vec_c_tp, event_axes);
  auto s_tp = ExtractPack(vec_s_tp, event_axes);
  auto cov_tp = ExtractEventAxes(tp_corr, event_axes);

  auto c_tn = ExtractPack(vec_c_tn, event_axes);
  auto s_tn = ExtractPack(vec_s_tn, event_axes);
  auto cov_tn = ExtractEventAxes(tn_corr, event_axes);

  auto c_p = ExtractPack(vec_c_p, event_axes);
  auto s_p = ExtractPack(vec_s_p, event_axes);
  auto cov_p = ExtractEventAxes(p_corr, event_axes);

  const auto correction_generator = []( 
    const vector3d<Qn::DataContainerStatCalculate>& vec_c, 
    const vector3d<Qn::DataContainerStatCalculate>& vec_s, 
    const vector2d<DataContainerMatrix>& vec_cov,
    const vector1d<Qn::AxisD>& axes 
  ){
    return [&vec_c, &vec_s, &vec_cov, &axes]( Qn::DataContainerQVector qvec, Double_t centrality, Double_t run_id ) -> Qn::DataContainerQVector {
      auto new_qvec = qvec;
      auto c_bin = axes.at(0).FindBin( centrality ); 
      auto r_bin = axes.at(1).FindBin( run_id ); 

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        if( fabs( qvec.At(i).sumweights()) < std::numeric_limits<double>::min() )
          continue;
        auto c1 = vec_c.at(0).at(c_bin).at(r_bin).At(i).Mean();
        auto c2 = vec_c.at(1).at(c_bin).at(r_bin).At(i).Mean();
        
        auto s1 = vec_s.at(0).at(c_bin).at(r_bin).At(i).Mean();
        auto s2 = vec_s.at(1).at(c_bin).at(r_bin).At(i).Mean();

        auto x1_old = qvec.At(i).x(1) - c1;
        auto x2_old = qvec.At(i).x(2) - c2;
        auto y1_old = qvec.At(i).y(1) - s1;
        auto y2_old = qvec.At(i).y(2) - s2;
    
        auto Minv = vec_cov.at(c_bin).at(r_bin).At(i);
        
        if( std::isnan(Minv(0, 0)) ){
          new_qvec.At(i).Reset();
        }

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
    .Filter( "1 < centrality && centrality < 60" )
    .Filter( "6700 < runId && runId < 8300" )
    .Define("F1_DECOMPOSED", correction_generator(c_f1, s_f1, cov_f1, event_axes), { f1_name, "centrality", "runId" } )
    .Define("F2_DECOMPOSED", correction_generator(c_f2, s_f2, cov_f2, event_axes), { f2_name, "centrality", "runId" } )
    .Define("F3_DECOMPOSED", correction_generator(c_f3, s_f3, cov_f3, event_axes), { f3_name, "centrality", "runId" } )

    .Define("Tpos_DECOMPOSED", correction_generator(c_tn, s_tn, cov_tn, event_axes), { tp_name, "centrality", "runId" } )
    .Define("Tneg_DECOMPOSED", correction_generator(c_tp, s_tp, cov_tp, event_axes), { tn_name, "centrality", "runId" } )
    
    .Define("proton_DECOMPOSED", correction_generator(c_p, s_p, cov_p, event_axes), { proton_name, "centrality", "runId" } )
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, tp{}, tn{}, p{};
  double cent{};
  double r_id{};
  tree->Branch( "centrality", &cent );
  tree->Branch( "runId", &r_id );

  tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  dd.Foreach( [tree, &cent, &r_id, &f1, &f2, &f3, &tp, &tn, &p]( 
    double centrality,
    double run_id,
    Qn::DataContainerQVector f1_ev, 
    Qn::DataContainerQVector f2_ev, 
    Qn::DataContainerQVector f3_ev, 
    Qn::DataContainerQVector tp_ev, 
    Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev ) mutable {
    cent = centrality;
    r_id = run_id;
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

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " min" << std::endl;
  auto elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  std::cout << "It is " << elapsed_s / tree->GetEntries() << " μs/ev." << std::endl;
}