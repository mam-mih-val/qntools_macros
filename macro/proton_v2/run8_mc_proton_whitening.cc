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
#include <tuple>
#include <vector>

constexpr size_t NDIM = 4;

using correction_matrix_t = Eigen::Matrix<double, NDIM, NDIM>;
using mixing_matrix_t = Eigen::Matrix<double, NDIM, NDIM>;
using column_t = Eigen::Matrix<double, NDIM, 1>;

template<typename T>
using vector1d = std::vector<T>;

template<typename T>
using vector2d = std::vector<std::vector<T>>;

template<typename T>
using vector3d = std::vector<std::vector<std::vector<T>>>;

using DataContainerMatrix = Qn::DataContainer< std::tuple<bool, correction_matrix_t, column_t>, Qn::AxisD>;

class Linearization{
public:
  Linearization(  vector1d<Qn::AxisD> a_vec ) : axis_vector{ std::move(a_vec) }, offset_vector( axis_vector.size() ) {
    std::fill(offset_vector.begin(), offset_vector.end(), 1);
    for( int i = offset_vector.size() - 2; i >= 0; --i ){
      offset_vector.at(i) = offset_vector.at(i+1);
      offset_vector.at(i) *= axis_vector.at(i+1).size();
    }
  }
  inline auto operator[]( const vector1d<size_t>& coordinates ) const -> size_t {
    assert( coordinates.size() == axis_vector.size() );
    auto l_idx = size_t{0};
    for( auto i=0; i<offset_vector.size(); ++i ){
      l_idx += coordinates[i] * offset_vector[i];
    }
    return l_idx;
  }
  inline auto operator[]( const vector1d<double>& values ) const -> size_t {
    assert( values.size() == axis_vector.size() );
    auto vec_idx = vector1d<size_t>( values.size() );
    for( auto i=0; i<axis_vector.size(); ++i ){
      vec_idx[i] = axis_vector[i].FindBin( values[i] );
    }
    return this->operator[]( vec_idx );
  }
  inline auto operator()( size_t l_idx ) const -> vector1d<size_t> {
    auto v_coord = vector1d<size_t>( offset_vector.size(), size_t{0} );
    for( auto i=0; i<offset_vector.size(); ++i ){
      v_coord[i] = l_idx / offset_vector[i];
      l_idx -= v_coord[i] * offset_vector[i];
    }
    return v_coord;
  }
  inline auto size() const -> size_t {
    return axis_vector.front().size() * offset_vector.front();
  }
  inline auto SelectionAxes( const vector1d<size_t>& coordinates ) const -> vector1d<Qn::AxisD>{
    auto axes = vector1d<Qn::AxisD>{};
    axes.reserve( axis_vector.size() );
    for( auto i = size_t{}; i<axis_vector.size(); ++i ){
      auto name = axis_vector[i].Name();
      auto lo = axis_vector[i].GetLowerBinEdge( coordinates[i] );
      auto hi = axis_vector[i].GetUpperBinEdge( coordinates[i] );
      axes.emplace_back( name, 1, lo, hi );
    }
    return axes;
  }
  inline auto SelectionAxes( size_t l_idx ) const -> vector1d<Qn::AxisD>{
    auto coordinates = this->operator()( l_idx );
    return SelectionAxes( coordinates );
  }
  auto Axes() const -> const vector1d<Qn::AxisD>& { return axis_vector; }

private:
  vector1d<Qn::AxisD> axis_vector{};
  vector1d<size_t> offset_vector{};
};

const auto whitening_mixing_matrix = [](const vector1d<double>& vec_mean, const vector1d<double>& vec_cov){
  auto x1 = vec_mean[0];
  auto y1 = vec_mean[1];
  auto x2 = vec_mean[2];
  auto y2 = vec_mean[3];

  auto x1x1 = vec_cov[0] - x1*x1;
  auto x1y1 = vec_cov[1] - x1*y1;
  auto y1x1 = vec_cov[2] - x1*y1;
  auto y1y1 = vec_cov[3] - y1*y1;

  auto x1x2 = vec_cov[4] - x1*x2;
  auto x1y2 = vec_cov[5] - x1*y2;
  auto y1x2 = vec_cov[6] - y1*x2;
  auto y1y2 = vec_cov[7] - y1*y2;

  auto x2x2 = vec_cov[8] - x2*x2;
  auto x2y2 = vec_cov[9] - x2*y2;
  auto y2x2 = vec_cov[10] - y2*x2;
  auto y2y2 = vec_cov[11] - y2*y2;
  
  auto M = mixing_matrix_t{};

  M << 
    //       x1    y1    x2    y2
  /*x1*/    x1x1, y1x1, x1x2, x1y2,
  /*y1*/    y1x1, y1y1, y1x2, y1y2,
  /*x2*/    x1x2, y1x2, x2x2, x2y2,
  /*y2*/    x1y2, y1y2, x2y2, y2y2
  ;

  return M;
};

std::tuple<bool, correction_matrix_t> PseudoInverse( const correction_matrix_t& M, double l ){
  auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
  auto singular_values = svd.singularValues();
  auto U = svd.matrixU();
  auto V = svd.matrixV();
  auto Splus = correction_matrix_t{};
  auto rank = size_t{0};
  for (auto i = size_t{0}; i < singular_values.size(); ++i) {
    auto s = singular_values(i);
    if( fabs(s) > l ){
      Splus(i, i) = sqrt(2.0 / s);
      rank++;
    }
  }
  auto is_valid = rank == NDIM;
  auto E = correction_matrix_t::Identity();
  auto Vr = V.leftCols(rank);
  auto Etilda = Vr * Vr.transpose() * E;
  auto Mpinv = correction_matrix_t{ V * Splus * U.transpose() };
  std::cout << "l: " << l << "\nMatrix M:\n" << M << "\nMatrix U:\n" << U << "\nS: " << singular_values.transpose() << "\nMatrix V:\n" << V << "\nInverse:\n" << Mpinv << "\nEtilda:\n" << Etilda << "\n";
  return {is_valid, Mpinv};
}


template<typename Func>
vector1d<DataContainerMatrix> MakeCorrectionMatrix(
  const vector2d<Qn::DataContainerStatCalculate>& vec_mean, 
  const vector2d<Qn::DataContainerStatCalculate>& vec_cov,
  const Func& func ){
  auto result = vector1d<DataContainerMatrix>{};
  result.reserve( vec_mean.front().size() );

  for(auto ev_bin=size_t{0}; ev_bin<vec_mean.front().size(); ++ev_bin ){
    auto axes = vec_mean[0][ev_bin].GetAxes();
    DataContainerMatrix corr_matrix;
    corr_matrix.AddAxes(axes);

    for( auto i = size_t{0}; i<vec_mean[0][ev_bin].size(); ++i ){
      auto vec_double_mean = std::vector<double>{}; vec_double_mean.reserve( vec_mean.size() );
      auto vec_double_cov = std::vector<double>{}; vec_double_cov.reserve( vec_cov.size() );

      std::for_each( vec_mean.begin(), vec_mean.end(), [&vec_double_mean, i, ev_bin]( const auto& corr ) mutable { vec_double_mean.push_back( corr[ev_bin].At(i).Mean()); } );
      std::for_each( vec_cov.begin(), vec_cov.end(), [&vec_double_cov, i, ev_bin]( const auto& corr ) mutable { vec_double_cov.push_back( corr[ev_bin].At(i).Mean()); } );

      auto sumw = vec_mean[0][ev_bin].At(i).SumWeights();
      auto M = func( vec_double_mean, vec_double_cov );
    
      auto c = column_t{};
      c << vec_double_mean[0], vec_double_mean[1], vec_double_mean[2], vec_double_mean[3];

      auto [is_valid, Minv] = PseudoInverse( M, 5e-3 );
      if( std::isinf( 1.0 / sqrt(sumw) ) )
        is_valid = false;
      std::cout << " 1 / sqrt(sumw) = " << 1.0 / sqrt(sumw) << "\n";

      corr_matrix.At(i) = std::tuple{is_valid, Minv, c};
    }
    result.push_back(corr_matrix);
  }
  return result;
}


std::tuple< vector1d<Qn::DataContainerStatCalculate>, vector1d<Qn::DataContainerStatCalculate> > ReadCnSn( std::string str_vec_name, TFile* calib_file ){
  std::cout << __func__ << std::endl;
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_mean = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_cov = std::vector<Qn::DataContainerStatCalculate>{};
  vec_mean.reserve(2);
  vec_cov.reserve(3);

  for( auto h_a = size_t{1}; h_a <= 4; ++h_a ){
    auto corr_name = str_vec_name+".x"+std::to_string(h_a)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_mean.emplace_back( *tmp );
    
    corr_name = str_vec_name+".y"+std::to_string(h_a)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_mean.emplace_back( *tmp );

    for( auto h_b = h_a; h_b <= 4; ++h_b ){
      corr_name = str_vec_name+".x"+std::to_string(h_a)+"x"s+std::to_string(h_b)+"centrality"s;
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );
      
      corr_name = str_vec_name+".x"+std::to_string(h_a)+"y"s+std::to_string(h_b)+"centrality"s;
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );

      corr_name = str_vec_name+".y"+std::to_string(h_a)+"x"s+std::to_string(h_b)+"centrality"s;
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );

      corr_name = str_vec_name+".y"+std::to_string(h_a)+"y"s+std::to_string(h_b)+"centrality"s;
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp ); 
    }
  }
  return {vec_mean, vec_cov};
}

template<typename T>
vector1d < Qn::DataContainer<T, Qn::AxisD> > ExtractEventAxes( const Qn::DataContainer<T, Qn::AxisD>& container, const Linearization& lin ){
  std::cout << __func__ << std::endl;
  vector1d < Qn::DataContainer<T, Qn::AxisD> > result;
  result.reserve( lin.size() );
  auto rebinned_container = container;
  for( const auto& a : lin.Axes() ){
    rebinned_container = rebinned_container.Rebin( a );
  }
  for( auto i=size_t{0}; i<lin.size(); ++i ){
    auto axes = lin.SelectionAxes( i );
    auto dc = rebinned_container;
    for( const auto& a : axes ){
      dc = dc.GetAxes().size() > 1 ? dc.Select( a ) : dc.Rebin( a );
    }
    result.push_back(dc);
  }
  return result;
}

template<typename T>
vector2d< Qn::DataContainer<T, Qn::AxisD> > ExtractPack( const vector1d< Qn::DataContainer<T, Qn::AxisD> >& vec_containers, const Linearization& lin ){
  std::cout << __func__ << std::endl;
  auto result = vector2d< Qn::DataContainer<T, Qn::AxisD> >{};
  result.reserve( vec_containers.size() );
  for( const auto& container : vec_containers ){
    result.emplace_back( ExtractEventAxes(container, lin) );
  }
  return result;
}

void run8_mc_proton_whitening(std::string in_file_name, std::string in_calib_file){

  auto event_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
  };

  const std::string f1_name {"F1_PLAIN"};
  const std::string f2_name {"F2_PLAIN"};
  const std::string f3_name {"F3_PLAIN"};
  const std::string f4_name {"F4_PLAIN"};
  const std::string tp_name {"Tpos_PLAIN"};
  const std::string tn_name {"Tneg_PLAIN"};
  const std::string proton_name {"proton_PLAIN"};

  auto calib_file = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( in_calib_file.c_str(), "READ" ), [](auto f){f ->Close(); } };
  assert(calib_file);

  // auto [vec_mean_f1, vec_cov_f1] = ReadCnSn(f1_name, calib_file.get());
  // auto [vec_mean_f2, vec_cov_f2] = ReadCnSn(f2_name, calib_file.get());
  // auto [vec_mean_f3, vec_cov_f3] = ReadCnSn(f3_name, calib_file.get());
  // auto [vec_mean_f4, vec_cov_f4] = ReadCnSn(f4_name, calib_file.get());
  // auto [vec_mean_tp, vec_cov_tp] = ReadCnSn(tp_name, calib_file.get());
  // auto [vec_mean_tn, vec_cov_tn] = ReadCnSn(tn_name, calib_file.get());
  auto [vec_mean_p, vec_cov_p] = ReadCnSn(proton_name, calib_file.get());

  auto lin = Linearization( event_axes );
  // auto v2_mean_f1 = ExtractPack( vec_mean_f1, lin );
  // auto v2_cov_f1 = ExtractPack( vec_cov_f1, lin );

  // auto v2_mean_f2 = ExtractPack( vec_mean_f2, lin );
  // auto v2_cov_f2 = ExtractPack( vec_cov_f2, lin );

  // auto v2_mean_f3 = ExtractPack( vec_mean_f3, lin );
  // auto v2_cov_f3 = ExtractPack( vec_cov_f3, lin );

  // auto v2_mean_f4 = ExtractPack( vec_mean_f4, lin );
  // auto v2_cov_f4 = ExtractPack( vec_cov_f4, lin );

  // auto v2_mean_tp = ExtractPack( vec_mean_tp, lin );
  // auto v2_cov_tp = ExtractPack( vec_cov_tp, lin );

  // auto v2_mean_tn = ExtractPack( vec_mean_tn, lin );
  // auto v2_cov_tn = ExtractPack( vec_cov_tn, lin );

  auto v2_mean_p = ExtractPack( vec_mean_p, lin );
  auto v2_cov_p = ExtractPack( vec_cov_p, lin );
  

  // auto f1_corr = MakeCorrectionMatrix(v2_mean_f1, v2_cov_f1, decomposition_mixing_matrix );
  // auto f2_corr = MakeCorrectionMatrix(v2_mean_f2, v2_cov_f2, decomposition_mixing_matrix);
  // auto f3_corr = MakeCorrectionMatrix(v2_mean_f3, v2_cov_f3, decomposition_mixing_matrix);
  // auto f4_corr = MakeCorrectionMatrix(v2_mean_f4, v2_cov_f4, decomposition_mixing_matrix);

  // auto tp_corr = MakeCorrectionMatrix(v2_mean_tp, v2_cov_tp, decomposition_mixing_matrix);
  // auto tn_corr = MakeCorrectionMatrix(v2_mean_tn, v2_cov_tn, decomposition_mixing_matrix);
  
  auto p_corr = MakeCorrectionMatrix(v2_mean_p, v2_cov_p, whitening_mixing_matrix );

  const auto correction_generator = []( 
    const vector1d<DataContainerMatrix>& vec_cor,
    const Linearization& lin 
  ){
    return [&vec_cor, &lin]( Qn::DataContainerQVector qvec, Double_t centrality ) -> Qn::DataContainerQVector {
      auto new_qvec = qvec;
      auto l_idx = lin[ std::vector<double>{ centrality } ];

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        if( fabs( qvec.At(i).sumweights()) < std::numeric_limits<double>::min() )
          continue;
        auto x1_old = qvec.At(i).x(1);
        auto y1_old = qvec.At(i).y(1);
    
        auto x2_old = qvec.At(i).x(2);
        auto y2_old = qvec.At(i).y(2);

        auto x3_old = qvec.At(i).x(3);
        auto y3_old = qvec.At(i).y(3);

        auto [is_valid, Minv, c] = vec_cor.at(l_idx).At(i);
        
        if( std::isnan(Minv(0, 0)) ){
          new_qvec.At(i).Reset();
          continue;
        }

        if( !is_valid ){
          new_qvec.At(i).Reset();
          continue;
        }

        auto X1old =  column_t{};
        X1old << x1_old, y1_old, x2_old, y2_old;
        // x3_old, y3_old;
        
        auto X1new =  Minv * ( X1old - c );
        
        auto x1_new = static_cast<double>(X1new(0));
        auto y1_new = static_cast<double>(X1new(1));
        
        auto x2_new = static_cast<double>(X1new(2));
        auto y2_new = static_cast<double>(X1new(3));

        // auto x3_new = static_cast<double>(X1new(4));
        // auto y3_new = static_cast<double>(X1new(5));

        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_new, y2_new );
        // new_qvec.At(i).SetQ( 3, x3_new, y3_new );
      }

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Filter( "1 < centrality && centrality < 60" )
    // .Define("F1_DECOMPOSED", correction_generator(f1_corr, lin), { f1_name, "centrality" } )
    // .Define("F2_DECOMPOSED", correction_generator(f2_corr, lin), { f2_name, "centrality" } )
    // .Define("F3_DECOMPOSED", correction_generator(f3_corr, lin), { f3_name, "centrality" } )
    // .Define("F4_DECOMPOSED", correction_generator(f4_corr, lin), { f4_name, "centrality" } )

    // .Define("Tpos_DECOMPOSED", correction_generator(tn_corr, lin), { tp_name, "centrality" } )
    // .Define("Tneg_DECOMPOSED", correction_generator(tp_corr, lin), { tn_name, "centrality" } )
    
    .Define("proton_DECOMPOSED", correction_generator(p_corr, lin), { proton_name, "centrality" } )
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, f4{}, tp{}, tn{}, p{}, tru_p{}, psi_rp{};
  double cent{};
  
  tree->Branch( "centrality", &cent );

  // tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  // tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  // tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  // tree->Branch( "F4_DECOMPOSED", "Qn::DataContainerQVector", &f4 );
  // tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  // tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );
  tree->Branch( "tru_proton_PLAIN", "Qn::DataContainerQVector", &tru_p );
  tree->Branch( "psi_rp_PLAIN", "Qn::DataContainerQVector", &psi_rp );

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  dd.Foreach( [tree, &cent, /* &f1, &f2, &f3, &f4, &tp, &tn, */ &p, &tru_p, &psi_rp]( 
    double centrality,
    // Qn::DataContainerQVector f1_ev, 
    // Qn::DataContainerQVector f2_ev, 
    // Qn::DataContainerQVector f3_ev, 
    // Qn::DataContainerQVector f4_ev, 
    // Qn::DataContainerQVector tp_ev, 
    // Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev,
    Qn::DataContainerQVector tru_p_ev,
    Qn::DataContainerQVector psi_rp_ev
   ) mutable {
    cent = centrality;
    // f1 = f1_ev;
    // f2 = f2_ev;
    // f3 = f3_ev;
    // f4 = f4_ev;
    // tp = tp_ev;
    // tn = tn_ev;
    p = p_ev;
    tru_p = tru_p_ev;
    psi_rp = psi_rp_ev;
    tree->Fill();
  }, 
  std::vector<std::string>{ 
    "centrality",
    // "F1_DECOMPOSED",
    // "F2_DECOMPOSED",
    // "F3_DECOMPOSED",
    // "F4_DECOMPOSED",
    // "Tpos_DECOMPOSED",
    // "Tneg_DECOMPOSED",
    "proton_DECOMPOSED",
    "tru_proton_PLAIN",
    "psi_rp_PLAIN",
  } );
  file_out->cd();
  tree->Write();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " min" << std::endl;
  auto elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  std::cout << "It is " << ( tree->GetEntries() > 0 ? elapsed_s / tree->GetEntries() : 0.0 ) << " μs/ev." << std::endl;
}