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

constexpr size_t NDIM = 2;
constexpr size_t N_HARM = 4;

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

const auto twist_rescaling_mixing_matrix = [](const vector1d<double>& vec_c, const vector1d<double>& vec_s){
  auto c1 = vec_c[0];
  auto c2 = vec_c[1];
  
  auto s1 = vec_s[0];
  auto s2 = vec_s[1];
  
  auto M = mixing_matrix_t{};

  M << 
    1+c2,    s2,
    s2,    1-c2
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
      Splus(i, i) = 1.0 / s;
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
vector2d<DataContainerMatrix> MakeCorrectionMatrix(
  const vector3d<Qn::DataContainerStatCalculate>& vec3_c, 
  const vector3d<Qn::DataContainerStatCalculate>& vec3_s,
  const Func& func, 
  size_t n_harm = 4 ){
  auto result = vector2d<DataContainerMatrix>{};
  result.reserve( n_harm );

  for( auto harm = size_t{}; harm < n_harm; ++harm ){
    result.emplace_back();
    auto vec_c = vec3_c[harm];
    auto vec_s = vec3_s[harm];
    result.back().reserve( vec_c.size() );
    for(auto ev_bin=size_t{0}; ev_bin<vec_c.front().size(); ++ev_bin ){
      auto axes = vec_c[0][ev_bin].GetAxes();
      DataContainerMatrix corr_matrix;
      corr_matrix.AddAxes(axes);

      for( auto i = size_t{0}; i<vec_c[0][ev_bin].size(); ++i ){
        auto vec_double_c = std::vector<double>(NDIM, 0);
        auto vec_double_s = std::vector<double>(NDIM, 0);

        vec_double_c[0] = vec_c[0][ev_bin].At(i).Mean();
        vec_double_c[1] = vec_c[1][ev_bin].At(i).Mean();

        vec_double_s[0] = vec_s[0][ev_bin].At(i).Mean();
        vec_double_s[1] = vec_s[1][ev_bin].At(i).Mean();

        auto sumw = vec_c[0][ev_bin].At(i).SumWeights();
        auto M = func( vec_double_c, vec_double_s );

        auto c = column_t{};
        c << vec_double_c[0], vec_double_s[0];
      
        auto [is_valid, Minv] = PseudoInverse( M, 5e-3 );
        if( std::isinf( 1.0 / sqrt(sumw) ) )
          is_valid = false;
        std::cout << " 1 / sqrt(sumw) = " << 1.0 / sqrt(sumw) << "\n";

        corr_matrix.At(i) = std::tuple{is_valid, Minv, c};
      }
      result.back().push_back(corr_matrix);
    }
  }
  return result;
}


std::tuple< vector2d<Qn::DataContainerStatCalculate>, vector2d<Qn::DataContainerStatCalculate> > ReadCnSn( 
  std::string str_vec_name, TFile* calib_file, const size_t n_harm = 4
){
  std::cout << __func__ << std::endl;
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = vector2d<Qn::DataContainerStatCalculate>{};
  auto vec_s = vector2d<Qn::DataContainerStatCalculate>{};
  vec_c.reserve( n_harm );
  vec_s.reserve( n_harm );
  for( auto i=size_t{0}; i < n_harm; ++i ){
    vec_c.emplace_back(); vec_c.back().reserve( 2 );
    vec_s.emplace_back(); vec_s.back().reserve( 2 );

    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.back().emplace_back( *tmp );
    
    corr_name = str_vec_name+".y"+std::to_string(i+1)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.back().emplace_back( *tmp );

    corr_name = str_vec_name+".x"+std::to_string(i+1)+std::to_string((i+1)*2)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.back().emplace_back( *tmp );

    corr_name = str_vec_name+".y"+std::to_string(i+1)+std::to_string((i+1)*2)+"centrality"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.back().emplace_back( *tmp );
  }
  return {vec_c, vec_s};
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
vector2d< Qn::DataContainer<T, Qn::AxisD> > ExtractVector( const vector1d< Qn::DataContainer<T, Qn::AxisD> >& vec_containers, const Linearization& lin ){
  std::cout << __func__ << std::endl;
  auto result = vector2d< Qn::DataContainer<T, Qn::AxisD> >{};
  result.reserve( vec_containers.size() );
  for( const auto& container : vec_containers ){
    result.emplace_back( ExtractEventAxes(container, lin) );
  }
  return result;
}

template<typename T>
vector3d< Qn::DataContainer<T, Qn::AxisD> > ExtractPack( const vector2d< Qn::DataContainer<T, Qn::AxisD> >& vec2_containers, const Linearization& lin ){
  std::cout << __func__ << std::endl;
  auto result = vector3d< Qn::DataContainer<T, Qn::AxisD> >{};
  result.reserve( vec2_containers.size() );
  for( auto harm = size_t{}; harm < vec2_containers.size(); ++harm ){
    result.emplace_back( ExtractVector( vec2_containers.at(harm), lin ) );
  }
  return result;
}

void run8_mc_proton_twist_rescale(std::string in_file_name, std::string in_calib_file){

  auto event_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
  };

  const std::string f1_name {"F1_RECENTERED"};
  const std::string f2_name {"F2_RECENTERED"};
  const std::string f3_name {"F3_RECENTERED"};
  const std::string f4_name {"F4_RECENTERED"};
  const std::string tp_name {"Tpos_RECENTERED"};
  const std::string tn_name {"Tneg_RECENTERED"};
  const std::string proton_name {"proton_RECENTERED"};

  auto calib_file = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( in_calib_file.c_str(), "READ" ), [](auto f){f ->Close(); } };
  assert(calib_file);

  auto [vec2_c_f1, vec2_s_f1] = ReadCnSn(f1_name, calib_file.get(), 2);
  auto [vec2_c_f2, vec2_s_f2] = ReadCnSn(f2_name, calib_file.get(), 2);
  auto [vec2_c_f3, vec2_s_f3] = ReadCnSn(f3_name, calib_file.get(), 2);
  auto [vec2_c_f4, vec2_s_f4] = ReadCnSn(f4_name, calib_file.get(), 2);
  auto [vec2_c_tp, vec2_s_tp] = ReadCnSn(tp_name, calib_file.get(), 2);
  auto [vec2_c_tn, vec2_s_tn] = ReadCnSn(tn_name, calib_file.get(), 2);
  auto [vec2_c_p, vec2_s_p] = ReadCnSn(proton_name, calib_file.get(), 2);

  auto lin = Linearization( event_axes );

  auto v2_c_f1 = ExtractPack( vec2_c_f1, lin );
  auto v2_s_f1 = ExtractPack( vec2_s_f1, lin );

  auto v2_c_f2 = ExtractPack( vec2_c_f2, lin );
  auto v2_s_f2 = ExtractPack( vec2_s_f2, lin );

  auto v2_c_f3 = ExtractPack( vec2_c_f3, lin );
  auto v2_s_f3 = ExtractPack( vec2_s_f3, lin );

  auto v2_c_f4 = ExtractPack( vec2_c_f4, lin );
  auto v2_s_f4 = ExtractPack( vec2_s_f4, lin );

  auto v2_c_tp = ExtractPack( vec2_c_tp, lin );
  auto v2_s_tp = ExtractPack( vec2_s_tp, lin );

  auto v2_c_tn = ExtractPack( vec2_c_tn, lin );
  auto v2_s_tn = ExtractPack( vec2_s_tn, lin );

  auto v2_c_p = ExtractPack( vec2_c_p, lin );
  auto v2_s_p = ExtractPack( vec2_s_p, lin );
  

  auto f1_corr = MakeCorrectionMatrix(v2_c_f1, v2_s_f1, twist_rescaling_mixing_matrix, 2);
  auto f2_corr = MakeCorrectionMatrix(v2_c_f2, v2_s_f2, twist_rescaling_mixing_matrix, 2);
  auto f3_corr = MakeCorrectionMatrix(v2_c_f3, v2_s_f3, twist_rescaling_mixing_matrix, 2);
  auto f4_corr = MakeCorrectionMatrix(v2_c_f4, v2_s_f4, twist_rescaling_mixing_matrix, 2);

  auto tp_corr = MakeCorrectionMatrix(v2_c_tp, v2_s_tp, twist_rescaling_mixing_matrix, 2);
  auto tn_corr = MakeCorrectionMatrix(v2_c_tn, v2_s_tn, twist_rescaling_mixing_matrix, 2);
  
  auto p_corr = MakeCorrectionMatrix(v2_c_p, v2_s_p, twist_rescaling_mixing_matrix, 2 );

  const auto correction_generator = []( 
    const vector2d<DataContainerMatrix>& vec_cor,
    const Linearization& lin,
    const size_t n_harm = 2
  ){
    return [&vec_cor, &lin, n_harm]( Qn::DataContainerQVector qvec, Double_t centrality ) -> Qn::DataContainerQVector {
      auto new_qvec = qvec;
      auto l_idx = lin[ std::vector<double>{ centrality } ];
      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        if( fabs( qvec.At(i).sumweights()) < std::numeric_limits<double>::min() )
          continue;
        for( auto harm = size_t{1}; harm <= n_harm; ++harm ){
          // std::cout << "Here: 1" << "\n";
          auto x_old = qvec.At(i).x(harm);
          auto y_old = qvec.At(i).y(harm);
          // std::cout << "Here: 2" << "\n";
          auto [is_valid, Minv, c] = vec_cor.at(harm-1).at(l_idx).At(i);
          // std::cout << "Here: 3" << "\n";
          if( !is_valid ){
            new_qvec.At(i).Reset();
            continue;
          }
          auto Xold =  column_t{};
          Xold << x_old, y_old;
          auto Xnew =  Minv * ( Xold - c );
          // std::cout << "Here: 4" << "\n";
          auto x_new = static_cast<double>(Xnew(0));
          auto y_new = static_cast<double>(Xnew(1));
          new_qvec.At(i).SetQ( harm, x_new, y_new );
        }
      }
      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Filter( "1 < centrality && centrality < 60" )
    .Define("F1_DECOMPOSED", correction_generator(f1_corr, lin, 2), { f1_name, "centrality" })
    .Define("F2_DECOMPOSED", correction_generator(f2_corr, lin, 2), { f2_name, "centrality" })
    .Define("F3_DECOMPOSED", correction_generator(f3_corr, lin, 2), { f3_name, "centrality" })
    .Define("F4_DECOMPOSED", correction_generator(f4_corr, lin, 2), { f4_name, "centrality" })

    .Define("Tpos_DECOMPOSED", correction_generator(tn_corr, lin, 2), { tp_name, "centrality" })
    .Define("Tneg_DECOMPOSED", correction_generator(tp_corr, lin, 2), { tn_name, "centrality" })
    
    .Define("proton_DECOMPOSED", correction_generator(p_corr, lin, 2), { proton_name, "centrality" })
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, f4{}, tp{}, tn{}, p{}, tru_p{}, psi_rp{};
  double cent{};
  
  tree->Branch( "centrality", &cent );

  tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  tree->Branch( "F4_DECOMPOSED", "Qn::DataContainerQVector", &f4 );
  tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );
  tree->Branch( "tru_proton_PLAIN", "Qn::DataContainerQVector", &tru_p );
  tree->Branch( "psi_rp_PLAIN", "Qn::DataContainerQVector", &psi_rp );

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  dd.Foreach( [tree, &cent, &f1, &f2, &f3, &f4, &tp, &tn, &p, &tru_p, &psi_rp]( 
    double centrality,
    Qn::DataContainerQVector f1_ev, 
    Qn::DataContainerQVector f2_ev, 
    Qn::DataContainerQVector f3_ev, 
    Qn::DataContainerQVector f4_ev, 
    Qn::DataContainerQVector tp_ev, 
    Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev,
    Qn::DataContainerQVector tru_p_ev,
    Qn::DataContainerQVector psi_rp_ev
   ) mutable {
    cent = centrality;
    f1 = f1_ev;
    f2 = f2_ev;
    f3 = f3_ev;
    f4 = f4_ev;
    tp = tp_ev;
    tn = tn_ev;
    p = p_ev;
    tru_p = tru_p_ev;
    psi_rp = psi_rp_ev;
    tree->Fill();
  }, 
  std::vector<std::string>{ 
    "centrality",
    "F1_DECOMPOSED",
    "F2_DECOMPOSED",
    "F3_DECOMPOSED",
    "F4_DECOMPOSED",
    "Tpos_DECOMPOSED",
    "Tneg_DECOMPOSED",
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