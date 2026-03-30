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

constexpr size_t NDIM = 13;

using correction_matrix_t = Eigen::Matrix<double, NDIM, NDIM>;
using mixing_matrix_t = Eigen::Matrix<double, NDIM, NDIM>;
using column_t = Eigen::Matrix<double, NDIM, 1>;

template<typename T>
using vector1d = std::vector<T>;

template<typename T>
using vector2d = std::vector<std::vector<T>>;

template<typename T>
using vector3d = std::vector<std::vector<std::vector<T>>>;

using DataContainerMatrix = Qn::DataContainer<column_t, Qn::AxisD>;

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

correction_matrix_t PseudoInverse( const correction_matrix_t& M, double l ){
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
  auto E = correction_matrix_t::Identity();
  auto Vr = V.leftCols(rank);
  auto Etilda = Vr * Vr.transpose() * E;
  auto Mpinv = correction_matrix_t{ V * Splus * U.transpose() };
  std::cout << "l: " << l << "\nMatrix M:\n" << M << "\nMatrix U:\n" << U << "\nS: " << singular_values.transpose() << "\nMatrix V:\n" << V << "\nInverse:\n" << Mpinv << "\nEtilda:\n" << Etilda << "\n";
  return Mpinv;
}


vector1d<DataContainerMatrix> MakeCorrectionMatrix(const vector2d<Qn::DataContainerStatCalculate>& vec_c, const vector2d<Qn::DataContainerStatCalculate>& vec_s ){
  auto result = vector1d<DataContainerMatrix>{};
  result.reserve( vec_c.front().size() );

  for(auto ev_bin=size_t{0}; ev_bin<vec_c.front().size(); ++ev_bin ){
    auto axes = vec_c[0][ev_bin].GetAxes();
    DataContainerMatrix corr_matrix;
    corr_matrix.AddAxes(axes);

    for( auto i = size_t{0}; i<vec_c[0][ev_bin].size(); ++i ){    
      auto c1 = vec_c[0][ev_bin].At(i).Mean();
      auto c2 = vec_c[1][ev_bin].At(i).Mean();
      auto c3 = vec_c[2][ev_bin].At(i).Mean();
      auto c4 = vec_c[3][ev_bin].At(i).Mean();
      auto c5 = vec_c[4][ev_bin].At(i).Mean();
      auto c6 = vec_c[5][ev_bin].At(i).Mean();
      // auto c7 = vec_c[6][ev_bin].At(i).Mean();
      // auto c8 = vec_c[7][ev_bin].At(i).Mean();

      auto s1 = vec_s[0][ev_bin].At(i).Mean();
      auto s2 = vec_s[1][ev_bin].At(i).Mean();
      auto s3 = vec_s[2][ev_bin].At(i).Mean();
      auto s4 = vec_s[3][ev_bin].At(i).Mean();
      auto s5 = vec_s[4][ev_bin].At(i).Mean();
      auto s6 = vec_s[5][ev_bin].At(i).Mean();
      // auto s7 = vec_s[6][ev_bin].At(i).Mean();
      // auto s8 = vec_s[7][ev_bin].At(i).Mean();

      auto sumw = vec_c[0][ev_bin].At(i).SumWeights();
      // auto M = mixing_matrix_t{ 
      //   { 1+c2,    s2,   c3+c1,  s3+s1 },
      //   { s2,    1-c2,   s3-s1,  c1-c3 },
      //   { c3+c1, s3-s1,  1+c4,      s4 },
      //   { s3+s1, c1-c3,    s4,    1-c4 }
      // };
      auto c = column_t{ 0.5, c1, s1, c2, s2, c3, s3, c4, s4, c5, s5, c6, s6 };
      c = c*2.0;

      // auto MTM = 2 * M.transpose() * M;
      // std::cout << "MTM\n"  << MTM << "\n\n";
      // auto C = correction_matrix_t{};
      // C.topLeftCorner(NDIM, NDIM) = MTM;
      // C(NDIM, 0) = 1;
      // C(0, NDIM) = 1;

      // auto Minv = PseudoInverse( M, 10.0 / sqrt(sumw) );      
      std::cout << " 1 / sqrt(sumw) = " << 1.0 / sqrt(sumw) << "\n";

      // corr_matrix.At(i).first = Minv;
      corr_matrix.At(i) = c;
    }
    result.push_back(corr_matrix);
  }
  return result;
}


std::tuple< vector1d<Qn::DataContainerStatCalculate>, vector1d<Qn::DataContainerStatCalculate> > ReadCnSn( std::string str_vec_name, TFile* calib_file ){
  std::cout << __func__ << std::endl;
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_s = std::vector<Qn::DataContainerStatCalculate>{};
  vec_c.reserve(6);
  vec_s.reserve(6);
  for( auto i=size_t{0}; i<5; ++i ){
    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centralityrunId"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.emplace_back( *tmp );
    
    corr_name = str_vec_name+".y"+std::to_string(i+1)+"centralityrunId"s;
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.emplace_back( *tmp );
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
vector2d< Qn::DataContainer<T, Qn::AxisD> > ExtractPack( const vector1d< Qn::DataContainer<T, Qn::AxisD> >& vec_containers, const Linearization& lin ){
  std::cout << __func__ << std::endl;
  auto result = vector2d< Qn::DataContainer<T, Qn::AxisD> >{};
  result.reserve( vec_containers.size() );
  for( const auto& container : vec_containers ){
    result.emplace_back( ExtractEventAxes(container, lin) );
  }
  return result;
}

void run8_proton_decompose(std::string in_file_name, std::string in_calib_file){

  auto event_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    // Qn::AxisD{ "vtxX", 5, -1.0, 1.0 },
    // Qn::AxisD{ "vtxY", 5, -1.0, 1.0 },
    Qn::AxisD{ "runId", 60, 7100, 8300 },
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

  auto [vec_c_f1, vec_s_f1] = ReadCnSn(f1_name, calib_file.get());
  auto [vec_c_f2, vec_s_f2] = ReadCnSn(f2_name, calib_file.get());
  auto [vec_c_f3, vec_s_f3] = ReadCnSn(f3_name, calib_file.get());
  auto [vec_c_f4, vec_s_f4] = ReadCnSn(f4_name, calib_file.get());
  auto [vec_c_tp, vec_s_tp] = ReadCnSn(tp_name, calib_file.get());
  auto [vec_c_tn, vec_s_tn] = ReadCnSn(tn_name, calib_file.get());
  auto [vec_c_p, vec_s_p] = ReadCnSn(proton_name, calib_file.get());

  auto lin = Linearization( event_axes );
  auto v2_c_f1 = ExtractPack( vec_c_f1, lin );
  auto v2_s_f1 = ExtractPack( vec_s_f1, lin );

  auto v2_c_f2 = ExtractPack( vec_c_f2, lin );
  auto v2_s_f2 = ExtractPack( vec_s_f2, lin );

  auto v2_c_f3 = ExtractPack( vec_c_f3, lin );
  auto v2_s_f3 = ExtractPack( vec_s_f3, lin );

  auto v2_c_f4 = ExtractPack( vec_c_f4, lin );
  auto v2_s_f4 = ExtractPack( vec_s_f4, lin );

  auto v2_c_tp = ExtractPack( vec_c_tp, lin );
  auto v2_s_tp = ExtractPack( vec_s_tp, lin );

  auto v2_c_tn = ExtractPack( vec_c_tn, lin );
  auto v2_s_tn = ExtractPack( vec_s_tn, lin );

  auto v2_c_p = ExtractPack( vec_c_p, lin );
  auto v2_s_p = ExtractPack( vec_s_p, lin );
  

  auto f1_corr = MakeCorrectionMatrix(v2_c_f1, v2_s_f1);
  auto f2_corr = MakeCorrectionMatrix(v2_c_f2, v2_s_f2);
  auto f3_corr = MakeCorrectionMatrix(v2_c_f3, v2_s_f3);
  auto f4_corr = MakeCorrectionMatrix(v2_c_f4, v2_s_f4);

  auto tp_corr = MakeCorrectionMatrix(v2_c_tp, v2_s_tp);
  auto tn_corr = MakeCorrectionMatrix(v2_c_tn, v2_s_tn);
  
  auto p_corr = MakeCorrectionMatrix(v2_c_p, v2_s_p);

  const auto correction_generator = []( 
    const vector1d<DataContainerMatrix>& vec_cor,
    const Linearization& lin 
  ){
    return [&vec_cor, &lin]( Qn::DataContainerQVector qvec, Double_t centrality, Double_t run_id ) -> Qn::DataContainerQVector {
      auto new_qvec = qvec;
      auto l_idx = lin[ std::vector<double>{ centrality, run_id } ];

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        if( fabs( qvec.At(i).sumweights()) < std::numeric_limits<double>::min() )
          continue;
        auto x1_old = qvec.At(i).x(1);
        auto y1_old = qvec.At(i).y(1);
    
        auto x2_old = qvec.At(i).x(2);
        auto y2_old = qvec.At(i).y(2);

        auto x3_old = qvec.At(i).x(3);
        auto y3_old = qvec.At(i).y(3);

        auto x4_old = qvec.At(i).x(4);
        auto y4_old = qvec.At(i).y(4);

        auto x5_old = qvec.At(i).x(5);
        auto y5_old = qvec.At(i).y(5);

        auto x6_old = qvec.At(i).x(6);
        auto y6_old = qvec.At(i).y(6);

        // auto x7_old = qvec.At(i).x(7);
        // auto y7_old = qvec.At(i).y(7);

        // auto x8_old = qvec.At(i).x(8);
        // auto y8_old = qvec.At(i).y(8);
    
        auto c = vec_cor.at(l_idx).At(i);
        // auto [Minv, c] = vec_cor.at(l_idx).At(i);
        
        // if( std::isnan(Minv(0, 0)) ){
        //   new_qvec.At(i).Reset();
        //   continue;
        // }

        auto X1old =  column_t{ 1, x1_old, y1_old, x2_old, y2_old, x3_old, y3_old, x4_old, y4_old, x5_old, y5_old, x6_old, y6_old };
        // auto b = 2 * M.transpose() * X1old;
        // auto b_tilda = Eigen::Matrix<double, NDIM+1, 1>{};
        // b_tilda << b, 1;
        
        auto X1new = 1.0 / ( c.transpose() * X1old ) * X1old;
        
        auto x1_new = static_cast<double>(X1new(1));
        auto y1_new = static_cast<double>(X1new(2));
        auto x2_new = static_cast<double>(X1new(3));
        auto y2_new = static_cast<double>(X1new(4));

        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_new, y2_new );
      }

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Filter( "1 < centrality && centrality < 60" )
    .Filter( "7100 < runId && runId < 8300" )
    .Filter( "-1.0 < vtxX && vtxX < 1.0" )
    .Filter( "-1.0 < vtxY && vtxY < 1.0" )
    .Define("F1_DECOMPOSED", correction_generator(f1_corr, lin), { f1_name, "centrality", "runId" } )
    .Define("F2_DECOMPOSED", correction_generator(f2_corr, lin), { f2_name, "centrality", "runId" } )
    .Define("F3_DECOMPOSED", correction_generator(f3_corr, lin), { f3_name, "centrality", "runId" } )
    .Define("F4_DECOMPOSED", correction_generator(f4_corr, lin), { f4_name, "centrality", "runId" } )

    .Define("Tpos_DECOMPOSED", correction_generator(tn_corr, lin), { tp_name, "centrality", "runId" } )
    .Define("Tneg_DECOMPOSED", correction_generator(tp_corr, lin), { tn_name, "centrality", "runId" } )
    
    .Define("proton_DECOMPOSED", correction_generator(p_corr, lin), { proton_name, "centrality", "runId" } )
  ;

  auto file_out = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( "decomposed_out.root", "RECREATE" ), [](auto f){f ->Close(); } };
  file_out->cd();
  auto tree = new TTree("tree", "tree");
  Qn::DataContainerQVector f1{}, f2{}, f3{}, f4{}, tp{}, tn{}, p{};
  double cent{};
  double r_id{};
  double v_x{};
  double v_y{};
  tree->Branch( "centrality", &cent );
  tree->Branch( "runId", &r_id );
  tree->Branch( "vtxX", &v_x );
  tree->Branch( "vtxY", &v_y );

  tree->Branch( "F1_DECOMPOSED", "Qn::DataContainerQVector", &f1 );
  tree->Branch( "F2_DECOMPOSED", "Qn::DataContainerQVector", &f2 );
  tree->Branch( "F3_DECOMPOSED", "Qn::DataContainerQVector", &f3 );
  tree->Branch( "F4_DECOMPOSED", "Qn::DataContainerQVector", &f4 );
  tree->Branch( "Tpos_DECOMPOSED", "Qn::DataContainerQVector", &tp );
  tree->Branch( "Tneg_DECOMPOSED", "Qn::DataContainerQVector", &tn );
  tree->Branch( "proton_DECOMPOSED", "Qn::DataContainerQVector", &p );

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  dd.Foreach( [tree, &cent, &r_id, &v_x, &v_y, &f1, &f2, &f3, &f4, &tp, &tn, &p]( 
    double centrality,
    double run_id,
    double vtx_x,
    double vtx_y,
    Qn::DataContainerQVector f1_ev, 
    Qn::DataContainerQVector f2_ev, 
    Qn::DataContainerQVector f3_ev, 
    Qn::DataContainerQVector f4_ev, 
    Qn::DataContainerQVector tp_ev, 
    Qn::DataContainerQVector tn_ev, 
    Qn::DataContainerQVector p_ev ) mutable {
    cent = centrality;
    r_id = run_id;
    v_x = vtx_x;
    v_y = vtx_y;
    f1 = f1_ev;
    f2 = f2_ev;
    f3 = f3_ev;
    f4 = f4_ev;
    tp = tp_ev;
    tn = tn_ev;
    p = p_ev;
    tree->Fill();
  }, 
  std::vector<std::string>{ 
    "centrality",
    "runId",
    "vtxX",
    "vtxY",
    "F1_DECOMPOSED",
    "F2_DECOMPOSED",
    "F3_DECOMPOSED",
    "F4_DECOMPOSED",
    "Tpos_DECOMPOSED",
    "Tneg_DECOMPOSED",
    "proton_DECOMPOSED",
  } );
  file_out->cd();
  tree->Write();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Elapsed time: " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " min" << std::endl;
  auto elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  std::cout << "It is " << ( tree->GetEntries() > 0 ? elapsed_s / tree->GetEntries() : 0.0 ) << " μs/ev." << std::endl;
}