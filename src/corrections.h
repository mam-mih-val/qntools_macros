#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <functional>
#include <iterator>
#include <limits>
#include <tuple>
#include <vector>
#include <type_traits>

#include <Eigen/Dense>

#include <DataContainer.hpp>
#include <Axis.hpp>
#include <QVector.hpp>

#include <vector_generators.h>

template<size_t NDIM>
using mixing_matrix_t = Eigen::Matrix<double, NDIM, NDIM>;

template<size_t NDIM>
using column_t = Eigen::Matrix<double, NDIM, 1>;

template<size_t NDIM>
using correction_container_t = Qn::DataContainer< std::tuple< mixing_matrix_t<NDIM>, column_t<NDIM> >, Qn::AxisD >;

template<size_t NHARM>
class InverseCovariance{
public:
  using correction_matrix_t = Eigen::Matrix<double, NHARM*2, NHARM*2>;

  InverseCovariance() = default;
  ~InverseCovariance() = default;
  
  static auto MixingMatrix() -> std::function< correction_matrix_t( std::vector<double>, std::vector<double> ) > {
    return [](const std::vector<double>& vec_mean, const std::vector<double>& vec_cov) -> correction_matrix_t {
      auto M = correction_matrix_t{ correction_matrix_t::Zero() };
      auto i = size_t {0};
      for( auto h_a = size_t{0}; h_a < NHARM; ++h_a ){
        auto x_a = vec_mean[2*h_a];
        auto y_a = vec_mean[2*h_a+1];
        for( auto h_b = h_a; h_b < NHARM; ++h_b ){
          auto x_b = vec_mean[2*h_b];
          auto y_b = vec_mean[2*h_b+1];
          auto cov = std::vector<double>{}; 
          cov.reserve(4);
          for( auto j=size_t{0}; j<4; ++j ){
            cov.push_back( vec_cov.at(i+j) );
          } i+=4;

          if( h_a != h_b )
            continue;

          cov[0] -= x_a*x_b;
          cov[3] -= y_a*y_b;

          M( 2*h_a, 2*h_b ) = cov[0];
          M( 2*h_a+1, 2*h_b+1 ) = cov[3];
        }
      }
      return M;
    };
  }

  static auto PseudoInverse() -> std::function< correction_matrix_t( correction_matrix_t, double ) > {
    return [](const correction_matrix_t& M, double l){
      auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
      auto singular_values = svd.singularValues();
      auto sv_sum = double{0.0};
      auto U = svd.matrixU();
      auto V = svd.matrixV();
      auto Splus = correction_matrix_t{ correction_matrix_t::Zero() };
      auto rank = size_t{0};
      for (auto i = size_t{0}; i < singular_values.size(); ++i) {
        auto s = singular_values(i);
        if( s / 0.5  < l )
          continue;
        Splus(i, i) = 0.5 / s;
        rank++;
        sv_sum += s;
      }
      auto Ur = U.leftCols(rank);
      auto Ur1 = correction_matrix_t{ correction_matrix_t::Zero() };
      auto UrUrT = Ur*Ur.transpose();
      for( auto r = size_t{0}; r < Ur.rows(); r++ ){
        for( auto c = size_t{0}; c < Ur.cols(); c++ ){
          if( fabs( UrUrT(r, r) ) < l ) 
            continue;
          Ur1(r, c) = Ur(r, c) / UrUrT(r, r);
        }
      }
      auto Mpinv = correction_matrix_t{ Ur1 * Splus * U.transpose() };
      std::cout << "l: " << l << "\nMatrix M:\n" << M << "\nMatrix U:\n" << Ur << "\nS: " << singular_values.transpose() << "\nMatrix S:\n" << Splus << "\nInverse:\n" << Mpinv << "\nE:\n" << Ur1 << "\nUrUr^T\n" << UrUrT << "\n\n";
      return Mpinv;
    };
  }
};


template<size_t NHARM>
class PrincipalComponents{
public:
  using correction_matrix_t = Eigen::Matrix<double, NHARM*2, NHARM*2>;

  PrincipalComponents() = default;
  ~PrincipalComponents() = default;
  
  static auto MixingMatrix() -> std::function< correction_matrix_t( std::vector<double>, std::vector<double> ) > {
    return [](const std::vector<double>& vec_mean, const std::vector<double>& vec_cov) -> correction_matrix_t {
      auto M = correction_matrix_t{ correction_matrix_t::Zero() };
      auto i = size_t {0};
      for( auto h_a = size_t{0}; h_a < NHARM; ++h_a ){
        auto x_a = vec_mean[2*h_a];
        auto y_a = vec_mean[2*h_a+1];
        for( auto h_b = h_a; h_b < NHARM; ++h_b ){
          auto x_b = vec_mean[2*h_b];
          auto y_b = vec_mean[2*h_b+1];
          auto cov = std::vector<double>{}; 
          cov.reserve(4);
          for( auto j=size_t{0}; j<4; ++j ){
            cov.push_back( vec_cov.at(i+j) );
          } i+=4;

          cov[0] -= x_a*x_b;
          cov[1] -= y_a*x_b;
          cov[2] -= x_a*y_b;
          cov[3] -= y_a*y_b;

          M( 2*h_a, 2*h_b ) = cov[0];
          M( 2*h_a+1, 2*h_b ) = cov[1];
          M( 2*h_a, 2*h_b+1 ) = cov[2];
          M( 2*h_a+1, 2*h_b+1 ) = cov[3];

          M( 2*h_b, 2*h_a ) = cov[0];
          M( 2*h_b, 2*h_a+1 ) = cov[1];
          M( 2*h_b+1, 2*h_a ) = cov[2];
          M( 2*h_b+1, 2*h_a+1 ) = cov[3];
        }
      }
      return M;
    };
  }

  static auto PseudoInverse() -> std::function< correction_matrix_t( correction_matrix_t, column_t<NHARM*2>, double ) > {
    return [](const correction_matrix_t& M, column_t<NHARM*2> sv_variance, double l){
      auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
      auto trace = M.trace();
      auto singular_values = svd.singularValues();
      auto sv_sum = double{0.0};
      auto U = svd.matrixU();
      auto V = svd.matrixV();
      auto Splus = correction_matrix_t{ correction_matrix_t::Zero() };
      auto rank = size_t{0};
      for (auto i = size_t{0}; i < singular_values.size(); ++i) {
        if( sv_sum > l * trace )
          continue;
        auto s = singular_values(i);
        Splus(i, i) = 0.5 / s;
        rank++;
        sv_sum+=s;
      }
      std::cout << "\nrank: " << rank << "\n";
      std::cout << "trace: " << trace << " used trace: " << sv_sum << "\n\n";
      auto Ur = U.leftCols(rank);
      auto Ur1 = correction_matrix_t{ correction_matrix_t::Zero() };
      auto UrUrT = Ur*Ur.transpose();
      for( auto r = size_t{0}; r < Ur.rows(); r++ ){
        for( auto c = size_t{0}; c < Ur.cols(); c++ ){
          if( fabs( UrUrT(r, r) ) < 5e-2 )
            continue;
          Ur1(r, c) = Ur(r, c) / UrUrT(r, r);
          // if( fabs( UrUrT(r, r) ) < 5e-2 )
          //   Ur1(r, c) = std::nan("");
        }
      }
      auto Mpinv = correction_matrix_t{ Ur1 * Splus * U.transpose() };
      std::cout << "l: " << l << "\nMatrix M:\n" << M << "\n";
      std::cout << "S: " << singular_values.transpose() << "\nV: " << sv_variance.transpose() << "\n";
      std::cout << "\nMatrix S:\n" << Splus << "\nInverse:\n" << Mpinv << "\nE:\n" << Ur1 << "\nUrUr^T\n" << UrUrT << "\n\n";
      return Mpinv;
    };
  }
};

template<size_t NHARM>
class TwistRescale{
public:
  using correction_matrix_t = Eigen::Matrix<double, NHARM*2, NHARM*2>;

  TwistRescale() = default;
  ~TwistRescale() = default;
  
  static auto MixingMatrix() -> std::function< correction_matrix_t( std::vector<double>, std::vector<double> ) > {
    return [](const std::vector<double>& vec_mean, const std::vector<double>& vec_cov) -> correction_matrix_t {
      auto M = correction_matrix_t{ correction_matrix_t::Zero() };
      auto i = size_t {0};
      for( auto h_a = size_t{0}; h_a < NHARM; ++h_a ){
        auto x_a = vec_mean[2*h_a];
        auto y_a = vec_mean[2*h_a+1];

        auto x_2a = vec_mean[4*h_a+2];
        auto y_2a = vec_mean[4*h_a+3];

        M( 2*h_a, 2*h_a ) = 1.0 + x_2a;
        M( 2*h_a+1, 2*h_a ) = y_2a;
        M( 2*h_a, 2*h_a+1 ) = y_2a;
        M( 2*h_a+1, 2*h_a+1 ) = 1.0 - x_2a;
      }
      return M;
    };
  }

  static auto PseudoInverse() -> std::function< correction_matrix_t( correction_matrix_t, column_t<NHARM*2>, double ) > {
    return [](const correction_matrix_t& M, column_t<NHARM*2> sv_variance, double l){
      auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
      auto singular_values = svd.singularValues();
      auto U = svd.matrixU();
      auto V = svd.matrixV();
      auto Splus = correction_matrix_t{ correction_matrix_t::Zero() };
      auto rank = size_t{0};
      for (auto i = size_t{0}; i < singular_values.size(); ++i) {
        auto s = singular_values(i);
        auto s_var = sv_variance(i);
        auto s_err = sqrt( s_var - s*s );
        if( s < l )
          continue;
        Splus(i, i) = 1.0 / s;
        rank++;
      }
      auto Ur = U.leftCols(rank);
      auto Mpinv = correction_matrix_t{ U * Splus * U.transpose() };
      std::cout << "l: " << l << "\nMatrix M:\n" << M << "\nMatrix U:\n" << Ur << "\nS: " << singular_values.transpose() << "\nMatrix S:\n" << Splus << "\nInverse:\n" << Mpinv << "\n\n";
      return Mpinv;
    };
  }
};

template<size_t NHARM>
class Recenter{
public:
  using correction_matrix_t = Eigen::Matrix<double, NHARM*2, NHARM*2>;

  Recenter() = default;
  ~Recenter() = default;
  
  static auto MixingMatrix() -> std::function< correction_matrix_t( std::vector<double>, std::vector<double> ) > {
    return [](const std::vector<double>& vec_mean, const std::vector<double>& vec_cov) -> correction_matrix_t {
      auto M = correction_matrix_t{ correction_matrix_t::Identity() };
      return M;
    };
  }

  static auto PseudoInverse() -> std::function< correction_matrix_t( correction_matrix_t, double ) > {
    return [](const correction_matrix_t& M, double l){
      return M;
    };
  }
};

template<size_t NHARM>
inline auto ReadMeanCov( std::string str_vec_name, TFile* calib_file ) -> std::tuple< std::vector<Qn::DataContainerStatCalculate>, std::vector<Qn::DataContainerStatCalculate> >{
  std::cout << __func__ << std::endl;
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_mean = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_cov = std::vector<Qn::DataContainerStatCalculate>{};
  vec_mean.reserve(16);
  vec_cov.reserve(13);

  for( auto h_a = size_t{1}; h_a <= NHARM; ++h_a ){
    auto corr_name = str_vec_name+"_x"+std::to_string(h_a);
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_mean.emplace_back( *tmp );
    
    corr_name = str_vec_name+"_y"+std::to_string(h_a);
    std::cout << "Extracting " << corr_name << "\n";
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_mean.emplace_back( *tmp );

    for( auto h_b = h_a; h_b <= NHARM; ++h_b ){
      corr_name = str_vec_name+"_x"+std::to_string(h_a)+"x"s+std::to_string(h_b);
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );

      corr_name = str_vec_name+"_y"+std::to_string(h_a)+"x"s+std::to_string(h_b);
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );

      corr_name = str_vec_name+"_x"+std::to_string(h_a)+"y"s+std::to_string(h_b);
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp );

      corr_name = str_vec_name+"_y"+std::to_string(h_a)+"y"s+std::to_string(h_b);
      std::cout << "Extracting " << corr_name << "\n";
      calib_file->GetObject( corr_name.c_str(), tmp );
      assert(tmp);
      vec_cov.emplace_back( *tmp ); 
    }
  }
  return {vec_mean, vec_cov};
}

template<size_t NHARM, typename CorrectionFunc> 
auto MakeMixingMatrix( const std::vector<Qn::StatCalculate>& vec_stat_mean, const std::vector<Qn::StatCalculate>& vec_stat_cov, const CorrectionFunc& correction_func ) -> mixing_matrix_t<NHARM*2> {
  if( vec_stat_mean.front().SumWeights() < std::numeric_limits<double>::epsilon() )
    return mixing_matrix_t<NHARM*2>{ mixing_matrix_t<NHARM*2>::Identity() };
  auto vec_means_double = std::vector<double>{};
  vec_means_double.reserve( vec_stat_mean.size() );
  std::for_each( vec_stat_mean.begin(), vec_stat_mean.end(), [&vec_means_double](const auto& stat) mutable { vec_means_double.push_back( stat.Mean() ); } );
  auto vec_cov_double = std::vector<double>{};
  vec_cov_double.reserve( vec_stat_cov.size() );
  std::for_each( vec_stat_cov.begin(), vec_stat_cov.end(), [&vec_cov_double](const auto& stat) mutable { vec_cov_double.push_back( stat.Mean() ); } );
  auto mixing_matrix_function = CorrectionFunc::MixingMatrix();
  auto M = mixing_matrix_function( vec_means_double, vec_cov_double );
  return M;
}

auto ExtractBootstrapSamplesMeanWeights(const std::vector<Qn::StatCalculate>& vec_stat ) -> std::tuple<  std::vector<std::vector<double>>,  std::vector<std::vector<double>> > {
  auto v2_means = std::vector<std::vector<double>>{};
  v2_means.reserve( vec_stat.size() );
  auto v2_weights = std::vector<std::vector<double>>{};
  v2_weights.reserve( vec_stat.size() );
  std::for_each( vec_stat.begin(), vec_stat.end(), [&v2_means, &v2_weights](const auto& stat) mutable { 
    v2_means.push_back( stat.GetSampleMeans() );
    v2_weights.push_back( stat.GetSampleWeights() );
  } );
  auto n_samples = v2_means.front().size();
  auto v2_means_T = std::vector<std::vector<double>>( n_samples );
  auto v2_weights_T = std::vector<std::vector<double>>( n_samples );
  for( auto i=size_t{}; i<n_samples; ++i ){
    v2_means_T.reserve( vec_stat.size() );
    v2_weights_T.reserve( vec_stat.size() );
    for( auto j=size_t{}; j<vec_stat.size(); ++j ){
      v2_means_T[i].push_back( v2_means[j][i] );
      v2_weights_T[i].push_back( v2_weights[j][i] );
    }
  }
  return std::tie( v2_means_T, v2_weights_T );
}

template<size_t NHARM, typename CorrectionFunc> 
auto EstimateVarianceEV( const std::vector<Qn::StatCalculate>& vec_stat_mean, const std::vector<Qn::StatCalculate>& vec_stat_cov, const CorrectionFunc& correction_func ) -> column_t<NHARM*2> {
  struct statistic{
    double sumw{0};
    double sumwx{0};
    double sumwx2{0};
    size_t entries{0};
  };
  using correction_matrix_t = Eigen::Matrix<double, NHARM*2, NHARM*2>;
  auto vec_statistics = std::vector<statistic>( 2*NHARM );
  auto mixing_matrix_function = CorrectionFunc::MixingMatrix();
  auto [v2_mean_v, v2_mean_w] = ExtractBootstrapSamplesMeanWeights(vec_stat_mean);
  auto [v2_cov_v, v2_cov_w] = ExtractBootstrapSamplesMeanWeights(vec_stat_cov);
  auto nsamples = v2_mean_v.size();

  for( auto i=size_t{0}; i<nsamples; ++i ){
    auto M = mixing_matrix_function( v2_mean_v[i], v2_cov_v[i] );
    auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
    auto singular_values = svd.singularValues();
    for( auto j=size_t{}; j<singular_values.size(); ++j ){
      vec_statistics[j].sumwx += singular_values(j);
      vec_statistics[j].sumwx2 += singular_values(j)*singular_values(j);
      vec_statistics[j].sumw += 1;
      vec_statistics[j].entries += 1;
    }
  }

  auto variance = column_t<NHARM*2>{};
  std::for_each( vec_statistics.begin(), vec_statistics.end(), [i=0, &variance]( auto& stat ) mutable { 
    auto sumwx = stat.sumwx;
    auto sumwx2 = stat.sumwx2;
    auto sumw = stat.sumw;
    auto entries = stat.entries;
    if( entries <= 1 ){
      variance(i) = 0;
      return;
    }
    auto mean = sumwx / sumw;
    auto mean2 = sumwx2 / sumw;
    variance(i) = std::sqrt( ( mean2 - mean*mean ) / (entries - 1.) ); 
    i++; 
  } );
  return variance;
}

template<size_t NHARM, typename CorrectionFunc> 
auto MakeCorrectionContainer( std::vector<Qn::DataContainerStatCalculate> vec_mean, std::vector<Qn::DataContainerStatCalculate> vec_cov, const CorrectionFunc& correction_func, const double l=5e-3 ) -> correction_container_t<NHARM*2> {
  auto n_bins = vec_mean.front().size();
  auto axes = vec_mean.front().GetAxes();
  auto correction_container = correction_container_t<NHARM*2>{ axes };
  auto pseudo_inverse_function = CorrectionFunc::PseudoInverse();
  for(auto bin=size_t{}; bin < n_bins; ++bin ){
    auto vec_means_stat = std::vector<Qn::StatCalculate>{};
    vec_means_stat.reserve( vec_mean.size() );
    std::for_each( vec_mean.begin(), vec_mean.end(), [bin, &vec_means_stat](const auto& c) mutable { vec_means_stat.push_back( c[bin] ); } );
    auto vec_cov_stat = std::vector<Qn::StatCalculate>{};
    vec_cov_stat.reserve( vec_cov.size() );
    std::for_each( vec_cov.begin(), vec_cov.end(), [bin, &vec_cov_stat](const auto& c) mutable { vec_cov_stat.push_back( c[bin] ); } );
    auto M = MakeMixingMatrix<NHARM, CorrectionFunc>( vec_means_stat, vec_cov_stat, correction_func );
    auto lambda_var = EstimateVarianceEV<NHARM, CorrectionFunc>( vec_means_stat, vec_cov_stat, correction_func );
    auto Minv = pseudo_inverse_function(M, lambda_var, l);
    auto c = column_t<NHARM*2>{};
    std::for_each( vec_means_stat.begin(), vec_means_stat.end(), [i=size_t{}, &c](auto m) mutable { c(i) = m.Mean(); ++i; } );
    correction_container[bin] = std::tie( Minv, c );
  }
  return correction_container;
}

template<size_t NHARM, typename... Args>
class QVectorCorrector{
public:
  QVectorCorrector(correction_container_t<NHARM*2> correction_container) :
  correction_container_{ correction_container } {}

  auto operator()( Args... args ) -> qvector_t {
    return Execute( args... );
  }

private:
  template<typename... Coord_t>
  auto Execute( qvector_t old_vectors, Coord_t... coordinates ) -> qvector_t{
    auto result = qvector_t{};
    auto coord = FormCoordinates( 0, coordinates... );
    auto bin = correction_container_.FindBin( coord );
    if( bin > correction_container_.size() )
        return result;
      if( bin < 0 )
        return result;
    auto Xold = column_t<NHARM*2>{};
    auto j=size_t{0};
    for( auto p : old_vectors ){
      auto x = p.second.x;
      auto y = p.second.y;
      Xold(2*j) = x;
      Xold(2*j+1) = y;
      j++;
    }
    auto [Minv, c] = correction_container_[bin];
    auto Xnew = Minv*( Xold - c );
    j=0;
    for( auto p : old_vectors ){
      auto harm = p.first;
      auto x = Xnew(2*j);
      auto y = Xnew(2*j+1);
      result[harm] = Qn::QVec{ static_cast<float>(x), static_cast<float>(y) };
      j++;
    }
    return result;
  }
  template<typename T, typename... ColumnTypes>
  std::vector<double> FormCoordinates( size_t i, T first, ColumnTypes... rest ){
    auto vec_coordinates = std::vector<double>{};
    if constexpr ( std::is_floating_point_v<T> || std::is_integral_v<T> ){
      vec_coordinates.push_back(static_cast<double>( first ) );
    } else {
      vec_coordinates.push_back(static_cast<double>( first.at(i) ) );
    }
    auto vec_rest_coord = FormCoordinates( i, rest... );
    vec_coordinates.insert( vec_coordinates.end(), vec_rest_coord.begin(), vec_rest_coord.end() );
    return vec_coordinates;
  }
  template<typename T, typename... ColumnTypes>
  std::vector<double> FormCoordinates( size_t i, T coordinate ){
    if constexpr ( std::is_floating_point_v<T> || std::is_integral_v<T> ){
      return std::vector<double>{ static_cast<double>( coordinate ) };
    }
    else{
      return std::vector<double>{ static_cast<double>( coordinate.at(i) ) };
    }
  }
  correction_container_t<NHARM*2> correction_container_;
};

template<size_t NHARM, typename... Args>
class UVectorCorrector{
public:
  UVectorCorrector(correction_container_t<NHARM*2> correction_container) :
  correction_container_{ correction_container } {}

  auto operator()( Args... args ) -> uvector_t {
    return Execute( args... );
  }

private:
  template<typename... Coord_t>
  auto Execute( uvector_t old_vectors, Coord_t... coordinates ) -> uvector_t{
    auto result = uvector_t{};
    result.reserve( old_vectors.size() );
    for( auto i=0; i<old_vectors.size(); ++i ){
      result.emplace_back();
      auto old = old_vectors.at(i);
      auto coord = FormCoordinates( i, coordinates... );
      auto bin = correction_container_.FindBin( coord );
      if( bin > correction_container_.size() )
        continue;
      if( bin < 0 )
        continue;
      auto Xold = column_t<NHARM*2>{};
      auto j=size_t{0};
      for( auto p : old ){
        auto x = p.second.x;
        auto y = p.second.y;
        Xold(2*j) = x;
        Xold(2*j+1) = y;
        j++;
      }
      auto [Minv, c] = correction_container_[bin];
      auto Xnew = Minv*( Xold - c );
      j=0;
      for( auto p : old ){
        auto harm = p.first;
        auto x = Xnew(2*j);
        auto y = Xnew(2*j+1);
        result.back()[harm] = Qn::QVec{static_cast<float>(x), static_cast<float>(y)};
        j++;
      }
    }
    std::cout << "here" << std::endl;
    return result;
  }
  template<typename T, typename... ColumnTypes>
  std::vector<double> FormCoordinates( size_t i, T first, ColumnTypes... rest ){
    auto vec_coordinates = std::vector<double>{};
    if constexpr ( std::is_floating_point_v<T> || std::is_integral_v<T> ){
      vec_coordinates.push_back(static_cast<double>( first ) );
    } else {
      vec_coordinates.push_back(static_cast<double>( first.at(i) ) );
    }
    auto vec_rest_coord = FormCoordinates( i, rest... );
    vec_coordinates.insert( vec_coordinates.end(), vec_rest_coord.begin(), vec_rest_coord.end() );
    return vec_coordinates;
  }
  template<typename T, typename... ColumnTypes>
  std::vector<double> FormCoordinates( size_t i, T coordinate ){
    if constexpr ( std::is_floating_point_v<T> || std::is_integral_v<T> ){
      return std::vector<double>{ static_cast<double>( coordinate ) };
    }
    else{
      return std::vector<double>{ static_cast<double>( coordinate.at(i) ) };
    }
  }
  correction_container_t<NHARM*2> correction_container_;
};

template<size_t NHARM>
class CorrectorBuilder{
public:
  CorrectorBuilder(correction_container_t<NHARM*2> correction_container) :
  correction_container_{ correction_container } {}
  
  template<typename... Args>
  auto IssueUVectorCorrector() -> UVectorCorrector<NHARM, Args...> {
    return UVectorCorrector<NHARM, Args...>( correction_container_ );
  }

  template<typename... Args>
  auto IssueQVectorCorrector() -> QVectorCorrector<NHARM, Args...> {
    return QVectorCorrector<NHARM, Args...>( correction_container_ );
  }

private:
  correction_container_t<NHARM*2> correction_container_;
};

#endif // CORRECTIONS_H