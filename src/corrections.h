#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <iostream>
#include <functional>
#include <iterator>
#include <vector>

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
auto MakeWhiteningMatrixFunc() -> mixing_matrix_t<NHARM*2>{
  return [](const std::vector<double>& vec_mean, const std::vector<double>& vec_cov){
    auto M = mixing_matrix_t<NHARM*2>{ mixing_matrix_t<NHARM*2>::Zero() };
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
          cov.push_back( vec_cov.at(i+j) * 2 );
        } i+=4;

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
template<typename correction_matrix_t>
correction_matrix_t PseudoInverse( const correction_matrix_t& M, double l ){
  auto svd = Eigen::JacobiSVD<correction_matrix_t> ( M, Eigen::ComputeThinU | Eigen::ComputeThinV );    
  auto singular_values = svd.singularValues();
  auto U = svd.matrixU();
  auto V = svd.matrixV();
  auto Splus = correction_matrix_t{ correction_matrix_t::Zero() };
  auto rank = size_t{0};
  for (auto i = size_t{0}; i < singular_values.size(); ++i) {
    auto s = singular_values(i);
    if( fabs(s) < l )
      continue;
    // Splus(i, i) = sqrt(0.5 ) / sqrt( s);
    Splus(i, i) = 1.0 / s;
    rank++;
  }
  auto E = correction_matrix_t::Identity();
  auto Vr = V.leftCols(rank);
  auto Etilda = Vr * Vr.transpose() * E;
  auto Mpinv = correction_matrix_t{ V * Splus * U.transpose() };
  std::cout << "l: " << l << "\nMatrix M:\n" << M << "\nMatrix U:\n" << U << "\nS: " << singular_values.transpose() << "\nMatrix V:\n" << V << "\nMatrix S:\n" << Splus << "\nInverse:\n" << Mpinv << "\nEtilda:\n" << Etilda << "\n";
  return Mpinv;
}

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

template<size_t NHARM, typename Func> 
auto MakeCorrectionContainer( std::vector<Qn::DataContainerStatCalculate> vec_mean, std::vector<Qn::DataContainerStatCalculate> vec_cov, const Func& mixing_matrix_generator, const double l=5e-3 ) -> correction_container_t<NHARM*2> {
  auto n_bins = vec_mean.front().size();
  auto axes = vec_mean.front().GetAxes();
  auto correction_container = correction_container_t<NHARM*2>{ axes };
  for(auto bin=size_t{}; bin < n_bins; ++bin ){
    auto vec_means_double = std::vector<double>{};
    vec_means_double.reserve( vec_mean.size() );
    std::for_each( vec_mean.begin(), vec_mean.end(), [bin, &vec_means_double](const auto& c) mutable { vec_means_double.push_back( c[bin].Mean() ); } );
    auto vec_cov_double = std::vector<double>{};
    vec_cov_double.reserve( vec_cov.size() );
    std::for_each( vec_cov.begin(), vec_cov.end(), [bin, &vec_cov_double](const auto& c) mutable { vec_cov_double.push_back( c[bin].Mean() ); } );
    auto M = mixing_matrix_generator( vec_means_double, vec_cov_double );
    auto Minv = PseudoInverse(M, l);
    auto c = column_t<NHARM*2>{};
    std::for_each( vec_means_double.begin(), vec_means_double.end(), [i=size_t{}, &c](auto m) mutable { c(i) = m; ++i; } );
    correction_container[bin] = std::tie( Minv, c );
  }
  return correction_container;
}


template<size_t NHAMR>
class CorrectorBuilder{
public:
  CorrectorBuilder(correction_container_t<NHARM*2> correction_container) :
  correction_container_{ correction_container } {}
  template<typename... Args>
  auto IssueUVectorCorrector() -> UVectorCorrector<NHARM, Args...> {
    return UVectorCorrector<NHAMR, Args...>( correction_container_ );
  }

private:
  correction_container_t<NHARM*2> correction_container_;
};

template<size_t NHARM, typename... Args>
class UVectorCorrector{
public:
  UVectorCorrector(correction_container_t<NHARM*2>&& correction_container) :
  correction_container_{ correction_container } {}

  auto operator()( Args... args ) -> std::vector<std::map<Qn::QVec>> {
    return Execute( args... );
  }

private:
  template<typename... Coord_t>
  auto Execute( uvector_t old_vectors, Coord_t... coordinates ) -> uvector_t{
    auto result = uvector_t{};
    result.reserve( old_vectors.size() );
    for( auto i=0; i<old_vectors.size(); ++i ){
      auto old = old_vectors.at(i);
      auto coord = FormCoordinates( i, coordinates... );
      auto Xold = column_t<NHARM*2>{};
      auto j=size_t{0};
      for( auto p : old ){
        auto x = p.second.x;
        auto y = p.second.y;
        Xold(2*j) = x;
        Xold(2*j+1) = y;
        j++;
      }
      auto [Minv, c] = correction_container_[coord];
      auto Xnew = Minv( Xold - c );
      j=0;
      result.emplace_back();
      for( auto p : old ){
        auto harm = p.first;
        auto x = Xnew(2*j);
        auto y = Xnew(2*j+1);
        result.back()[harm] = Qn::QVec{static_cast<float>(x), static_cast<float>(y)};
        j++;
      }
    }
    return result;
  }
  template<typename T, typename... ColumnTypes>
  std::vector<double> FormCoordinates( size_t i, T first, ColumnTypes... rest ){
    auto vec_coordinates = std::vector<double>{};
    if constexpr ( std::is_floating_point_v<T> ){
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
    if constexpr ( std::is_floating_point_v<T> ){
      return std::vector<double>{ static_cast<double>( coordinate ) };
    }
    else{
      return std::vector<double>{ static_cast<double>( coordinate.at(i) ) };
    }
  }
  correction_container_t<NHARM*2> correction_container_;
};

#endif // CORRECTIONS_H