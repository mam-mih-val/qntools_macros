#ifndef VECTOR_GENERATORS_H
#define VECTOR_GENERATORS_H

#include <string>
#include <vector>
#include <cmath>


const auto u_generator( size_t harmonic, std::function< double(double) > component ){
  return [harmonic, component]( std::vector<float> vec_phi ){
    auto vec_results = std::vector<double>{};
    vec_results.reserve(vec_phi.size());
    for( auto phi : vec_phi ){
      vec_results.push_back( component( harmonic*phi ) );
    }
    return vec_results;
  };
}

const auto cov_generator( size_t h_a, std::function< double(double) > c_a, size_t h_b, std::function< double(double) > c_b ){
  return [h_a, h_b, c_a, c_b]( std::vector<float> vec_phi ){
    auto vec_results = std::vector<double>{};
    vec_results.reserve(vec_phi.size());
    for( auto phi : vec_phi ){
      vec_results.push_back( c_a( h_a*phi ) * c_b( h_b*phi )  );
    }
    return vec_results;
  };
}

template<typename DataFrame>
auto AddUVector( DataFrame& df, const std::string& vector_name, const std::vector<size_t>& harmonics, const std::string& phi_variable_name ) -> std::vector<std::string> {
  auto vec_defined = std::vector<std::string>{};
  vec_defined.reserve( harmonics.size() *2 );
  for( const auto& harm : harmonics ){
    auto x_name = std::string{vector_name}.append(".x").append(std::to_string(harm));
    auto y_name = std::string{vector_name}.append(".x").append(std::to_string(harm));
    df = df.Define( x_name, u_generator(harm, []( double phi ) { return std::cos(phi); } ), std::vector<std::string>{ phi_variable_name } );
    df = df.Define( y_name, u_generator(harm, []( double phi ) { return std::sin(phi); } ), std::vector<std::string>{ phi_variable_name } );
    vec_defined.push_back(x_name);
    vec_defined.push_back(y_name);
  }
  return vec_defined;
}

template<typename DataFrame>
auto AddUVectorCovariance( DataFrame& df, const std::string& vector_name, const std::vector<size_t>& harmonics, const std::string& phi_variable_name ) -> std::vector<std::string> {
  auto vec_defined = std::vector<std::string>{};
  vec_defined.reserve( harmonics.size() * harmonics.size() );
  for( auto i=size_t{0}; i < harmonics.size(); ++i ){
    auto h_a = harmonics[i];
    for( auto j=size_t{i}; j < harmonics.size(); ++j  ){
      auto h_b = harmonics[j];
      
      auto xx_name = std::string{vector_name}.append(".x").append(std::to_string(h_a)).append("x").append(std::to_string(h_b));
      auto yx_name = std::string{vector_name}.append(".y").append(std::to_string(h_a)).append("x").append(std::to_string(h_b));
      auto xy_name = std::string{vector_name}.append(".x").append(std::to_string(h_a)).append("y").append(std::to_string(h_b));
      auto yy_name = std::string{vector_name}.append(".y").append(std::to_string(h_a)).append("y").append(std::to_string(h_b));
      
      df = df.Define( xx_name, cov_generator(h_a, []( double phi ) { return std::cos(phi); }, h_b, []( double phi ) { return std::cos(phi); } ), std::vector<std::string>{ phi_variable_name } );
      df = df.Define( yx_name, cov_generator(h_a, []( double phi ) { return std::sin(phi); }, h_b, []( double phi ) { return std::cos(phi); } ), std::vector<std::string>{ phi_variable_name } );
      df = df.Define( xy_name, cov_generator(h_a, []( double phi ) { return std::cos(phi); }, h_b, []( double phi ) { return std::sin(phi); } ), std::vector<std::string>{ phi_variable_name } );
      df = df.Define( yy_name, cov_generator(h_a, []( double phi ) { return std::sin(phi); }, h_b, []( double phi ) { return std::sin(phi); } ), std::vector<std::string>{ phi_variable_name } );

      vec_defined.push_back( xx_name );
      vec_defined.push_back( yx_name );
      vec_defined.push_back( xy_name );
      vec_defined.push_back( yy_name );
    }
    return vec_defined;
  }
}


#endif // VECTOR_GENERATORS_H