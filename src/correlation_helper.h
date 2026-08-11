#ifndef CORRELATION_HELPER_H
#define CORRELATION_HELPER_H

#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <type_traits>
#include <functional>
#include <DataContainer.hpp>
#include <StatCalculate.hpp>
#include <Axis.hpp>

#include <ROOT/RDataFrame.hxx>
#include "vector_generators.h"


class CorrelationHelper :  public ROOT::Detail::RDF::RActionImpl<CorrelationHelper>{
public:
   /// This type is a requirement for every helper.
  using Result_t = Qn::DataContainerStatCollect;

  CorrelationHelper( std::vector< Qn::AxisD > vec_axes, size_t n_samples=100, Qn::Stat::WeightType weight_type = Qn::Stat::WeightType::OBSERVABLE ) : 
  final_result_{ new Qn::DataContainerStatCollect(vec_axes) },
  thread_results_{ std::vector<Qn::DataContainerStatCollect>{} }{
    for( auto i=size_t{0}; i<final_result_->size(); ++i ){
      final_result_->operator[](i).SetNumberOfSamples(n_samples);
      final_result_->operator[](i).SetWeightType( weight_type );
    }
    const auto n_slots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
    thread_results_.reserve( n_slots );
    for( auto i=size_t{}; i<n_slots; ++i ){
      thread_results_.emplace_back( *final_result_ );
      for( auto& bin : thread_results_.back()  ){
        bin.SetNumberOfSamples(n_samples);
        bin.SetWeightType( weight_type );
      }
    }
  }
  ~CorrelationHelper() = default;

  CorrelationHelper( CorrelationHelper&& ) = default;
  CorrelationHelper( const CorrelationHelper& ) = delete;
  CorrelationHelper& operator=( CorrelationHelper&& ) = default;
  CorrelationHelper& operator=( const CorrelationHelper& ) = delete;

  std::shared_ptr<Qn::DataContainerStatCollect> GetResultPtr() const { return final_result_; }
  void Initialize() {}
  void InitTask(TTreeReader *, unsigned int) {}

  template <typename... ColumnTypes>
  void Exec(unsigned int slot, ColumnTypes... values){ 
    Execute(slot, values...); 
  }

  void Finalize(){
    auto list = TList();
    for( const auto& obj : thread_results_ ){
      auto obj_ptr = new Qn::DataContainerStatCollect( obj );
      list.Add( obj_ptr );
    }
    final_result_->Merge( dynamic_cast<TCollection*>( &list ) );
  }
 
  std::string GetActionName(){
    return std::string{"CorrelationHelper"};
  }

private:
  template <typename T, typename V, typename... ColumnTypes>
  void Execute(unsigned int slot, T vec_val, V vec_weights, ROOT::RVec<ULong64_t> vec_samples, ColumnTypes... coordinates){
    if constexpr ( std::is_floating_point_v<T> ){
      auto coord = FormCoordinates( 0, coordinates... );
      auto bin = thread_results_[slot].FindBin( coord );
      auto weight = static_cast<double>(vec_weights);
      thread_results_[slot][ bin ].Fill( vec_val, weight, vec_samples );
    } 
    else {
      for( auto i=size_t{}; i<vec_val.size(); ++i ){
        auto val = vec_val.at(i);
        auto weight = static_cast<double>(vec_weights.at(i));
        auto coord = FormCoordinates( i, coordinates... );
        auto bin = thread_results_[slot].FindBin( coord );
        if( bin > thread_results_[slot].size() )
          continue;
        if( bin < 0 )
          continue;
        thread_results_[slot][ bin ].Fill( val, weight, vec_samples );
      }
    }
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


  std::shared_ptr<Qn::DataContainerStatCollect> final_result_;
  std::vector< Qn::DataContainerStatCollect > thread_results_;
};

template<typename DataFrame, typename Func>
auto Define2PartCorrelation( DataFrame& df, Func corr_func, const std::string& first_name, const std::string& second_name, const std::vector< std::pair<size_t, size_t> >& harmonics  ) -> std::vector<std::string> {
  auto vec_res_names = std::vector<std::string>{};
  vec_res_names.reserve( 4*harmonics.size() );
  for( const auto& h_pair : harmonics ){
    auto h1 = h_pair.first;
    auto h2 = h_pair.second;
    auto correlation_name = std::string{first_name}.append("_").append(second_name);
    auto component_names = std::vector<std::string>(4);
    component_names[0] = std::string{correlation_name}.append("_x").append(std::to_string(h1)).append("x").append(std::to_string(h2));
    component_names[1] = std::string{correlation_name}.append("_y").append(std::to_string(h1)).append("x").append(std::to_string(h2));
    component_names[2] = std::string{correlation_name}.append("_x").append(std::to_string(h1)).append("y").append(std::to_string(h2));
    component_names[3] = std::string{correlation_name}.append("_y").append(std::to_string(h1)).append("y").append(std::to_string(h2));

    if constexpr ( std::is_same_v<qvector_t, typename Func::First_t> ){
      df = df.Define( component_names[0], [h1, h2]( qvector_t first, qvector_t second ){ return first[h1].x*second[h2].x; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[1], [h1, h2]( qvector_t first, qvector_t second ){ return first[h1].y*second[h2].x; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[2], [h1, h2]( qvector_t first, qvector_t second ){ return first[h1].x*second[h2].y; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[3], [h1, h2]( qvector_t first, qvector_t second ){ return first[h1].y*second[h2].y; }, std::vector{first_name, second_name} );
    } 
    if constexpr ( std::is_same_v<uvector_t, typename Func::First_t> ) {
      df = df.Define( component_names[0], [h1, h2]( uvector_t first, qvector_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * second[h2].x ); } return res; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[1], [h1, h2]( uvector_t first, qvector_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * second[h2].x ); } return res; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[2], [h1, h2]( uvector_t first, qvector_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * second[h2].y ); } return res; }, std::vector{first_name, second_name} );
      df = df.Define( component_names[3], [h1, h2]( uvector_t first, qvector_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * second[h2].y ); } return res; }, std::vector{first_name, second_name} );
    }
    vec_res_names.insert( vec_res_names.end(), component_names.begin(), component_names.end() );
  }
  return vec_res_names;
}

template<typename U, typename V>
struct CorrFunc2Part{
  using First_t = U;
  using Second_t = U;
};

template<typename DataFrame, typename Func>
auto DefineVectorCovariance( DataFrame& df, Func corr_func, const std::string& first_name, const std::vector<size_t>& harmonics  ) -> std::vector<std::string> {
  auto vec_res_names = std::vector<std::string>{};
  vec_res_names.reserve( 4*harmonics.size() );
  for( auto i=size_t{0}; i<harmonics.size(); ++i ){
    auto h1 = harmonics.at(i);
    for( auto j=i; j<harmonics.size(); ++j ){
      auto h2 = harmonics.at(j);
      auto component_names = std::vector<std::string>(4);
      component_names[0] = std::string{first_name}.append("_x").append(std::to_string(h1)).append("x").append(std::to_string(h2));
      component_names[1] = std::string{first_name}.append("_y").append(std::to_string(h1)).append("x").append(std::to_string(h2));
      component_names[2] = std::string{first_name}.append("_x").append(std::to_string(h1)).append("y").append(std::to_string(h2));
      component_names[3] = std::string{first_name}.append("_y").append(std::to_string(h1)).append("y").append(std::to_string(h2));

      if constexpr ( std::is_same_v<qvector_t, typename Func::First_t> ){
        df = df.Define( component_names[0], [h1, h2]( qvector_t first ){ return first[h1].x*first[h2].x; }, std::vector{first_name} );
        df = df.Define( component_names[1], [h1, h2]( qvector_t first ){ return first[h1].y*first[h2].x; }, std::vector{first_name} );
        df = df.Define( component_names[2], [h1, h2]( qvector_t first ){ return first[h1].x*first[h2].y; }, std::vector{first_name} );
        df = df.Define( component_names[3], [h1, h2]( qvector_t first ){ return first[h1].y*first[h2].y; }, std::vector{first_name} );
      } 
      if constexpr ( std::is_same_v<uvector_t, typename Func::First_t> ) {
        df = df.Define( component_names[0], [h1, h2]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * f[h2].x ); } return res; }, std::vector{first_name} );
        df = df.Define( component_names[1], [h1, h2]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * f[h2].x ); } return res; }, std::vector{first_name} );
        df = df.Define( component_names[2], [h1, h2]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * f[h2].y ); } return res; }, std::vector{first_name} );
        df = df.Define( component_names[3], [h1, h2]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * f[h2].y ); } return res; }, std::vector{first_name} );
      }
  
      vec_res_names.insert( vec_res_names.end(), component_names.begin(), component_names.end() );
    }
  }
  
  return vec_res_names;
}

template<typename DataFrame, typename Func>
auto DefineVectorMeans( DataFrame& df, Func corr_func, const std::string& first_name, const std::vector<size_t>& harmonics  ) -> std::vector<std::string> {
  auto vec_res_names = std::vector<std::string>{};
  vec_res_names.reserve( 2*harmonics.size() );
  for( auto i=size_t{0}; i<harmonics.size(); ++i ){
    auto h1 = harmonics.at(i);
    
    auto component_names = std::vector<std::string>(4);
    component_names[0] = std::string{first_name}.append("_x").append(std::to_string(h1));
    component_names[1] = std::string{first_name}.append("_y").append(std::to_string(h1));

    if constexpr ( std::is_same_v<qvector_t, typename Func::First_t> ){
      df = df.Define( component_names[0], [h1]( qvector_t first ){ return first[h1].x; }, std::vector{first_name} );
      df = df.Define( component_names[1], [h1]( qvector_t first ){ return first[h1].y; }, std::vector{first_name} );
    } 
    if constexpr ( std::is_same_v<uvector_t, typename Func::First_t> ) {
      df = df.Define( component_names[0], [h1]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x ); } return res; }, std::vector{first_name} );
      df = df.Define( component_names[1], [h1]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y ); } return res; }, std::vector{first_name} );
    }

    vec_res_names.insert( vec_res_names.end(), component_names.begin(), component_names.end() );
    
  }
  
  return vec_res_names;
}

template<typename U, typename V>
struct CorrFunc1Part{
  using First_t = U;
};

#endif // CORRELATION_HELPER_H