#ifndef CORRELATION_HELPER_H
#define CORRELATION_HELPER_H

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <type_traits>
#include <DataContainer.hpp>
#include <Axis.hpp>

class CorrelationHelper :  public ROOT::Detail::RDF::RActionImpl<CorrelationHelper>{
public:
   /// This type is a requirement for every helper.
  using Result_t = Qn::DataContainerStatCollect;

  CorrelationHelper( std::vector< Qn::AxisD > vec_axes, size_t n_samples=100 ) : 
  final_result_{ new Qn::DataContainerStatCollect(vec_axes) },
  thread_results_{ std::vector<Qn::DataContainerStatCollect>{} }{
    for( auto i=size_t{0}; i<final_result_->size(); ++i ){
      final_result_->operator[](i).SetNumberOfSamples(n_samples);
    }
    const auto n_slots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
    thread_results_.reserve( n_slots );
    for( auto i=size_t{}; i<n_slots; ++i ){
      thread_results_.emplace_back( *final_result_ );
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
  template <typename T, typename... ColumnTypes>
  void Execute(unsigned int slot, T vec_val, T vec_weights, ROOT::RVec<ULong64_t> vec_samples, ColumnTypes... coordinates){
    if constexpr ( std::is_floating_point_v<T> ){
      auto coord = FormCoordinates( 0, coordinates... );
      auto bin = thread_results_[slot].FindBin( coord );
      thread_results_[slot][ bin ].Fill( vec_val, vec_weights, vec_samples );
    } 
    else {
      for( auto i=size_t{}; i<vec_val.size(); ++i ){
        auto val = vec_val.at(i);
        auto weight = vec_weights.at(i);
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

    if constexpr ( std::is_floating_point_v<typename std::decay_t<Func>::First_t> ){
      df = df.Define( component_names[0], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ return first[h1].x*second[h2].x; } );
      df = df.Define( component_names[1], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ return first[h1].y*second[h2].x; } );
      df = df.Define( component_names[2], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ return first[h1].x*second[h2].y; } );
      df = df.Define( component_names[3], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ return first[h1].y*second[h2].y; } );
    } else {
      df = df.Define( component_names[0], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * second[h2].x ); return res; } } );
      df = df.Define( component_names[1], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * second[h2].x ); return res; } } );
      df = df.Define( component_names[2], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x * second[h2].y ); return res; } } );
      df = df.Define( component_names[3], [h1, h2]( typename std::decay_t<Func>::First_t first, typename std::decay_t<Func>::Second_t second ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y * second[h2].y ); return res; } } );
    }
    vec_res_names.insert( vec_res_names.back(), component_names.begin(), component_names.end() );
  }
  return vec_res_names;
}

template<typename U, typename V>
struct CorrFunc2Part(){
  using First_t = U;
  using Second_t = V;  
};

#endif // CORRELATION_HELPER_H