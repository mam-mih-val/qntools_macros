#ifndef CORRELATION_HELPER_H
#define CORRELATION_HELPER_H

#include <cstddef>
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
      if( bin < 0 )
          return;
      if( bin > thread_results_[slot].size() )
          return;
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
      df = df.Define( component_names[0], [h1, h2]( qvector_t first, qvector_t second ){ return static_cast<double>(first[h1].x*second[h2].x); }, std::vector{first_name, second_name} );
      df = df.Define( component_names[1], [h1, h2]( qvector_t first, qvector_t second ){ return static_cast<double>(first[h1].y*second[h2].x); }, std::vector{first_name, second_name} );
      df = df.Define( component_names[2], [h1, h2]( qvector_t first, qvector_t second ){ return static_cast<double>(first[h1].x*second[h2].y); }, std::vector{first_name, second_name} );
      df = df.Define( component_names[3], [h1, h2]( qvector_t first, qvector_t second ){ return static_cast<double>(first[h1].y*second[h2].y); }, std::vector{first_name, second_name} );
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
        df = df.Define( component_names[0], [h1, h2]( qvector_t first ){ return static_cast<double>(first[h1].x*first[h2].x); }, std::vector{first_name} );
        df = df.Define( component_names[1], [h1, h2]( qvector_t first ){ return static_cast<double>(first[h1].y*first[h2].x); }, std::vector{first_name} );
        df = df.Define( component_names[2], [h1, h2]( qvector_t first ){ return static_cast<double>(first[h1].x*first[h2].y); }, std::vector{first_name} );
        df = df.Define( component_names[3], [h1, h2]( qvector_t first ){ return static_cast<double>(first[h1].y*first[h2].y); }, std::vector{first_name} );
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
    
    auto component_names = std::vector<std::string>(2);
    component_names[0] = std::string{first_name}.append("_x").append(std::to_string(h1));
    component_names[1] = std::string{first_name}.append("_y").append(std::to_string(h1));

    if constexpr ( std::is_same_v<qvector_t, typename Func::First_t> ){
      df = df.Define( component_names[0], [h1]( qvector_t first ){ return static_cast<double>(first[h1].x); }, std::vector{first_name} );
      df = df.Define( component_names[1], [h1]( qvector_t first ){ return static_cast<double>(first[h1].y); }, std::vector{first_name} );
    } 
    if constexpr ( std::is_same_v<uvector_t, typename Func::First_t> ) {
      df = df.Define( component_names[0], [h1]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].x ); } return res; }, std::vector{first_name} );
      df = df.Define( component_names[1], [h1]( uvector_t first ){ std::vector<double> res{}; res.reserve( first.size() ); for( auto f : first ){ res.push_back( f[h1].y ); } return res; }, std::vector{first_name} );
    }

    vec_res_names.insert( vec_res_names.end(), component_names.begin(), component_names.end() );
    
  }
  
  return vec_res_names;
}

template<typename U>
struct CorrFunc1Part{
  using First_t = U;
};

struct x{
  float operator()(Qn::QVec vec){
    return vec.x;
  }
};

struct y{
  float operator()(Qn::QVec vec){
    return vec.x;
  }
};

template<typename T>
auto MakeComponent( T component ){
  return [component](Qn::QVec vec){ return component(vec); };
}

template<typename RetType, typename... Args>
class Correlator{
public:
  Correlator( std::vector<size_t> harmonics, std::vector< std::function< float(Qn::QVec) > > components ) : 
    harmonics_( std::move(harmonics) ), components_( std::move(components_) ) {}
  auto operator()( Args... args ) -> RetType {
    counter_=0;
    Exec(args...);
  }
private:
  std::vector<size_t> harmonics_{0};
  std::vector< std::function< float(Qn::QVec) > > components_{0};
  size_t counter_{0};
  template< typename First, typename... Last >
  auto Exec( First first, Last... last ) -> RetType {
    auto result = RetType{};
    if constexpr( std::is_floating_point_v<RetType> ){
      result = components[counter_](first[ harmonics_[counter_] ]);
      counter_++;
      result *= Exec( last... );
    } else {
      for( auto i = 0; i<first.size(); ++i ){
        result.push_back( components[counter_]( first[i][ harmonics_[counter_] ] ) );
      }
      counter_++;
      auto rest_result = Exec( last... );
      for( auto i = 0; i<first.size(); ++i ){
        result[i] *= rest_result;
      }
    }

    return result;
  }
  template< typename First>
  auto Exec( First first ) -> RetType {
    auto result = RetType{};
    if constexpr( std::is_floating_point_v<RetType> ){
      result = components[counter_](first[ harmonics_[counter_] ]);
    } else {
      for( auto i = 0; i<first.size(); ++i ){
        result.push_back( components[counter_]( first[i][ harmonics_[counter_] ] ) );
      }
    }
    return result;
  }
};

template<typename First, typename... Args>
class CorrelationDecorator{
public:
  using First_t = First;
  CorrelationDecorator( std::vector<std::string> vector_names, std::vector<size_t> harmonics ) : 
    vector_names_(std::move(vector_names)), 
    harmonics_(std::move(harmonics)) {}
  template<typename DF>
  auto operator()( DF& df ) const -> std::vector<std::string> {
    auto general_correlation_name = std::string{};
    std::for_each( vector_names_.begin(), vector_names_.end(), [&general_correlation_name]( const auto& name ) mutable { general_correlation_name.append(name).append("_"); } );
    general_correlation_name.pop_back();
    auto vec_components = std::vector< std::vector< std::function<float(Qn::QVec)> > >{ std::vector< std::function<float(Qn::QVec)> >{} };
    auto vec_corr_names = std::vector< std::string >{ general_correlation_name };
    for( size_t i=0; i<vector_names_.size(); ++i ){
      auto upd_vec_corr_names = std::vector<std::string>{};
      auto upd_vec_components = std::vector< std::vector< std::function<float(Qn::QVec)> > >{};
      for( auto j=0; j < vec_corr_names.size(); ++j ){
        upd_vec_corr_names.push_back( vec_corr_names[j]+"_x"+std::to_string(harmonics_[i]) );
        upd_vec_corr_names.push_back( vec_corr_names[j]+"_y"+std::to_string(harmonics_[i]) );

        auto curr_component_layout = vec_components.at(j);
        upd_vec_components.push_back( curr_component_layout.push_back( MakeComponent(x{}) ) );
        curr_component_layout = vec_components.at(j);
        upd_vec_components.push_back( curr_component_layout.push_back( MakeComponent(y{}) ) );
      }
      vec_corr_names = std::move(upd_vec_corr_names);
      vec_components = std::move(upd_vec_components);
    }

    if constexpr( std::is_same_v<qvector_t, First> ){
      for( auto i=size_t{0}; i<vec_corr_names.size(); ++i ){
        df = df.Define( vec_corr_names[i], Correlator<double, First, Args...>{ harmonics_, vec_components[i] }, vector_names_ );
      }
    } else {
      for( auto i=size_t{0}; i<vec_corr_names.size(); ++i ){
        df = df.Define( vec_corr_names[i], Correlator< std::vector<double>, First, Args...>{ harmonics_, vec_components[i] }, vector_names_ );
      }
    }
    return vec_corr_names;
  }
private:
  std::vector<std::string> vector_names_{};
  std::vector<size_t> harmonics_{};
};

template<typename... Column_t>
struct CorrelationAxes{
  std::vector<std::string> axes_columns{};
  std::vector<Qn::AxisD> axes{};
}

template<typename Column_t>
struct Weight{
  std::string weight_column{};
}

template<typename DF>
class CorrelationHandler{
public:
  CorrelationHandler(DF& df, size_t n_samples=100) : dataframe_(df), n_samples_(n_samples) {}
  
  template<typename Decorator_t, typename Weight_t typename... Axes_t>
  auto AddCorrelation( const Decorator_t& decorator, const Weight<Weight_t>& weight, const CorrelationAxes<Axes_t...>& axes, Qn::Stat::WeightType correlation_weight_type = Qn::Stat::WeightType::OBSERVABLE ) -> CorrelationHandler& {
    auto vec_corr_names = decorator( df );
    for( const auto& name : vec_corr_names ){
      auto vec_columns = std::vector< std::string > { name, weight.weight_column, "samples" };
      vec_columns.insert( vec_columns.end(), axes.begin(), axes.end() );
      if constexpr ( std::is_same_v<qvector_t, typename Decorator_t::First_t>  ) {
        result_ptrs_.emplace_back(
          dataframe_.Book< double, Weight_t, ROOT::RVec<ULong64_t>, Axes_t >( CorrelationHelper( axes.axes, n_samples_, correlation_weight_type ), vec_columns )
        );
      } else {
        result_ptrs_.emplace_back(
          dataframe_.Book< std::vector<double>, Weight_t, ROOT::RVec<ULong64_t>, Axes_t >( CorrelationHelper( axes.axes, n_samples_, correlation_weight_type ), vec_columns )
        );
      }
    }
    result_names_.insert( result_names_.end(), vec_corr_names.begin(), vec_corr_names.end() );
    return *this;
  }

  auto DumpCorrelations( TFile* file_out ){
    file_out->cd();
    std::for_each( result_ptrs_.begin(), result_ptrs_.end(), [this, i=0]( const auto& p ) mutable { p->Write( result_names_[i] ); } );
  }

private:
  DF& dataframe_;
  size_t n_samples_{};
  std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > > result_ptrs_{};
  std::vector< std::string > result_names_{};

};

#endif // CORRELATION_HELPER_H  