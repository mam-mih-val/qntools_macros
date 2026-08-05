#ifndef CORRELATION_HELPER_H
#define CORRELATION_HELPER_H

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include <vector>
#include <cmath>
#include <type_traits>
#include <DataContainer.hpp>
#include <Axis.hpp>

class CorrelationHelper :  public ROOT::Detail::RDF::RActionImpl<CorrelationHelper>{
public:
   /// This type is a requirement for every helper.
  using Result_t = Qn::DataContainerStatCalculate;

  CorrelationHelper( std::vector< Qn::AxisD > vec_axes ) : 
  final_result_{ new Qn::DataContainerStatCollect(vec_axes) },
  thread_results_{ std::vector<Qn::DataContainerStatCollect>{} }{
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
  void Exec(unsigned int slot, ColumnTypes... values){ Execute(slot, values...); }

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
  template <typename... ColumnTypes>
  void Execute(unsigned int slot, std::vector<double> vec_val, std::vector<double> vec_weights, ROOT::RVec<ULong64_t> vec_samples, ColumnTypes... coordinates){
    for( auto i=size_t{}; i<vec_val.size(); ++i ){
      auto val = vec_val.at(i);
      auto weight = vec_weights.at(i);
      auto coord = FormCoordinates( i, coordinates... );
      auto lin_idx = thread_results_[slot].GetLinearIndex( coord );
      thread_results_[slot][lin_idx].Fill( val, weight, vec_samples );
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

#endif // CORRELATION_HELPER_H