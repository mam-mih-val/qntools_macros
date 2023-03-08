//
// Created by Misha on 3/7/2023.
//

#ifndef QNTOOLSINTERFACE_VARIABLE_MANAGER_H
#define QNTOOLSINTERFACE_VARIABLE_MANAGER_H

#include <vector>
#include <string>
#include <regex>

#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <CorrectionManager.hpp>


class VariableManager {
public:
  VariableManager() = default;
  ~VariableManager() = default;

  void SetCorrectionManager(const std::shared_ptr<Qn::CorrectionManager> &correction_manager) {
    correction_manager_ = correction_manager;
  }

  void SetEventVariablePattern(const std::regex &event_variable_pattern) {
    event_variable_pattern_ = event_variable_pattern;
  }

  void SetChannelVariablePatterns(const std::vector<std::regex> &channel_variable_patterns) {
    channel_variable_patterns_ = channel_variable_patterns;
  }

  void SetTrackVariablePatterns(const std::vector<std::regex> &track_variable_patterns) {
    track_variable_patterns_ = track_variable_patterns;
  }
  void InitVariables( ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>& filtered_data_frame ){
    channel_variable_names_ = std::vector<std::vector<std::string>>( channel_variable_patterns_.size() );
    channel_variable_sizes_ = std::vector<std::vector<int>>( channel_variable_patterns_.size() );
    channel_variable_positions_ = std::vector<std::vector<int>>( channel_variable_patterns_.size() );
    track_variable_names_ = std::vector<std::vector<std::string>>( track_variable_patterns_.size() );
    track_variable_positions_ = std::vector<std::vector<int>>( track_variable_patterns_.size() );
    for(const auto& name : filtered_data_frame.GetColumnNames() ){
      auto column_type = filtered_data_frame.GetColumnType(name);
      if( std::regex_match( name, event_variable_pattern_ ) ){
        event_variable_names_.push_back( name );
        continue;
      }
      for( size_t i=0; i<channel_variable_patterns_.size(); ++i ){
        if( !std::regex_match( name, channel_variable_patterns_.at(i) ) )
          continue;
        auto size=static_cast<int>(*filtered_data_frame.Range(0, 0).Define("n", name + ".size()").Mean("n"));
        channel_variable_names_.at(i).push_back(name);
        channel_variable_sizes_.at(i).push_back(size);
      }
      for( size_t i=0; i<track_variable_patterns_.size(); ++i ){
        if( !std::regex_match( name, track_variable_patterns_.at(i) ) )
          continue;
        track_variable_names_.at(i).push_back(name);
      }
    }
    int idx = 0;
    for( const auto& name : event_variable_names_ ){
      correction_manager_->AddVariable(name, idx, 1 );
      correction_manager_->AddEventVariable( name );
      event_variable_positions_.push_back(idx);
      idx++;
    }
    for( size_t i=0; i<channel_variable_names_.size(); ++i ){
      auto names = channel_variable_names_.at(i);
      auto sizes = channel_variable_sizes_.at(i);
      for( size_t j=0; j<names.size(); j++ ){
        std::cout << names.at(j) << " " << sizes.at(j) << " " << idx << std::endl;
        correction_manager_->AddVariable(names.at(j), idx, sizes.at(j) );
        channel_variable_positions_.back().push_back(idx);
        idx+=sizes.at(j);
      }
    }
    for( size_t i=0; i<track_variable_names_.size(); ++i ){
      auto names = track_variable_names_.at(i);
      channel_variable_sizes_.emplace_back();
      for(const auto & name : names){
        correction_manager_->AddVariable(name, idx, 1 );
        track_variable_positions_.at(i).push_back(idx);
        idx++;
      }
    }
    correction_manager_->AddVariable("track_type", idx, 1 );
    track_type_position_=idx;
  }

  void ConfigureRDF( ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>& filtered_data_frame ){
    std::string event_variable_code{ "std::vector<float> event_variables;\n" };
    for( const auto& name : event_variable_names_ ){
      std::string line = "event_variables.push_back("+name+");\n";
      event_variable_code+=line;
    }
    event_variable_code+="return event_variables;\n";
    std::cout << event_variable_code << std::endl;

    std::string channel_variables_code{ "std::vector< std::vector<std::vector<float>> > channel_variables;\n" };
    for(const auto& channels : channel_variable_names_){
      channel_variables_code+="channel_variables.emplace_back();\n";
      for( const auto& name : channels ){
        channel_variables_code+="channel_variables.back().emplace_back();\n";
        channel_variables_code+="for(const auto& val : "+name+") channel_variables.back().back().push_back(val);\n";
      }
    }
    channel_variables_code+="return channel_variables;\n";
    std::cout << channel_variables_code << std::endl;

    std::string track_variable_code{ "std::vector< std::vector<std::vector<float>> > track_variables;\n" };
    for(const auto& track : track_variable_names_){
      track_variable_code+="track_variables.emplace_back();\n";
      track_variable_code+="for(size_t i=0; i< " + track.front() + ".size(); ++i){\n";
      track_variable_code+="\ttrack_variables.back().emplace_back();\n";
      for( const auto& name : track ){
        track_variable_code+="\ttrack_variables.back().back().push_back("+name+".at(i));\n";
      }
      track_variable_code+="}\n";
    }
    track_variable_code+="return track_variables;\n";
    std::cout << track_variable_code << std::endl;
    filtered_data_frame = filtered_data_frame
            .Define("event_variables", event_variable_code)
            .Define("channel_variables", channel_variables_code)
            .Define("track_variables", track_variable_code);
  }

  void operator()( const ULong64_t event_id,
                   const std::vector<float>& event_variables,
                   const std::vector<std::vector<std::vector<float>>>& module_variables,
                   const std::vector<std::vector<std::vector<float>>>& track_variables ){
    correction_manager_->Reset();
    double *var_container = correction_manager_->GetVariableContainer();
    for( int i=0; i<event_variables.size(); ++i ){
      var_container[ event_variable_positions_.at(i) ] = event_variables.at(i);
    }
    for( int i=0; i<module_variables.size(); ++i ){
      for( int j=0; j<module_variables.at(i).size(); ++j ){
        for( int k=0; k<module_variables.at(i).at(j).size(); ++k ) {
          var_container[channel_variable_positions_.at(i).at(j) + k] = module_variables.at(i).at(j).at(k);
        }
      }
    }
    correction_manager_->ProcessEvent();
    correction_manager_->FillChannelDetectors();
    for( int i=0; i<track_variables.size(); ++i ){
      auto track_type=i;
      for( int j=0; j<track_variables.at(i).size(); ++j ){
        for( int k=0; k<track_variables.at(i).at(j).size(); ++k )
          var_container[ track_variable_positions_.at(i).at(k) ] = track_variables.at(i).at(j).at(k);
        var_container[ track_type_position_ ] = track_type;
        correction_manager_->FillTrackingDetectors();
      }
    }
    correction_manager_->ProcessCorrections();
  }
  void PrintVCLayout(){
    std::cout << "Event Variables" << std::endl;
    for( int i=0; i<event_variable_names_.size(); ++i )
      std::cout << "name: " << event_variable_names_.at(i) << " id: " << event_variable_positions_.at(i) << std::endl;

    std::cout << "Channel Variables" << std::endl;
    for( int i=0; i<channel_variable_names_.size(); ++i )
      for( int j=0; j<channel_variable_names_.at(i).size(); ++j )
        std::cout << "name: " << channel_variable_names_.at(i).at(j) << " id: " << channel_variable_positions_.at(i).at(j) << "+" << channel_variable_sizes_.at(i).at(j) << std::endl;

    std::cout << "Track Variables" << std::endl;
    for( int i=0; i<track_variable_names_.size(); ++i )
      for( int j=0; j<track_variable_names_.at(i).size(); ++j )
        std::cout << "name: " << track_variable_names_.at(i).at(j) << " id: " << track_variable_positions_.at(i).at(j) << std::endl;
  }

private:
  std::shared_ptr<Qn::CorrectionManager> correction_manager_;

  std::regex event_variable_pattern_;
  std::vector<std::regex> channel_variable_patterns_;
  std::vector<std::regex> track_variable_patterns_;

  std::vector<std::string> event_variable_names_;
  std::vector<int> event_variable_positions_;

  std::vector<std::vector<std::string>> channel_variable_names_;
  std::vector<std::vector<int>> channel_variable_sizes_;
  std::vector<std::vector<int>> channel_variable_positions_;

  std::vector<std::vector<std::string>> track_variable_names_;
  std::vector<std::vector<int>> track_variable_positions_;
  int track_type_position_=-1;
};


#endif //QNTOOLSINTERFACE_VARIABLE_MANAGER_H
