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

  void ConfigureRDF( ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>& filtered_data_frame );

  void InitVariables( ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>& filtered_data_frame );

  [[nodiscard]] const std::regex &GetEventVariablePattern() const {
    return event_variable_pattern_;
  }

  [[nodiscard]] const std::vector<std::regex> &GetChannelVariablePatterns() const {
    return channel_variable_patterns_;
  }

  [[nodiscard]] const std::vector<std::regex> &GetTrackVariablePatterns() const {
    return track_variable_patterns_;
  }

  void operator()( ULong64_t event_id,
                   const std::vector<float>& event_variables,
                   const std::vector<std::vector<std::vector<float>>>& module_variables,
                   const std::vector<std::vector<std::vector<float>>>& track_variables );

  void PrintVCLayout();

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
