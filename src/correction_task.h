//
// Created by Misha on 2/1/2023.
//

#ifndef QNTOOLSINTERFACE_CORRECTION_TASK_H
#define QNTOOLSINTERFACE_CORRECTION_TASK_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <regex>

#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <CorrectionManager.hpp>

#include "variable_manager.h"
#include "vector_config.h"

class CorrectionTask {
public:
  CorrectionTask(
          ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> data_frame,
          const std::string& outFilePath,
          const std::string& calibFilePath);

  ~CorrectionTask() = default;

  void InitVariables();

  [[nodiscard]] const std::shared_ptr<Qn::CorrectionManager> &GetCorrectionManager() const {
    return correction_manager_;
  }

  void Run();

  void SetEventVariables(const std::regex& event_variable_pattern){
    varible_manager_.SetEventVariablePattern(event_variable_pattern);
  }

  void SetChannelVariables(const std::vector<std::regex>& channel_variable_patterns){
    varible_manager_.SetChannelVariablePatterns(channel_variable_patterns);
  }

  void SetTrackVariables(const std::vector<std::regex>& track_variable_patterns){
    varible_manager_.SetTrackVariablePatterns(track_variable_patterns);
  }

  void SpecifyVariables(
          const std::regex& event_variable_pattern,
          const std::vector<std::regex>& channel_variable_patterns,
          const std::vector<std::regex>& track_variable_patterns);

  void AddVector(const VectorConfig& vector_config ){
    vector_configs_.push_back(vector_config);
  }

  void SetVectorConfigs(const std::vector<VectorConfig> &vector_configs) {
    vector_configs_ = vector_configs;
  }

  void SetEventAxes(const std::vector<Qn::AxisD> &event_axes) {
    event_axes_ = event_axes;
  }
  void AddEventAxis(const Qn::AxisD &event_axis){
    event_axes_.push_back(event_axis);
  }

private:
  TFile *out_file_;
  TTree *out_tree_;
  std::string qa_file_path_;
  VariableManager varible_manager_;
  std::vector<VectorConfig> vector_configs_;
  std::shared_ptr<Qn::CorrectionManager> correction_manager_;
  std::vector<Qn::AxisD > event_axes_;
  ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> filtered_data_frame_;
};


#endif //QNTOOLSINTERFACE_CORRECTION_TASK_H
