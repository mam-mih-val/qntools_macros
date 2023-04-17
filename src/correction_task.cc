//
// Created by Misha on 2/1/2023.
//

#include "correction_task.h"

CorrectionTask::CorrectionTask(
        ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> data_frame,
        const std::string& outFilePath,
        const std::string& calibFilePath) : filtered_data_frame_{std::move(data_frame)} {
  out_file_ = TFile::Open(outFilePath.c_str(), "recreate");
  out_file_->cd();
  out_tree_ = new TTree("tree", "tree");
  qa_file_path_ = outFilePath;
  correction_manager_ = std::make_shared<Qn::CorrectionManager>();
  correction_manager_->SetCalibrationInputFileName(calibFilePath);
  correction_manager_->SetFillOutputTree(true);
  correction_manager_->SetFillCalibrationQA(true);
  correction_manager_->SetFillValidationQA(true);
  correction_manager_->ConnectOutputTree(out_tree_);
  varible_manager_.SetCorrectionManager(correction_manager_);
}

void CorrectionTask::InitVariables(){
  varible_manager_.InitVariables(filtered_data_frame_);
  varible_manager_.ConfigureRDF(filtered_data_frame_);
  varible_manager_.PrintVCLayout();
}

void CorrectionTask::SpecifyVariables(
        const std::regex& event_variable_pattern,
        const std::vector<std::regex>& channel_variable_patterns,
        const std::vector<std::regex>& track_variable_patterns){
  varible_manager_.SetEventVariablePattern(event_variable_pattern);
  varible_manager_.SetChannelVariablePatterns(channel_variable_patterns);
  varible_manager_.SetTrackVariablePatterns(track_variable_patterns);
}

void CorrectionTask::Run() {
  std::for_each( event_axes_.begin(), event_axes_.end(), [this](const Qn::AxisD& ax){ this->correction_manager_->AddCorrectionAxis(ax); } );
  for( const auto& conf : vector_configs_ ){
    conf.Decorate( correction_manager_ );
  }
  correction_manager_->InitializeOnNode();
  correction_manager_->SetCurrentRunName("test");
  filtered_data_frame_.Foreach(varible_manager_, {"rdfentry_", "event_variables", "channel_variables", "track_variables"});
  std::cout << std::endl;

  correction_manager_->Finalize();
  out_file_->cd();
  out_tree_->Write("tree");
  out_file_->Close();

  TFile qaFile(qa_file_path_.c_str(), "RECREATE");
  correction_manager_->GetCorrectionQAList()->Write("CorrectionQAHistograms", TObject::kSingleKey);
  correction_manager_->GetCorrectionList()->Write("CorrectionHistograms", TObject::kSingleKey);
  qaFile.Close();
}