//
// Created by Misha on 4/6/2023.
//

#include "vector_config.h"

#include <type_traits>
#include <utility>

VectorConfig::VectorConfig(std::string name, std::string phi_field, std::string weight_field,
                           VECTOR_TYPE type, NORMALIZATION normalization) : name_(std::move(name)), phi_field_(std::move(phi_field)), weight_field_(std::move(weight_field)),
                                               type_(type), normalization_(normalization) {}

void VectorConfig::AddCut( const std::string& field, const std::function<bool(double)>& function, const std::string& description ){
  cuts_.push_back( {field, function, description} );
}

void VectorConfig::AddHisto1D( const Qn::AxisD& axis, const std::string& weight ){
  vec_histo1d_.push_back( {axis, weight} );
}

void VectorConfig::AddHisto2D( const std::vector<Qn::AxisD>& axes, const std::string& weight ){
  vec_histo2d_.push_back( {axes, weight} );
}

void VectorConfig::Decorate(const std::shared_ptr<Qn::CorrectionManager>& man) const {

  std::map<VECTOR_TYPE, Qn::DetectorType> detector_types{
          { VECTOR_TYPE::TRACK, Qn::DetectorType::TRACK },
          { VECTOR_TYPE::CHANNEL, Qn::DetectorType::CHANNEL }
  };
  std::map<NORMALIZATION, Qn::QVector::Normalization> normalizations{
          { NORMALIZATION::M, Qn::QVector::Normalization::M },
          { NORMALIZATION::MAG, Qn::QVector::Normalization::MAGNITUDE },
          { NORMALIZATION::NONE, Qn::QVector::Normalization::NONE },
  };

  auto plain=Qn::QVector::CorrectionStep::PLAIN;
  auto recentered=Qn::QVector::CorrectionStep::RECENTERED;
  auto twist=Qn::QVector::CorrectionStep::TWIST;
  auto rescaled=Qn::QVector::CorrectionStep::RESCALED;

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(recentering_width_equalization_);

  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(apply_rescaling_);
  twistRescale.SetApplyTwist(apply_twist_);
  if( twis_rescaling_method_ == TWIST_RESCALING_METHOD::DOUBLE_HARMONIC )
    twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
  else{
    twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::CORRELATIONS );
    twistRescale.SetReferenceConfigurationsForTwistAndRescale( twist_rescaling_reference_.at(0), twist_rescaling_reference_.at(1) );
  }

  std::bitset<Qn::QVector::kmaxharmonics> harmonics_bitset{};
  for (int harm : harmonic_array_) {
    harmonics_bitset[harm - 1] = true;
  }
  man->AddDetector( name_,
                    detector_types.at(type_),
                    phi_field_, weight_field_,
                    correction_axes_,
                    harmonics_bitset,
                    normalizations.at(normalization_) );

  std::vector<Qn::QVector::CorrectionStep> correction_output;
  for( const auto& corr : corrections_ ){
    if( corr == CORRECTION::PLAIN ){
      correction_output.push_back( plain );
    }
    if( corr == CORRECTION::RECENTERING ){
      man->AddCorrectionOnQnVector( name_, recentering );
      correction_output.push_back( recentered );
    }
    if( corr == CORRECTION::TWIST_RESCALING ){
      man->AddCorrectionOnQnVector( name_, twistRescale );
      correction_output.push_back( twist );
      correction_output.push_back( rescaled );
    }
  }
  man->SetOutputQVectors( name_, correction_output );
  for( const auto & cut : cuts_){
    man->AddCutOnDetector( name_, { cut.field.c_str() }, cut.function, cut.description );
  }
  for( const auto& hist : vec_histo1d_ ){
    man->AddHisto1D( name_, hist.axis, hist.weight );
  }
  for( const auto& hist : vec_histo2d_ ){
    man->AddHisto2D( name_, hist.axes, hist.weight );
  }
}
