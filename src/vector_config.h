//
// Created by Misha on 4/6/2023.
//

#ifndef QNTOOLSINTERFACE_VECTOR_CONFIG_H
#define QNTOOLSINTERFACE_VECTOR_CONFIG_H

#include <Axis.hpp>
#include <CorrectionManager.hpp>

enum class VECTOR_TYPE{
  CHANNEL,
  TRACK
};

enum class NORMALIZATION{
  M,
  MAG,
  NONE,
};

enum class CORRECTION{
  PLAIN,
  RECENTERING,
  TWIST_RESCALING
};

enum class TWIST_RESCALING_METHOD{
  DOUBLE_HARMONIC,
  CORRELATION,
};


struct vector_cut;
struct histo1d;
struct histo2d;

struct vector_cut{
  std::string field;
  std::function<bool(double)> function;
  std::string description;
};

struct histo1d{
  Qn::AxisD axis;
  std::string weight{"Ones"};
};

struct histo2d{
  std::vector<Qn::AxisD> axes;
  std::string weight{"Ones"};
};


class VectorConfig {
public:
  VectorConfig(std::string name, std::string phi_field, std::string weight_field,
               VECTOR_TYPE type, NORMALIZATION normalization);
  VectorConfig(const VectorConfig&) = default;
  VectorConfig(VectorConfig&&) = default;
  VectorConfig& operator=(const VectorConfig&) = default;
  VectorConfig& operator=(VectorConfig&&) = default;
  ~VectorConfig() = default;

  const std::string &GetName() const {
    return name_;
  }

  const std::string &GetPhiField() const {
    return phi_field_;
  }

  const std::string &GetWeightField() const {
    return weight_field_;
  }

  VECTOR_TYPE GetType() const {
    return type_;
  }

  VectorConfig& SetHarmonicArray(const std::vector<int> &harmonic_array) {
    harmonic_array_ = harmonic_array;
    return *this;
  }
  VectorConfig& SetCorrectionAxes(const std::vector<Qn::AxisD> &correction_axes) {
    correction_axes_ = correction_axes;
    return *this;
  }
  VectorConfig& SetCorrections(const std::vector<CORRECTION> &corrections) {
    corrections_ = corrections;
    return *this;
  }
  VectorConfig& SetTwistRescalingMethod(TWIST_RESCALING_METHOD method){ twis_rescaling_method_ = method; return *this; }
  VectorConfig& SetTwistRescalingReference(std::array<std::string, 2> reference){ twist_rescaling_reference_ = reference; return *this; }
  VectorConfig& SetRecenteringWidthEqualization( bool value ){ recentering_width_equalization_ = value; return *this; }
  VectorConfig& SetApplyTwist( bool value ){ apply_twist_ = value; return *this; }
  VectorConfig& SetRescaling( bool value ){ apply_rescaling_ = value; return *this; }
  void AddCut( const std::string& field, const std::function<bool(double)>& function, const std::string& description );
  void AddHisto1D( const Qn::AxisD& axis, const std::string& weight= "Ones" );
  void AddHisto2D( const std::vector<Qn::AxisD>& axis, const std::string& weight= "Ones" );
  void Decorate( const std::shared_ptr<Qn::CorrectionManager>& man ) const;

private:
  std::string name_;
  std::string phi_field_;
  std::string weight_field_;
  VECTOR_TYPE type_;
  NORMALIZATION normalization_;
  TWIST_RESCALING_METHOD twis_rescaling_method_{TWIST_RESCALING_METHOD::DOUBLE_HARMONIC};
  std::array<std::string, 2> twist_rescaling_reference_{};
  std::vector<int> harmonic_array_{1};
  std::vector<CORRECTION> corrections_;
  bool recentering_width_equalization_{false};
  bool apply_twist_{true};
  bool apply_rescaling_{true};
  std::vector<Qn::AxisD> correction_axes_{};
  std::vector<vector_cut> cuts_{};
  std::vector<histo1d> vec_histo1d_{};
  std::vector<histo2d> vec_histo2d_{};
};


#endif //QNTOOLSINTERFACE_VECTOR_CONFIG_H
