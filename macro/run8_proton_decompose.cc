#include "matrix.h"
#include "DataContainer.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"
#include "QVector.hpp"
#include <array>
#include <cassert>
#include <cstddef>
#include <string>

Matrix< Qn::DataContainerStatCalculate, 4 > CorrectionMatrix( std::string str_vec_name, TFile* calib_file ){
  Qn::DataContainerStatCollect* tmp{nullptr};

  auto vec_c = std::vector<Qn::DataContainerStatCalculate>{};
  auto vec_s = std::vector<Qn::DataContainerStatCalculate>{};
  vec_c.reserve(4);
  vec_s.reserve(4);
  for( auto i=size_t{0}; i<4; ++i ){
    auto corr_name = str_vec_name+".x"+std::to_string(i+1)+"centrality"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_c.emplace_back( *tmp );

    corr_name = str_vec_name+".y"+std::to_string(i+1)+"centrality"s;
    calib_file->GetObject( corr_name.c_str(), tmp );
    assert(tmp);
    vec_s.emplace_back( *tmp );
  }

  auto M = Matrix< Qn::DataContainerStatCalculate, 4 >{};
  M[0] = std::array<Qn::DataContainerStatCalculate, 4>{ vec_c[1] + 1, vec_s[1], vec_c[0] + vec_c[2], vec_s[0] + vec_s[2] };
  M[1] = std::array<Qn::DataContainerStatCalculate, 4>{ vec_s[1], 1 - vec_c[1], vec_s[2] - vec_c[0], vec_c[0] - vec_c[2] };
  M[2] = std::array<Qn::DataContainerStatCalculate, 4>{ vec_c[0] + vec_c[2], vec_s[2] - vec_s[0], vec_c[3] + 1, vec_s[3] };
  M[3] = std::array<Qn::DataContainerStatCalculate, 4>{ vec_s[0] + vec_s[2], vec_c[0] - vec_c[2], vec_s[3], 1 - vec_c[3] };

  auto invM = Inverse(M);

  return invM;
}

void run8_proton_decompose(std::string in_file_name, std::string in_calib_file){

  const std::string f1_name {"F1_PLAIN"};
  const std::string f2_name {"F2_PLAIN"};
  const std::string f3_name {"F3_PLAIN"};
  const std::string tp_name {"Tpos_PLAIN"};
  const std::string tn_name {"Tneg_PLAIN"};
  const std::string proton_name {"proton_PLAIN"};

  auto calib_file = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( in_calib_file.c_str(), "READ" ), [](auto f){f ->Close(); } };
  assert(calib_file);

  auto corrM_f1 = CorrectionMatrix(f1_name, calib_file.get());
  auto corrM_f2 = CorrectionMatrix(f2_name, calib_file.get());
  auto corrM_f3 = CorrectionMatrix(f3_name, calib_file.get());
  auto corrM_tp = CorrectionMatrix(tp_name, calib_file.get());
  auto corrM_tn = CorrectionMatrix(tn_name, calib_file.get());
  auto corrM_p = CorrectionMatrix(proton_name, calib_file.get());

  const auto correction_generator = []( const Matrix< Qn::DataContainerStatCalculate, 4 >& corrM ){
    return [corrM&]( Double_t centrality, Qn::DataContainerQVector qvec ){
      auto new_qvec = qvec;
      auto corrM_c = corrM;
      auto c_axis = corrM_c[0][0].GetAxis( "centrality" ).FindBin( centrality ); 
      auto c_bin = c_axis.FindBin( centrality );
      auto bin_lo = c_axis.GetLowerBinEdge( c_bin );
      auto bin_hi = c_axis.GetUpperBinEdge( c_bin );
      auto new_c_axis = Qn::AxisD{ "centrality", 1, bin_lo, bin_hi };
      for( auto i = size_t{0}; i<4; ++i ){
        for( auto j = size_t{0}; j<4; ++i ){
          auto c_axis = corrM_c[i][j].Select( new_c_axis );
        }
      }

      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        auto x1_new = qvec.At(i).x(1) * corrM_c[0][0].At(i) + qvec.At(i).y(1) * corrM_c[1][0].At(i) + qvec.At(i).x(2) * corrM_c[2][0].At(i) + qvec.At(i).y(2) * corrM_c[3][0].At(i);
        auto y1_new = qvec.At(i).x(1) * corrM_c[0][1].At(i) + qvec.At(i).y(1) * corrM_c[1][1].At(i) + qvec.At(i).x(2) * corrM_c[2][1].At(i) + qvec.At(i).y(2) * corrM_c[3][1].At(i);
        auto x2_new = qvec.At(i).x(1) * corrM_c[0][2].At(i) + qvec.At(i).y(1) * corrM_c[1][2].At(i) + qvec.At(i).x(2) * corrM_c[2][2].At(i) + qvec.At(i).y(2) * corrM_c[3][2].At(i);
        auto y2_new = qvec.At(i).x(1) * corrM_c[0][3].At(i) + qvec.At(i).y(1) * corrM_c[1][3].At(i) + qvec.At(i).x(2) * corrM_c[2][3].At(i) + qvec.At(i).y(2) * corrM_c[3][3].At(i);
        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_new, y2_new );
      }

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Define("F1_DECOMPOSED", correction_generator(corrM_f1), { f1_name, "centrality" } )
    .Define("F2_DECOMPOSED", correction_generator(corrM_f2), { f2_name, "centrality" } )
    .Define("F3_DECOMPOSED", correction_generator(corrM_f3), { f3_name, "centrality" } )

    .Define("Tpos_DECOMPOSED", correction_generator(corrM_tp), { tp_name, "centrality" } )
    .Define("Tneg_DECOMPOSED", correction_generator(corrM_tn), { tn_name, "centrality" } )
    
    .Define("proton_DECOMPOSED", correction_generator(corrM_tp), { corrM_p, "centrality" } )
  ;

  dd.Snapshot("tree", "decomposed_out.root", dd.GetColumnNames() );

}