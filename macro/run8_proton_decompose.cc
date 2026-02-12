#include "matrix.h"
#include "DataContainer.hpp"
#include "StatCollect.hpp"
#include "StatCalculate.hpp"
#include "QVector.hpp"
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

std::tuple< Matrix< Qn::DataContainerStatCalculate, 4 >, std::vector<Qn::DataContainerStatCalculate> > CorrectionMatrix( std::string str_vec_name, TFile* calib_file ){
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

  auto vec_recentering = std::vector< Qn::DataContainerStatCalculate >{ vec_c[0], vec_s[0], vec_c[1], vec_s[1] };

  return {invM, vec_recentering};
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

  auto [corrM_f1, rec_f1] = CorrectionMatrix(f1_name, calib_file.get());
  auto [corrM_f2, rec_f2] = CorrectionMatrix(f2_name, calib_file.get());
  auto [corrM_f3, rec_f3] = CorrectionMatrix(f3_name, calib_file.get());
  auto [corrM_tp, rec_tp] = CorrectionMatrix(tp_name, calib_file.get());
  auto [corrM_tn, rec_tn] = CorrectionMatrix(tn_name, calib_file.get());
  auto [corrM_p, rec_p] = CorrectionMatrix(proton_name, calib_file.get());

  const auto correction_generator = []( const Matrix< Qn::DataContainerStatCalculate, 4 >& corrM, const std::vector<Qn::DataContainerStatCalculate>& vec_rec ){
    return [corrM, vec_rec]( Qn::DataContainerQVector qvec, Double_t centrality ){
      auto new_qvec = qvec;
      auto corrM_c = corrM;
      auto rec_c = vec_rec;
      auto c_axis = corrM_c[0][0].GetAxis( "centrality" ); 
      auto c_bin = c_axis.FindBin( centrality );
      auto bin_lo = c_axis.GetLowerBinEdge( c_bin );
      auto bin_hi = c_axis.GetUpperBinEdge( c_bin );
      auto new_c_axis = Qn::AxisD{ "centrality", 1, bin_lo, bin_hi };
      std::cout << "here1" << std::endl;
      for( auto i = size_t{0}; i<4; ++i ){
      std::cout << "here1" << " dimensions = " << corrM_c[i][j].GetAxes().size() << std::endl;
        
        for( auto j = size_t{0}; j<4; ++i ){
          if( corrM_c[i][j].GetAxes().size() > 1 )
            corrM_c[i][j] = corrM_c[i][j].Select( new_c_axis );
          if( corrM_c[i][j].GetAxes().size() == 1 )
            corrM_c[i][j] = corrM_c[i][j].Rebin( new_c_axis );
        }
      }
      std::cout << "here2" << std::endl;
      for( auto& el : rec_c ) {
        if( el.GetAxes().size() > 1 )
            el = el.Select( new_c_axis );
          if( el.GetAxes().size() == 1 )
            el = el.Rebin( new_c_axis );
      }
      std::cout << "here3" << std::endl;
      for( auto i=size_t{0}; i<qvec.size(); ++i ){
        auto x1_old = qvec.At(i).x(1) - rec_c[0].At(i).Mean();
        auto y1_old = qvec.At(i).y(1) - rec_c[1].At(i).Mean();
        auto x2_old = qvec.At(i).x(2) - rec_c[2].At(i).Mean();
        auto y2_old = qvec.At(i).y(2) - rec_c[3].At(i).Mean();

        auto x1_new = x1_old * corrM_c[0][0].At(i).Mean() + y1_old * corrM_c[1][0].At(i).Mean() + x2_old * corrM_c[2][0].At(i).Mean() + y2_old * corrM_c[3][0].At(i).Mean();
        auto y1_new = x1_old * corrM_c[0][1].At(i).Mean() + y1_old * corrM_c[1][1].At(i).Mean() + x2_old * corrM_c[2][1].At(i).Mean() + y2_old * corrM_c[3][1].At(i).Mean();
        auto x2_new = x1_old * corrM_c[0][2].At(i).Mean() + y1_old * corrM_c[1][2].At(i).Mean() + x2_old * corrM_c[2][2].At(i).Mean() + y2_old * corrM_c[3][2].At(i).Mean();
        auto y2_new = x1_old * corrM_c[0][3].At(i).Mean() + y1_old * corrM_c[1][3].At(i).Mean() + x2_old * corrM_c[2][3].At(i).Mean() + y2_old * corrM_c[3][3].At(i).Mean();
        
        new_qvec.At(i).SetQ( 1, x1_new, y1_new );
        new_qvec.At(i).SetQ( 2, x2_new, y2_new );
      }
      std::cout << "here4" << std::endl;

      return new_qvec;
    };
  };

  auto d = ROOT::RDataFrame( "tree", in_file_name );
  auto dd = d
    .Define("F1_DECOMPOSED", correction_generator(corrM_f1, rec_f1), { f1_name, "centrality" } )
    .Define("F2_DECOMPOSED", correction_generator(corrM_f2, rec_f2), { f2_name, "centrality" } )
    .Define("F3_DECOMPOSED", correction_generator(corrM_f3, rec_f3), { f3_name, "centrality" } )

    .Define("Tpos_DECOMPOSED", correction_generator(corrM_tp, rec_tp), { tp_name, "centrality" } )
    .Define("Tneg_DECOMPOSED", correction_generator(corrM_tn, rec_tn), { tn_name, "centrality" } )
    
    .Define("proton_DECOMPOSED", correction_generator(corrM_p, rec_p), { proton_name, "centrality" } )
  ;

  dd.Snapshot("tree", "decomposed_out.root", std::vector<std::string>{ "F1_PLAIN", "F1_DECOMPOSED" } );
  std::cout << "here5" << std::endl;


}