#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <vector>
#include <functional>

#include <DataContainer.hpp>
#include <QnDataFrame.hpp>
#include "correlation_helper.h"
#include "bmn_env.h"
#include "vector_generators.h"
#include "corrections.h"

void eigenvalues( std::string str_calib_file ){

  std::cout << "starting execution" << std::endl;

  constexpr size_t NHARM = 6;
  constexpr size_t Q_NHARM = 1;

  const auto l = double{5e-2};
  auto harmonics = std::vector<size_t>(NHARM);
  auto harmonics_q = std::vector<size_t>(Q_NHARM);
  std::iota( harmonics.begin(), harmonics.end(), 1 );
  std::iota( harmonics_q.begin(), harmonics_q.end(), 1 );

  auto calib_file = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( str_calib_file.c_str(), "READ"), [](auto f){ f->Close(); } };
  auto [vec_pca_mean, vec_pca_cov] = ReadMeanCov<NHARM>("proton", calib_file.get());
  auto p_pca_container = MakeCorrectionContainer<NHARM>( vec_pca_mean, vec_pca_cov, PrincipalComponents<NHARM>{}, 0.95 );
  auto p_pca_builder = CorrectorBuilder<NHARM>( p_pca_container );
  
  auto [vec_f1_mean, vec_f1_cov] = ReadMeanCov<2*Q_NHARM>("F1", calib_file.get());
  auto f1_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f1_mean, vec_f1_cov, TwistRescale<Q_NHARM>{}, l );
  auto f1_corr_builder = CorrectorBuilder<Q_NHARM>( f1_correction_container );

  auto [vec_f2_mean, vec_f2_cov] = ReadMeanCov<2*Q_NHARM>("F2", calib_file.get());
  auto f2_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f2_mean, vec_f2_cov, TwistRescale<Q_NHARM>{}, l );
  auto f2_corr_builder = CorrectorBuilder<Q_NHARM>( f2_correction_container );

  auto [vec_f3_mean, vec_f3_cov] = ReadMeanCov<2*Q_NHARM>("F3", calib_file.get());
  auto f3_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f3_mean, vec_f3_cov, TwistRescale<Q_NHARM>{}, l );
  auto f3_corr_builder = CorrectorBuilder<Q_NHARM>( f3_correction_container );
  
  auto [vec_f4_mean, vec_f4_cov] = ReadMeanCov<2*Q_NHARM>("F4", calib_file.get());
  auto f4_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f4_mean, vec_f4_cov, TwistRescale<Q_NHARM>{}, l );
  auto f4_corr_builder = CorrectorBuilder<Q_NHARM>( f4_correction_container );
  
  auto [vec_tp_mean, vec_tp_cov] = ReadMeanCov<2*Q_NHARM>("Tpos", calib_file.get());
  auto tp_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tp_mean, vec_tp_cov, TwistRescale<Q_NHARM>{}, l );
  auto tp_corr_builder = CorrectorBuilder<Q_NHARM>( tp_correction_container );
  
  auto [vec_tn_mean, vec_tn_cov] = ReadMeanCov<2*Q_NHARM>("Tneg", calib_file.get());
  auto tn_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tn_mean, vec_tn_cov, TwistRescale<Q_NHARM>{}, l );
  auto tn_corr_builder = CorrectorBuilder<Q_NHARM>( tn_correction_container );
}