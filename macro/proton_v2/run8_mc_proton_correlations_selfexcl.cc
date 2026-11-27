//
// Created by Misha on 3/7/2023.
//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <vector>
#include <functional>

#include <TFileCollection.h>
#include <THashList.h>

#include <DataContainer.hpp>
#include <QnDataFrame.hpp>
#include "correlation_helper.h"
#include "bmn_env.h"
#include "vector_generators.h"
#include "corrections.h"

namespace {

// SplitMix64 finalizer.  Used to build deterministic per-track random
// numbers for the proton-subsampling scan.
inline std::uint64_t mix_hash( std::uint64_t x ){
  x += 0x9e3779b97f4a7c15ULL;
  x = ( x ^ ( x >> 30 ) ) * 0xbf58476d1ce4e5b9ULL;
  x = ( x ^ ( x >> 27 ) ) * 0x94d049bb133111ebULL;
  return x ^ ( x >> 31 );
}

inline double uniform_from_hash( std::uint64_t x ){
  return static_cast<double>( mix_hash( x ) >> 11 ) * 0x1.0p-53;
}

}  // namespace

void run8_mc_proton_correlations_selfexcl( std::string list, std::string str_effieciency_file, std::string str_calib_file ){

  std::cout << "starting execution" << std::endl;

  constexpr size_t NHARM = 8;
  constexpr size_t Q_NHARM = 1;

  const auto l = double{ 5e-2 };
  auto harmonics = std::vector<size_t>(NHARM);
  auto harmonics_q = std::vector<size_t>(Q_NHARM);
  
  std::iota( harmonics.begin(), harmonics.end(), 1 );
  std::iota( harmonics_q.begin(), harmonics_q.end(), 1 );

  auto f1_mod = std::vector<int>(10);
  std::iota( f1_mod.begin(), f1_mod.end(), 45 );
  auto f2_mod = std::vector<int>(34);
  std::iota( f2_mod.begin(), f2_mod.end(), 1 );
  auto f3_mod = std::vector<int>(10);
  std::iota( f3_mod.begin(), f3_mod.end(), 35 );

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("h2_efficiency_2212_good", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No was found in file " << str_effieciency_file << "\n";

  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  
  auto dd = GenerateBmnExtendedTreeMC(d, efficiency_histo);

  auto sampled_d = Qn::Correlation::Resample(dd, 100);

  sampled_d = sampled_d.Define( "F1w", fhcal_weight_generator(f1_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F2w", fhcal_weight_generator(f2_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F3w", fhcal_weight_generator(f3_mod), { "fhcalModId", "fhcalModE" } );

  sampled_d = sampled_d.Define( "F1phi", subevent_centered_phi(f1_mod), { "fhcalModId", "fhcalModX", "fhcalModY" } );
  sampled_d = sampled_d.Define( "F2phi", subevent_centered_phi(f2_mod), { "fhcalModId", "fhcalModX", "fhcalModY" } );
  sampled_d = sampled_d.Define( "F3phi", subevent_centered_phi(f3_mod), { "fhcalModId", "fhcalModX", "fhcalModY" } );

  // Q-vector from all protons, weighted with pT.
  sampled_d = sampled_d.Define( "trProtonPtW",
    []( const ROOT::VecOps::RVec<float>& vec_pT, const std::vector<int>& vec_is_proton ){
      std::vector<double> vec_w;
      const auto n = std::min( vec_pT.size(), vec_is_proton.size() );
      vec_w.reserve( n );
      for( size_t i=0; i<n; ++i )
        vec_w.push_back( vec_pT[i] * vec_is_proton[i] );
      return vec_w;
    }, { "trPt", "trIsProton" } );

  DefineVector( sampled_d, "ini_F1", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"F1phi", "F1w"} );
  DefineVector( sampled_d, "ini_F2", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"F2phi", "F2w"} );
  DefineVector( sampled_d, "ini_F3", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"F3phi", "F3w"} );
  
  DefineVector( sampled_d, "ini_Tpos", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"trPhi", "trTposW"} );
  DefineVector( sampled_d, "ini_Tneg", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"trPhi", "trTnegW"} );

  DefineVector(sampled_d, "ini_proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );

  // For each track, the pT-weighted average of all other tracks:
  //   P_{-i} = (W P - w_i u_i)/(W - w_i).
  DefineVector( sampled_d, "P_excl",
    u_vector_excluding_self<uvector_t, std::vector<double>>(harmonics),
    std::vector<std::string>{"ini_proton", "trProtonPtW"} );

  sampled_d = sampled_d.Define( "P_excl_sq",
    []( const uvector_t& vec_p ){
      std::vector<double> vec_sq;
      vec_sq.reserve( vec_p.size() );
      for( const auto& q : vec_p ){
        const auto it = q.find( 1 );
        if( it != q.end() )
          vec_sq.push_back( it->second.x * it->second.x + it->second.y * it->second.y );
        else
          vec_sq.push_back( 0. );
      }
      return vec_sq;
    }, { "P_excl" } );
  DefineVector( sampled_d, "psi_rp", psi_rp_vector< double >( harmonics ), std::vector<std::string>{"psiRP"s} );

  auto calib_file = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( str_calib_file.c_str(), "READ"), [](auto f){ f->Close(); } };
  auto [vec_pca_mean, vec_pca_cov] = ReadMeanCov<NHARM>("proton", calib_file.get());
  auto p_pca_container = MakeCorrectionContainer<NHARM>( vec_pca_mean, vec_pca_cov, PrincipalComponents<NHARM>{}, 0.99 );
  auto p_pca_builder = CorrectorBuilder<NHARM>( p_pca_container );
  sampled_d = sampled_d.Define( "pca_proton", p_pca_builder.IssueUVectorCorrector<uvector_t, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >(), { "ini_proton", "centrality", "trProtonY", "trPt" } );

  // Alias of the corrected u-vector for the centrality-only correlations.
  auto uvector_copy = []( const uvector_t& vec_u ){ return vec_u; };
  sampled_d = sampled_d.Define( "pca_proton_c", uvector_copy, { "pca_proton" } );

  auto [vec_f1_mean, vec_f1_cov] = ReadMeanCov<2*Q_NHARM>("F1", calib_file.get());
  auto f1_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f1_mean, vec_f1_cov, TwistRescale<Q_NHARM>{}, 0.005 );
  auto f1_corr_builder = CorrectorBuilder<Q_NHARM>( f1_correction_container );
  sampled_d = sampled_d.Define( "F1", f1_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F1", "centrality" } );

  auto [vec_f2_mean, vec_f2_cov] = ReadMeanCov<2*Q_NHARM>("F2", calib_file.get());
  auto f2_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f2_mean, vec_f2_cov, TwistRescale<Q_NHARM>{}, 0.005 );
  auto f2_corr_builder = CorrectorBuilder<Q_NHARM>( f2_correction_container );
  sampled_d = sampled_d.Define( "F2", f2_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F2", "centrality" } );

  auto [vec_f3_mean, vec_f3_cov] = ReadMeanCov<2*Q_NHARM>("F3", calib_file.get());
  auto f3_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f3_mean, vec_f3_cov, TwistRescale<Q_NHARM>{}, 0.005 );
  auto f3_corr_builder = CorrectorBuilder<Q_NHARM>( f3_correction_container );
  sampled_d = sampled_d.Define( "F3", f3_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F3", "centrality" } );

  auto [vec_tp_mean, vec_tp_cov] = ReadMeanCov<2*Q_NHARM>("Tpos", calib_file.get());
  auto tp_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tp_mean, vec_tp_cov, TwistRescale<Q_NHARM>{}, 0.005 );
  auto tp_corr_builder = CorrectorBuilder<Q_NHARM>( tp_correction_container );
  sampled_d = sampled_d.Define( "Tpos", tp_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tpos", "centrality" } );

  auto [vec_tn_mean, vec_tn_cov] = ReadMeanCov<2*Q_NHARM>("Tneg", calib_file.get());
  auto tn_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tn_mean, vec_tn_cov, TwistRescale<Q_NHARM>{}, 0.005 );
  auto tn_corr_builder = CorrectorBuilder<Q_NHARM>( tn_correction_container );
  sampled_d = sampled_d.Define( "Tneg", tn_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tneg", "centrality" } );

  auto qvector_axes = CorrelationAxes<float>{
    .axes_columns={ "centrality" },
    .axes=std::vector<Qn::AxisD>{
      Qn::AxisD{ "centrality", 6, 0, 60 },
    }
  };

  auto proton_weight = Weight<std::vector<double>>{ "trProtonWeight" };
  auto qvector_weight = Weight<double>{ "One" };

  auto begin = std::chrono::steady_clock::now();
  
  auto handler = CorrelationHandler{ sampled_d, 100 };
  handler

    


    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "F2"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "F3"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "F3"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"Tpos", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CovarianceDecorator<qvector_t>{ "Tpos", { 1 } }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CovarianceDecorator<qvector_t>{ "Tneg", { 1 } }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    // Centrality-only correlations for the self-exclusion solve.
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"pca_proton_c", "F1"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"pca_proton_c", "F2"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"pca_proton_c", "F3"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"pca_proton_c", "psi_rp"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, uvector_t>{ std::vector<std::string>{"pca_proton_c", "P_excl"}, {1, 1} }, proton_weight, qvector_axes)

    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"P_excl", "F1"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"P_excl", "F2"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"P_excl", "F3"}, {1, 1} }, proton_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"P_excl", "psi_rp"}, {1, 1} }, proton_weight, qvector_axes)
    .AddVectorCorrelation( "P_excl_sq", proton_weight, qvector_axes )

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"F1"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"F2"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"F3"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CovarianceDecorator<qvector_t>{ "F1", { 1 } }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CovarianceDecorator<qvector_t>{ "F2", { 1 } }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CovarianceDecorator<qvector_t>{ "F3", { 1 } }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"ini_F1"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"ini_F2"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddCorrelation( CorrelationDecorator<qvector_t>{ std::vector<std::string>{"ini_F3"}, {1} }, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
  ;

  auto histo1d = std::vector<ROOT::RDF::RResultPtr<::TH1D>>{};
  auto histo2d = std::vector<ROOT::RDF::RResultPtr<::TH2D>>{};

  histo2d.push_back( sampled_d.Histo2D( {"h2_modE_modId", ";Module ID; E_{loss}", 55, 0, 55, 300, 0, 30 }, "fhcalModId", "fhcalModE" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_qx_qy_f1", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "F1_x1", "F1_y1" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_qx_qy_f2", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "F2_x1", "F2_y1" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_qx_qy_f3", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "F3_x1", "F3_y1" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_ini_qx_qy_f1", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "ini_F1_x1", "ini_F1_y1" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_ini_qx_qy_f2", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "ini_F2_x1", "ini_F2_y1" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_ini_qx_qy_f3", ";Q_{x}; Q_{y}", 300, -1.5, 1.5, 300, -1.5, 1.5 }, "ini_F3_x1", "ini_F3_y1" ) );

  histo2d.push_back( sampled_d.Histo2D( {"h2_fhcal_x_y_f1", ";x (cm); y (cm)", 200, -100., 100., 200, -100., 100. }, "fhcalModX", "fhcalModY", "F1w" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_fhcal_x_y_f2", ";x (cm); y (cm)", 200, -100., 100., 200, -100., 100. }, "fhcalModX", "fhcalModY", "F2w" ) );
  histo2d.push_back( sampled_d.Histo2D( {"h2_fhcal_x_y_f3", ";x (cm); y (cm)", 200, -100., 100., 200, -100., 100. }, "fhcalModX", "fhcalModY", "F3w" ) );

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr_selfexcl.root", "RECREATE"), [](auto f){ f->Close(); } };
  handler.DumpCorrelations( file_out.get() );
  file_out->cd();
  file_out->mkdir("qa");
  file_out->cd("qa");
  std::for_each( histo1d.begin(), histo1d.end(), []( auto h1 ){ h1->Write(); } );
  std::for_each( histo2d.begin(), histo2d.end(), []( auto h1 ){ h1->Write(); } );

  auto n_events_filtered = static_cast<double>(*(sampled_d.Count())) / 1E+3;
  auto end = std::chrono::steady_clock::now();
  auto elapsed_s = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
  auto speed = n_events_filtered / elapsed_s * 3600;
  std::cout << "Elapsed time: " << elapsed_s << " sec" << std::endl;
  std::cout << "Processing speed: " << std::setprecision(3) << speed << " kev/h" << std::endl;
}
