//
// Created by Misha on 3/7/2023.
//

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

void run8_mc_proton_correlations( std::string list, std::string str_effieciency_file, std::string str_calib_file ){

  std::cout << "starting execution" << std::endl;

  constexpr size_t NHARM = 8;
  constexpr size_t Q_NHARM = 1;

  const auto l = double{5e-2};
  auto harmonics = std::vector<size_t>(NHARM);
  auto harmonics_q = std::vector<size_t>(Q_NHARM);
  std::iota( harmonics.begin(), harmonics.end(), 1 );
  std::iota( harmonics_q.begin(), harmonics_q.end(), 1 );

  std::vector<int> f1_mod = {
    6,  7,  8,
    11, 12, 13,
    16,     17,
    20, 21, 22, 
    25, 26, 27
  };
  std::vector<int> f2_mod = {
    0,  1,  2,  3,  4,
    5,              9,
    10,             14,
    15,             18,
    19,             23,
    24,             28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_mod = {
    35,                 44,
    37,                 46, 
    39,                 48, 
    41,                 50,
    43,                 52
  };

  std::vector<int> f4_mod = {
    34,                     45,
    36,                     47, 
    38,                     49, 
    40,                     51,
    42,                     53
  };

  std::for_each( f1_mod.begin(), f1_mod.end(), [](auto& m){ m += 1; } );
  std::for_each( f2_mod.begin(), f2_mod.end(), [](auto& m){ m += 1; } );
  std::for_each( f3_mod.begin(), f3_mod.end(), [](auto& m){ m += 1; } );
  std::for_each( f4_mod.begin(), f4_mod.end(), [](auto& m){ m += 1; } );

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
  sampled_d = sampled_d.Define( "F4w", fhcal_weight_generator(f4_mod), { "fhcalModId", "fhcalModE" } );

  DefineVector( sampled_d, "ini_F1", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"fhcalModPhi", "F1w"} );
  DefineVector( sampled_d, "ini_F2", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"fhcalModPhi", "F2w"} );
  DefineVector( sampled_d, "ini_F3", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"fhcalModPhi", "F3w"} );
  DefineVector( sampled_d, "ini_F4", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"fhcalModPhi", "F4w"} );
  DefineVector( sampled_d, "ini_Tpos", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"trPhi", "trTposW"} );
  DefineVector( sampled_d, "ini_Tneg", q_vector< std::vector<float>, std::vector<double> >(harmonics_q), std::vector<std::string>{"trPhi", "trTnegW"} );

  DefineVector(sampled_d, "ini_proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );
  DefineVector(sampled_d, "psi_rp", psi_rp_vector< double >( harmonics ), std::vector<std::string>{"psiRP"s} );

  auto calib_file = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( str_calib_file.c_str(), "READ"), [](auto f){ f->Close(); } };
  auto [vec_p_mean, vec_p_cov] = ReadMeanCov<NHARM>("proton", calib_file.get());
  auto p_correction_container = MakeCorrectionContainer<NHARM>( vec_p_mean, vec_p_cov, PrincipalComponents<NHARM>{}, l );
  auto p_corr_builder = CorrectorBuilder<NHARM>( p_correction_container );
  sampled_d = sampled_d.Define( "proton", p_corr_builder.IssueUVectorCorrector<uvector_t, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >(), { "ini_proton", "centrality", "trProtonY", "trPt" } );

  auto [vec_f1_mean, vec_f1_cov] = ReadMeanCov<2*Q_NHARM>("F1", calib_file.get());
  auto f1_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f1_mean, vec_f1_cov, TwistRescale<Q_NHARM>{}, l );
  auto f1_corr_builder = CorrectorBuilder<Q_NHARM>( f1_correction_container );
  sampled_d = sampled_d.Define( "F1", f1_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F1", "centrality" } );

  auto [vec_f2_mean, vec_f2_cov] = ReadMeanCov<2*Q_NHARM>("F2", calib_file.get());
  auto f2_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f2_mean, vec_f2_cov, TwistRescale<Q_NHARM>{}, l );
  auto f2_corr_builder = CorrectorBuilder<Q_NHARM>( f2_correction_container );
  sampled_d = sampled_d.Define( "F2", f2_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F2", "centrality" } );

  auto [vec_f3_mean, vec_f3_cov] = ReadMeanCov<2*Q_NHARM>("F3", calib_file.get());
  auto f3_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f3_mean, vec_f3_cov, TwistRescale<Q_NHARM>{}, l );
  auto f3_corr_builder = CorrectorBuilder<Q_NHARM>( f3_correction_container );
  sampled_d = sampled_d.Define( "F3", f3_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F3", "centrality" } );

  auto [vec_f4_mean, vec_f4_cov] = ReadMeanCov<2*Q_NHARM>("F4", calib_file.get());
  auto f4_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_f4_mean, vec_f4_cov, TwistRescale<Q_NHARM>{}, l );
  auto f4_corr_builder = CorrectorBuilder<Q_NHARM>( f4_correction_container );
  sampled_d = sampled_d.Define( "F4", f4_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F4", "centrality" } );

  auto [vec_tp_mean, vec_tp_cov] = ReadMeanCov<2*Q_NHARM>("Tpos", calib_file.get());
  auto tp_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tp_mean, vec_tp_cov, TwistRescale<Q_NHARM>{}, l );
  auto tp_corr_builder = CorrectorBuilder<Q_NHARM>( tp_correction_container );
  sampled_d = sampled_d.Define( "Tpos", tp_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tpos", "centrality" } );

  auto [vec_tn_mean, vec_tn_cov] = ReadMeanCov<2*Q_NHARM>("Tneg", calib_file.get());
  auto tn_correction_container = MakeCorrectionContainer<Q_NHARM>( vec_tn_mean, vec_tn_cov, TwistRescale<Q_NHARM>{}, l );
  auto tn_corr_builder = CorrectorBuilder<Q_NHARM>( tn_correction_container );
  sampled_d = sampled_d.Define( "Tneg", tn_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tneg", "centrality" } );

  auto proton_axes = CorrelationAxes<float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>>{
    .axes_columns={ "centrality", "trProtonY", "trPt" },
    .axes=std::vector<Qn::AxisD>{
      Qn::AxisD{ "centrality", 6, 0, 60 },
      Qn::AxisD{ "y", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2 } },
      Qn::AxisD{ "pT", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0 } },
    }
  };

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
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"proton", "psi_rp"}, {1, 1} }, proton_weight, proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"proton", "psi_rp"}, {2, 2} }, proton_weight, proton_axes)
    
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t, qvector_t>{ std::vector<std::string>{"proton", "F1", "F3"}, {2, 1, 1} }, proton_weight, proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t, qvector_t>{ std::vector<std::string>{"proton", "F1", "F4"}, {2, 1, 1} }, proton_weight, proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t, qvector_t>{ std::vector<std::string>{"proton", "F2", "F4"}, {2, 1, 1} }, proton_weight, proton_axes)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "F2"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "F3"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "F4"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "F3"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "F4"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "F4"}, {1, 1} }, qvector_weight, qvector_axes)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F4", "Tpos"}, {1, 1} }, qvector_weight, qvector_axes)

    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F1", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F2", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F3", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"F4", "Tneg"}, {1, 1} }, qvector_weight, qvector_axes)
  ;

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  handler.DumpCorrelations( file_out.get() );
  
  auto n_events_filtered = static_cast<double>(*(sampled_d.Count())) / 1E+3;
  auto end = std::chrono::steady_clock::now();
  auto elapsed_s = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
  auto speed = n_events_filtered / elapsed_s * 3600;
  std::cout << "Elapsed time: " << elapsed_s << " sec" << std::endl;
  std::cout << "Processing speed: " << std::setprecision(3) << speed << " kev/h" << std::endl;
}