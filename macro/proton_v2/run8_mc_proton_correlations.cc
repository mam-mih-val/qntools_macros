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

  auto proton_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "y", 6, 0.0, 1.2 },
    Qn::AxisD{ "pT", 5, 0.0, 2.0 },
  };
  
  auto qvector_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
  };

  constexpr size_t NHARM = 2;
  const auto l = double{5e-2};
  auto harmonics = std::vector<size_t>(NHARM);
  std::iota( harmonics.begin(), harmonics.end(), 1 );

  std::vector<int> f1_mod = {
    34, 35, 
    36, 37, 
    38, 39, 
    40, 41, 
    42, 43, 
  };
  std::vector<int> f2_mod = {
     0,  1,  2,  3,  4,
     5,  6,  7,  8,  9,
    10, 11, 12, 13, 14,
    15, 16,     17, 18, 
    19, 20, 21, 22, 23,
    24, 25, 26, 27, 28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_mod = {
    44, 45,
    46, 47,
    48, 49,
    50, 51, 
    52, 53,
  };

  std::for_each( f1_mod.begin(), f1_mod.end(), [](auto& m){ m += 1; } );
  std::for_each( f2_mod.begin(), f2_mod.end(), [](auto& m){ m += 1; } );
  std::for_each( f3_mod.begin(), f3_mod.end(), [](auto& m){ m += 1; } );

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("h2_efficiency_2212_good", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No efficiency for both tof was found in file " << str_effieciency_file << "\n";

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  
  auto dd = GenerateBmnExtendedTree(d, efficiency_histo);

  auto sampled_d = Qn::Correlation::Resample(dd, 100);

  sampled_d = sampled_d.Define( "F1w", fhcal_weight_generator(f1_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F2w", fhcal_weight_generator(f2_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F3w", fhcal_weight_generator(f3_mod), { "fhcalModId", "fhcalModE" } );

  DefineVector( sampled_d, "ini_F1", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F1w"} );
  DefineVector( sampled_d, "ini_F2", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F2w"} );
  DefineVector( sampled_d, "ini_F3", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F3w"} );
  DefineVector( sampled_d, "ini_Tpos", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTposW"} );
  DefineVector( sampled_d, "ini_Tneg", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTnegW"} );

  DefineVector(sampled_d, "ini_proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );
  DefineVector(sampled_d, "tru_proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"simPhi"s} );
  DefineVector(sampled_d, "psi_rp", psi_rp_vector< double >( harmonics ), std::vector<std::string>{"psiRP"s} );

  auto calib_file = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( str_calib_file.c_str(), "READ"), [](auto f){ f->Close(); } };
  auto [vec_p_mean, vec_p_cov] = ReadMeanCov<NHARM>("proton", calib_file.get());
  auto p_correction_container = MakeCorrectionContainer<NHARM>( vec_p_mean, vec_p_cov, PrincipalComponents<NHARM>{}, l );
  auto p_corr_builder = CorrectorBuilder<NHARM>( p_correction_container );
  sampled_d = sampled_d.Define( "proton", p_corr_builder.IssueUVectorCorrector<uvector_t, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >(), { "ini_proton", "centrality", "trProtonY", "trPt" } );

  auto [vec_f1_mean, vec_f1_cov] = ReadMeanCov<NHARM>("F1", calib_file.get());
  auto f1_correction_container = MakeCorrectionContainer<NHARM>( vec_f1_mean, vec_f1_cov, PrincipalComponents<NHARM>{}, l );
  auto f1_corr_builder = CorrectorBuilder<NHARM>( f1_correction_container );
  sampled_d = sampled_d.Define( "F1", f1_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F1", "centrality" } );

  auto [vec_f2_mean, vec_f2_cov] = ReadMeanCov<NHARM>("F2", calib_file.get());
  auto f2_correction_container = MakeCorrectionContainer<NHARM>( vec_f2_mean, vec_f2_cov, PrincipalComponents<NHARM>{}, l );
  auto f2_corr_builder = CorrectorBuilder<NHARM>( f2_correction_container );
  sampled_d = sampled_d.Define( "F2", f2_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F2", "centrality" } );

  auto [vec_f3_mean, vec_f3_cov] = ReadMeanCov<NHARM>("F3", calib_file.get());
  auto f3_correction_container = MakeCorrectionContainer<NHARM>( vec_f3_mean, vec_f3_cov, PrincipalComponents<NHARM>{}, l );
  auto f3_corr_builder = CorrectorBuilder<NHARM>( f3_correction_container );
  sampled_d = sampled_d.Define( "F3", f3_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_F3", "centrality" } );

  auto [vec_tp_mean, vec_tp_cov] = ReadMeanCov<NHARM>("Tpos", calib_file.get());
  auto tp_correction_container = MakeCorrectionContainer<NHARM>( vec_tp_mean, vec_tp_cov, PrincipalComponents<NHARM>{}, l );
  auto tp_corr_builder = CorrectorBuilder<NHARM>( tp_correction_container );
  sampled_d = sampled_d.Define( "Tpos", tp_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tpos", "centrality" } );

  auto [vec_tn_mean, vec_tn_cov] = ReadMeanCov<NHARM>("Tneg", calib_file.get());
  auto tn_correction_container = MakeCorrectionContainer<NHARM>( vec_tn_mean, vec_tn_cov, PrincipalComponents<NHARM>{}, l );
  auto tn_corr_builder = CorrectorBuilder<NHARM>( tn_correction_container );
  sampled_d = sampled_d.Define( "Tneg", tn_corr_builder.IssueQVectorCorrector<qvector_t, float>(), { "ini_Tneg", "centrality" } );
  
  auto p_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  // auto tru_p_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "tru_proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  // auto ini_p_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "ini_proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  auto p_f1_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "proton", "F1", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  auto p_f2_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "proton", "F2", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  auto p_f3_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "proton", "F3", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  
  auto f1_f2_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F1", "F2", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f2_f3_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F2", "F3", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f1_f3_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F1", "F3", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  
  auto f1_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F1", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f2_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F2", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f3_psi_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F3", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1} } );

  auto f1_tp_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F1", "Tpos", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f1_tn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F1", "Tneg", std::vector{ std::pair<size_t, size_t>{1, 1} } );

  auto f2_tp_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F2", "Tpos", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f2_tn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F2", "Tneg", std::vector{ std::pair<size_t, size_t>{1, 1} } );

  auto f3_tp_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F3", "Tpos", std::vector{ std::pair<size_t, size_t>{1, 1} } );
  auto f3_tn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< qvector_t, qvector_t >{}, "F3", "Tneg", std::vector{ std::pair<size_t, size_t>{1, 1} } );

  auto p_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  // auto ini_p_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  // auto tru_p_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto p_f1_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto p_f2_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto p_f3_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  
  auto f1_f2_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f1_f3_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f2_f3_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  
  auto f1_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f2_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f3_psi_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};

  auto f1_tp_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f1_tn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};

  auto f2_tp_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f2_tn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};

  auto f3_tp_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto f3_tn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};

  sampled_d = sampled_d.Define( "is_corrected", [&p_correction_container]( std::vector<double> weights,  float centrality, ROOT::VecOps::RVec<float> vec_y, ROOT::VecOps::RVec<float> vec_pT ){
    auto vec_corrected = std::vector<double>( weights.size(), 0 );
    for( auto i=size_t{}; i<weights.size(); ++i ){
      auto pT = vec_pT[i];
      auto y = vec_y[i];
      auto coord = std::vector<double>{ centrality, y, pT };
      auto bin = p_correction_container.FindBin( coord );
      if( bin > p_correction_container.size() )
        continue;
      if( bin < 0 )
        continue;
      auto [Mpinv, c] = p_correction_container[bin];
      if( std::isnan(Mpinv(0,0)) )
        continue;
      vec_corrected[i] = weights[i];
    }
    return vec_corrected;
  }, {"trProtonWeight", "centrality", "trProtonY", "trPt"} );

  for( const auto& name : p_psi_names ){
    p_psi_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "is_corrected", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  // for( const auto& name : ini_p_psi_names ){
  //   ini_p_psi_ptr.emplace_back(
  //     sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "is_corrected", "samples", "centrality", "trProtonY", "trPt" } )
  //   ); 
  // }

  // for( const auto& name : tru_p_psi_names ){
  //   tru_p_psi_ptr.emplace_back(
  //     sampled_d.Book< std::vector<double>, std::vector<int>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, std::vector<float> >( CorrelationHelper(proton_axes), {name, "simIsProton", "samples", "centrality", "simProtonY", "simPt" } )
  //   ); 
  // }

  for( const auto& name : p_f1_names ){
    p_f1_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "is_corrected", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : p_f2_names ){
    p_f2_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "is_corrected", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : p_f3_names ){
    p_f3_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "is_corrected", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : f1_f2_names ){
    f1_f2_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f1_f3_names ){
    f1_f3_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_f3_names ){
    f2_f3_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f1_psi_names ){
    f1_psi_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_psi_names ){
    f2_psi_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_psi_names ){
    f3_psi_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f1_tp_names ){
    f1_tp_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f1_tn_names ){
    f1_tn_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_tp_names ){
    f2_tp_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_tn_names ){
    f2_tn_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_tp_names ){
    f3_tp_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_tn_names ){
    f3_tn_ptr.emplace_back(
      sampled_d.Book< double, double,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes, 100, Qn::Stat::WeightType::REFERENCE), {name, "One", "samples", "centrality" } )
    ); 
  }


  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  
  std::for_each( p_psi_ptr.begin(), p_psi_ptr.end(), [i=0, &p_psi_names]( auto& p ) mutable { p->Write( p_psi_names.at(i).c_str() ); ++i; } );
  // std::for_each( ini_p_psi_ptr.begin(), ini_p_psi_ptr.end(), [i=0, &ini_p_psi_names]( auto& p ) mutable { p->Write( ini_p_psi_names.at(i).c_str() ); ++i; } );
  // std::for_each( tru_p_psi_ptr.begin(), tru_p_psi_ptr.end(), [i=0, &tru_p_psi_names]( auto& p ) mutable { p->Write( tru_p_psi_names.at(i).c_str() ); ++i; } );
  std::for_each( p_f1_ptr.begin(), p_f1_ptr.end(), [i=0, &p_f1_names]( auto& p ) mutable { p->Write( p_f1_names.at(i).c_str() ); ++i; } );
  std::for_each( p_f2_ptr.begin(), p_f2_ptr.end(), [i=0, &p_f2_names]( auto& p ) mutable { p->Write( p_f2_names.at(i).c_str() ); ++i; } );
  std::for_each( p_f3_ptr.begin(), p_f3_ptr.end(), [i=0, &p_f3_names]( auto& p ) mutable { p->Write( p_f3_names.at(i).c_str() ); ++i; } );
  
  std::for_each( f1_f2_ptr.begin(), f1_f2_ptr.end(), [i=0, &f1_f2_names]( auto& p ) mutable { p->Write( f1_f2_names.at(i).c_str() ); ++i; } );
  std::for_each( f1_f3_ptr.begin(), f1_f3_ptr.end(), [i=0, &f1_f3_names]( auto& p ) mutable { p->Write( f1_f3_names.at(i).c_str() ); ++i; } );
  std::for_each( f2_f3_ptr.begin(), f2_f3_ptr.end(), [i=0, &f2_f3_names]( auto& p ) mutable { p->Write( f2_f3_names.at(i).c_str() ); ++i; } );

  std::for_each( f1_psi_ptr.begin(), f1_psi_ptr.end(), [i=0, &f1_psi_names]( auto& p ) mutable { p->Write( f1_psi_names.at(i).c_str() ); ++i; } );
  std::for_each( f2_psi_ptr.begin(), f2_psi_ptr.end(), [i=0, &f2_psi_names]( auto& p ) mutable { p->Write( f2_psi_names.at(i).c_str() ); ++i; } );
  std::for_each( f3_psi_ptr.begin(), f3_psi_ptr.end(), [i=0, &f3_psi_names]( auto& p ) mutable { p->Write( f3_psi_names.at(i).c_str() ); ++i; } );
  
  std::for_each( f1_tp_ptr.begin(), f1_tp_ptr.end(), [i=0, &f1_tp_names]( auto& p ) mutable { p->Write( f1_tp_names.at(i).c_str() ); ++i; } );
  std::for_each( f1_tn_ptr.begin(), f1_tn_ptr.end(), [i=0, &f1_tn_names]( auto& p ) mutable { p->Write( f1_tn_names.at(i).c_str() ); ++i; } );

  std::for_each( f2_tp_ptr.begin(), f2_tp_ptr.end(), [i=0, &f2_tp_names]( auto& p ) mutable { p->Write( f2_tp_names.at(i).c_str() ); ++i; } );
  std::for_each( f2_tn_ptr.begin(), f2_tn_ptr.end(), [i=0, &f2_tn_names]( auto& p ) mutable { p->Write( f2_tn_names.at(i).c_str() ); ++i; } );

  std::for_each( f3_tp_ptr.begin(), f3_tp_ptr.end(), [i=0, &f3_tp_names]( auto& p ) mutable { p->Write( f3_tp_names.at(i).c_str() ); ++i; } );
  std::for_each( f3_tn_ptr.begin(), f3_tn_ptr.end(), [i=0, &f3_tn_names]( auto& p ) mutable { p->Write( f3_tn_names.at(i).c_str() ); ++i; } );


  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}