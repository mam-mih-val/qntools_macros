//
// Created by Misha on 3/7/2023.
//

#include <algorithm>
#include <cassert>
#include <cmath>
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

void run8_mc_proton_fill( std::string list, std::string str_effieciency_file ){

  std::cout << "starting execution" << std::endl;

  auto proton_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "y", 12, 0.0, 1.2 },
    Qn::AxisD{ "pT", 10, 0.0, 2.0 },
  };

  auto qvector_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
  };

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

  auto harmonics = std::vector<size_t>( 5 );
  std::iota( harmonics.begin(), harmonics.end(), 1 );

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
  ;

  auto sampled_d = Qn::Correlation::Resample(dd, 100);

  sampled_d = sampled_d.Define( "F1w", fhcal_weight_generator(f1_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F2w", fhcal_weight_generator(f2_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F3w", fhcal_weight_generator(f3_mod), { "fhcalModId", "fhcalModE" } );

  DefineQnToolsLikeWeights(sampled_d, "proton_qntools_like", QnToolsLikeWeightsFunc<std::vector<double>, ROOT::VecOps::RVec<ULong64_t>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> > {
    std::vector<Qn::AxisD>{
      Qn::AxisD{ "y", 12, 0.0, 1.2 },
      Qn::AxisD{ "pT", 10, 0.0, 2.0 },
    }
  } , std::vector<std::string>{
    "trProtonWeight", "samples", "trProtonY", "trPt"
  });

  // DefineVector( sampled_d, "F1", q_vector< ROOT::VecOps::RVec<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F1w"} );
  // DefineVector( sampled_d, "F2", q_vector< ROOT::VecOps::RVec<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F2w"} );
  // DefineVector( sampled_d, "F3", q_vector< ROOT::VecOps::RVec<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F3w"} );
  
  // DefineVector( sampled_d, "Tpos", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTposW"} );
  // DefineVector( sampled_d, "Tneg", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTnegW"} );

  auto p_components_names = AddUVectorComponents(sampled_d, "proton", harmonics, "trPhi" );
  auto p_cov_names = AddUVectorCovariance(sampled_d, "proton", harmonics, "trPhi" );

  auto f1_means_str = AddUVectorComponents( sampled_d, "F1", harmonics, "fhcalModPhi" );
  auto f2_means_str = AddUVectorComponents( sampled_d, "F2", harmonics, "fhcalModPhi" );
  auto f3_means_str = AddUVectorComponents( sampled_d, "F3", harmonics, "fhcalModPhi" );
  auto tp_means_str = AddUVectorComponents( sampled_d, "Tpos", harmonics, "trPhi" );
  auto tn_means_str = AddUVectorComponents( sampled_d, "Tneg", harmonics, "trPhi" );

  auto f1_cov_str = AddUVectorCovariance( sampled_d, "F1", harmonics, "fhcalModPhi" );
  auto f2_cov_str = AddUVectorCovariance( sampled_d, "F2", harmonics, "fhcalModPhi" );
  auto f3_cov_str = AddUVectorCovariance( sampled_d, "F3", harmonics, "fhcalModPhi" );
  auto tp_cov_str = AddUVectorCovariance( sampled_d, "Tpos", harmonics, "trPhi" );
  auto tn_cov_str = AddUVectorCovariance( sampled_d, "Tneg", harmonics, "trPhi" );

  auto p_components_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  p_components_ptr.reserve( p_components_names.size() );
  auto p_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  p_cov_ptr.reserve( p_cov_names.size() );

  auto f1_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f1_means_ptr.reserve( f1_means_str.size() );
  auto f1_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f1_cov_ptr.reserve( f1_cov_str.size() );

  auto f2_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f2_means_ptr.reserve( f2_means_str.size() );
  auto f2_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f2_cov_ptr.reserve( f2_cov_str.size() );

  auto f3_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f3_means_ptr.reserve( f3_means_str.size() );
  auto f3_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f3_cov_ptr.reserve( f3_cov_str.size() );

  auto tp_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tp_means_ptr.reserve( tp_means_str.size() );
  auto tp_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tp_cov_ptr.reserve( tp_cov_str.size() );

  auto tn_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tn_means_ptr.reserve( tn_means_str.size() );
  auto tn_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tn_cov_ptr.reserve( tn_cov_str.size() );

  for( const auto& name : p_components_names ){
    p_components_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), std::vector<std::string>{name, "proton_qntools_like", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : p_cov_names ){
    p_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), std::vector<std::string>{name, "proton_qntools_like", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : f1_means_str ){
    f1_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F1w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_means_str ){
    f2_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_means_str ){
    f3_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F3w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tp_means_str ){
    tp_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTposW", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tn_means_str ){
    tn_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTnegW", "samples", "centrality" } )
    ); 
  }


  for( const auto& name : f1_cov_str ){
    f1_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F1w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_cov_str ){
    f2_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_cov_str ){
    f3_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tp_cov_str ){
    tp_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTposW", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tn_cov_str ){
    tn_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTnegW", "samples", "centrality" } )
    ); 
  }

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  std::for_each( p_components_ptr.begin(), p_components_ptr.end(), [i=0, &p_components_names]( auto& p ) mutable { p->Write( p_components_names.at(i).c_str() ); ++i; } );
  std::for_each( p_cov_ptr.begin(), p_cov_ptr.end(), [i=0, &p_cov_names]( auto& p ) mutable { p->Write( p_cov_names.at(i).c_str() ); ++i; } );
  
  std::for_each( f1_means_ptr.begin(), f1_means_ptr.end(), [i=0, &f1_means_str]( auto& p ) mutable { p->Write( f1_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f1_cov_ptr.begin(), f1_cov_ptr.end(), [i=0, &f1_cov_str]( auto& p ) mutable { p->Write( f1_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( f2_means_ptr.begin(), f2_means_ptr.end(), [i=0, &f2_means_str]( auto& p ) mutable { p->Write( f2_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f2_cov_ptr.begin(), f2_cov_ptr.end(), [i=0, &f2_cov_str]( auto& p ) mutable { p->Write( f2_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( f3_means_ptr.begin(), f3_means_ptr.end(), [i=0, &f3_means_str]( auto& p ) mutable { p->Write( f3_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f3_cov_ptr.begin(), f3_cov_ptr.end(), [i=0, &f3_cov_str]( auto& p ) mutable { p->Write( f3_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( tp_means_ptr.begin(), tp_means_ptr.end(), [i=0, &tp_means_str]( auto& p ) mutable { p->Write( tp_means_str.at(i).c_str() ); ++i; } );
  std::for_each( tp_cov_ptr.begin(), tp_cov_ptr.end(), [i=0, &tp_cov_str]( auto& p ) mutable { p->Write( tp_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( tn_means_ptr.begin(), tn_means_ptr.end(), [i=0, &tn_means_str]( auto& p ) mutable { p->Write( tn_means_str.at(i).c_str() ); ++i; } );
  std::for_each( tn_cov_ptr.begin(), tn_cov_ptr.end(), [i=0, &tn_cov_str]( auto& p ) mutable { p->Write( tn_cov_str.at(i).c_str() ); ++i; } );

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}