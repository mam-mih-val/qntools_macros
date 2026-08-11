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
#include "corrections.h"

void run8_mc_proton_correlations( std::string list, std::string str_effieciency_file, std::string str_calib_file ){

  std::cout << "starting execution" << std::endl;

  auto proton_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "y", 6, 0.0, 1.2 },
    Qn::AxisD{ "pT", 5, 0.0, 2.0 },
  };
  constexpr size_t NHARM = 10;
  auto harmonics = std::vector<size_t>(NHARM);
  std::iota( harmonics.begin(), harmonics.end(), 1 );

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH3* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("h3_efficiency_2212_good", efficiency_histo);
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

  DefineVector(sampled_d, "ini_proton", "trPhi", u_vector< std::vector<float> >( harmonics ) );
  DefineVector(sampled_d, "tru_proton", "simPhi", u_vector< std::vector<float> >( harmonics ) );
  DefineVector(sampled_d, "psi_rp", "psiRP", psi_rp_vector< double >( harmonics ) );

  auto calib_file = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( str_calib_file.c_str(), "READ"), [](auto f){ f->Close(); } };
  auto [vec_p_mean, vec_p_cov] = ReadMeanCov<NHARM>("proton", calib_file.get());
  auto correction_container = MakeCorrectionContainer<NHARM>( vec_p_mean, vec_p_cov, MakeWhiteningMatrixFunc<NHARM>() );
  auto corr_builder = CorrectorBuilder<NHARM>( correction_container );

  sampled_d = sampled_d.Define( "proton", corr_builder.IssueUVectorCorrector<uvector_t, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >(), { "ini_proton", "centrality", "trProtonY", "trPt" } );
  
  auto vn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  auto ini_vn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "ini_proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2} } );
  // auto tru_vn_names = Define2PartCorrelation( sampled_d, CorrFunc2Part< uvector_t, qvector_t >{}, "tru_proton", "psi_rp", std::vector{ std::pair<size_t, size_t>{1, 1}, std::pair<size_t, size_t>{2, 2}, std::pair<size_t, size_t>{3, 3} } );

  auto vn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  vn_ptr.reserve( vn_names.size() );

  auto ini_vn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  ini_vn_ptr.reserve( ini_vn_names.size() );

  // auto tru_vn_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  // tru_vn_ptr.reserve( tru_vn_names.size() );

  for( const auto& name : vn_names ){
    vn_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : ini_vn_names ){
    ini_vn_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  // for( const auto& name : tru_vn_names ){
  //   tru_vn_ptr.emplace_back(
  //     sampled_d.Book< std::vector<double>, std::vector<int>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, std::vector<float> >( CorrelationHelper(proton_axes), {name, "simIsProton", "samples", "centrality", "simProtonY", "simPt" } )
  //   ); 
  // }

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  std::for_each( vn_ptr.begin(), vn_ptr.end(), [i=0, &vn_names]( auto& p ) mutable { p->Write( vn_names.at(i).c_str() ); ++i; } );
  std::for_each( ini_vn_ptr.begin(), ini_vn_ptr.end(), [i=0, &ini_vn_names]( auto& p ) mutable { p->Write( ini_vn_names.at(i).c_str() ); ++i; } );
  // std::for_each( tru_vn_ptr.begin(), tru_vn_ptr.end(), [i=0, &tru_vn_names]( auto& p ) mutable { p->Write( tru_vn_names.at(i).c_str() ); ++i; } );

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}