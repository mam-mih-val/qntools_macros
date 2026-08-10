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
    Qn::AxisD{ "y", 6, 0.0, 1.2 },
    Qn::AxisD{ "pT", 5, 0.0, 2.0 },
  };
  auto harmonics = std::vector<size_t>(4);
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
  ;

  auto sampled_d = Qn::Correlation::Resample(dd, 100);

  auto p_components_names = AddUVectorComponents(sampled_d, "proton", harmonics, "trPhi" );
  auto p_cov_names = AddUVectorCovariance(sampled_d, "proton", harmonics, "trPhi" );

  auto p_components_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  auto p_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  p_components_ptr.reserve( p_components_names.size() );
  p_cov_ptr.reserve( p_cov_names.size() );

  for( const auto& name : p_components_names ){
    p_components_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : p_cov_names ){
    p_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), {name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  std::for_each( p_components_ptr.begin(), p_components_ptr.end(), [i=0, &p_components_names]( auto& p ) mutable { p->Write( p_components_names.at(i).c_str() ); ++i; } );
  std::for_each( p_cov_ptr.begin(), p_cov_ptr.end(), [i=0, &p_cov_names]( auto& p ) mutable { p->Write( p_cov_names.at(i).c_str() ); ++i; } );

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}