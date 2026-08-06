//
// Created by Misha on 3/7/2023.
//

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

const auto u_generator( size_t harmonic, std::function< double(double) > component ){
  return [harmonic, component]( std::vector<float> vec_phi ){
    auto vec_results = std::vector<double>{};
    vec_results.reserve(vec_phi.size());
    for( auto phi : vec_phi ){
      vec_results.push_back( component( harmonic*phi ) );
    }
    return vec_results;
  };
}

const auto cov_generator( size_t h_a, std::function< double(double) > c_a, size_t h_b, std::function< double(double) > c_b ){
  return [h_a, h_b, c_a, c_b]( std::vector<float> vec_phi ){
    auto vec_results = std::vector<double>{};
    vec_results.reserve(vec_phi.size());
    for( auto phi : vec_phi ){
      vec_results.push_back( c_a( h_a*phi ) * c_b( h_b*phi )  );
    }
    return vec_results;
  };
}

void run8_mc_proton_fill( std::string list, std::string str_effieciency_file ){

  std::cout << "starting execution" << std::endl;

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
  auto x1_corr = sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "y", 6, 0.0, 1.2 },
    Qn::AxisD{ "pT", 5, 0.0, 2.0 },

  }), {"x1", "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } );

  
  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  x1_corr->Write( "proton.x1" );

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}