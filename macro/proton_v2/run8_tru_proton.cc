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

  constexpr size_t NHARM = 4;
  constexpr size_t Q_NHARM = 1;

  const auto l = double{5e-2};
  auto harmonics = std::vector<size_t>(NHARM);
  auto harmonics_q = std::vector<size_t>(Q_NHARM);
  std::iota( harmonics.begin(), harmonics.end(), 1 );
  std::iota( harmonics_q.begin(), harmonics_q.end(), 1 );

  
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

  DefineVector(sampled_d, "tru_proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"simPhi"s} );
  DefineVector(sampled_d, "simF1", q_vector< std::vector<float>, std::vector<double> >( harmonics_q ), std::vector<std::string>{ "simPhi", "simF1w" } );
  DefineVector(sampled_d, "simF2", q_vector< std::vector<float>, std::vector<double> >( harmonics_q ), std::vector<std::string>{ "simPhi", "simF2w" } );
  DefineVector(sampled_d, "simF3", q_vector< std::vector<float>, std::vector<double> >( harmonics_q ), std::vector<std::string>{ "simPhi", "simF3w" } );  
  DefineVector(sampled_d, "psi_rp", psi_rp_vector< double >( harmonics ), std::vector<std::string>{"psiRP"s} );

  auto tru_proton_axes = CorrelationAxes<float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>>{
    .axes_columns={ "centrality", "simProtonY", "simPt" },
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

  auto tru_proton_weight = Weight<std::vector<double>>{ "simIsPrimProton" };
  auto qvector_weight = Weight<double>{ "One" };

  auto begin = std::chrono::steady_clock::now();
  
  auto handler = CorrelationHandler{ sampled_d, 100 };
  handler
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"tru_proton", "psi_rp"}, {1, 1} }, tru_proton_weight, tru_proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"tru_proton", "psi_rp"}, {2, 2} }, tru_proton_weight, tru_proton_axes)

    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"tru_proton", "simF1"}, {1, 1} }, tru_proton_weight, tru_proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"tru_proton", "simF2"}, {1, 1} }, tru_proton_weight, tru_proton_axes)
    .AddCorrelation( CorrelationDecorator<uvector_t, qvector_t>{ std::vector<std::string>{"tru_proton", "simF3"}, {1, 1} }, tru_proton_weight, tru_proton_axes)
    
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"simF1", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"simF2", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes)
    .AddCorrelation( CorrelationDecorator<qvector_t, qvector_t>{ std::vector<std::string>{"simF3", "psi_rp"}, {1, 1} }, qvector_weight, qvector_axes)
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