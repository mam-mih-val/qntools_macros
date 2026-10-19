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
    Qn::AxisD{ "y", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2 } },
    // Qn::AxisD{ "y", std::vector<double>{ 0, 0.4, 0.6, 0.8, 1, 1.2 } },
    Qn::AxisD{ "pT", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5,  2.0 } },
    // Qn::AxisD{ "pT", std::vector<double>{ 0, 0.4, 0.6, 1.0, 1.5, 2.0 } },
  };

  auto qvector_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
  };

  // std::vector<int> f1_mod = {
  //   6,  7,  8,
  //   11, 12, 13,
  //   16,     17,
  //   20, 21, 22, 
  //   25, 26, 27
  // };
  // std::vector<int> f2_mod = {
  //   0,  1,  2,  3,  4,
  //   5,              9,
  //   10,             14,
  //   15,             18,
  //   19,             23,
  //   24,             28,
  //   29, 30, 31, 32, 33,
  // };
  // std::vector<int> f3_mod = {
  //   35,                 44,
  //   37,                 46, 
  //   39,                 48, 
  //   41,                 50,
  //   43,                 52
  // };

  // std::vector<int> f4_mod = {
  //   34,                     45,
  //   36,                     47, 
  //   38,                     49, 
  //   40,                     51,
  //   42,                     53
  // };

  // std::for_each( f1_mod.begin(), f1_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f2_mod.begin(), f2_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f3_mod.begin(), f3_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f4_mod.begin(), f4_mod.end(), [](auto& m){ m += 1; } );

  std::vector<int> f1_mod = {
    12, 13, 14,
    17,         18,
    21, 22, 23, 
  };

  std::vector<int> f2_mod = {
    1, 2,
    6, 7 
  };

  std::vector<int> f3_mod = {
    25, 26,
    30, 31 
  };

  std::vector<int> f4_mod = {
    35, 36,
    37, 38 
  };

  std::vector<int> f5_mod = {
    41, 42,
    43, 44 
  };

  auto harmonics = std::vector<size_t>( 4 );
  std::iota( harmonics.begin(), harmonics.end(), 1 );

  auto qvector_harmonics = std::vector<size_t>( 2 );
  std::iota( qvector_harmonics.begin(), qvector_harmonics.end(), 1 );

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("h2_efficiency_2212_good", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No was found in file " << str_effieciency_file << "\n";

  TStopwatch timer;
  timer.Start();
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
  sampled_d = sampled_d.Define( "F5w", fhcal_weight_generator(f5_mod), { "fhcalModId", "fhcalModE" } );

  DefineVector( sampled_d, "F1", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F1w"} );
  DefineVector( sampled_d, "F2", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F2w"} );
  DefineVector( sampled_d, "F3", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F3w"} );
  DefineVector( sampled_d, "F4", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F4w"} );
  DefineVector( sampled_d, "F5", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"fhcalModPhi", "F5w"} );
  
  DefineVector( sampled_d, "T1", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trT1W"} );
  DefineVector( sampled_d, "T2", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trT2W"} );
  DefineVector( sampled_d, "T3", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trT3W"} );

  DefineVector( sampled_d, "Tpos", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTposW"} );
  DefineVector( sampled_d, "Tneg", q_vector< std::vector<float>, std::vector<double> >(harmonics), std::vector<std::string>{"trPhi", "trTnegW"} );

  auto proton_axes = CorrelationAxes<float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>>{
    .axes_columns={ "centrality", "trProtonY", "trPt" },
    .axes=std::vector<Qn::AxisD>{
      Qn::AxisD{ "centrality", 6, 0, 60 },
      Qn::AxisD{ "y", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 } },
      // Qn::AxisD{ "y", std::vector<double>{ 0, 0.4, 0.6, 0.8, 1, 1.2 } },
      Qn::AxisD{ "pT", std::vector<double>{ 0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5 } },
      // Qn::AxisD{ "pT", std::vector<double>{ 0, 0.4, 0.6, 1.0, 1.5, 2.0 } },
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

  auto start = std::chrono::steady_clock::now();
  auto handler = CorrelationHandler{ sampled_d, 100 };
  handler
    .AddComponentMeans( Vector<uvector_t>{ "proton" }, harmonics, proton_weight, proton_axes)
    .AddComponentCovariance( Vector<uvector_t>{ "proton" }, harmonics, proton_weight, proton_axes)

    .AddComponentMeans( Vector<qvector_t>{ "F1" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "F1" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "F2" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "F2" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "F3" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "F3" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "F4" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "F4" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "F5" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "F5" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "T1" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "T1" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "T2" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "T2" }, qvector_axes, qvector_weight, qvector_axes)

    .AddComponentMeans( Vector<qvector_t>{ "T3" }, qvector_axes, qvector_weight, qvector_axes)
    .AddComponentCovariance( Vector<qvector_t>{ "T3" }, qvector_axes, qvector_weight, qvector_axes)
    
  ;

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  handler.DumpCorrelations( file_out.get() );
  
  auto n_events_filtered = static_cast<double>(*(sampled_d.Count())) / 1E+3;
  auto finish = std::chrono::steady_clock::now();
  auto elapsed_s = std::chrono::duration_cast<std::chrono::seconds>(finish - start).count();
  auto speed = n_events_filtered / elapsed_s * 3600;
  std::cout << "Elapsed time: " << elapsed_s << " sec" << std::endl;
  std::cout << "Processing speed: " << std::setprecision(3) << speed << " kev/h" << std::endl;
}