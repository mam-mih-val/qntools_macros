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

  auto f1_mod = std::vector<int>(10);
  std::iota( f1_mod.begin(), f1_mod.end(), 45 );
  auto f2_mod = std::vector<int>(34);
  std::iota( f2_mod.begin(), f2_mod.end(), 1 );
  auto f3_mod = std::vector<int>(10);
  std::iota( f3_mod.begin(), f3_mod.end(), 35 );

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

  DefineVector( sampled_d, "F1", q_vector< std::vector<float>, std::vector<double> >(qvector_harmonics), std::vector<std::string>{"fhcalModPhi", "F1w"} );
  DefineVector( sampled_d, "F2", q_vector< std::vector<float>, std::vector<double> >(qvector_harmonics), std::vector<std::string>{"fhcalModPhi", "F2w"} );
  DefineVector( sampled_d, "F3", q_vector< std::vector<float>, std::vector<double> >(qvector_harmonics), std::vector<std::string>{"fhcalModPhi", "F3w"} );
  
  // DefineVector( sampled_d, "Tpos", q_vector< std::vector<float>, std::vector<double> >(qvector_harmonics), std::vector<std::string>{"trPhi", "trTposW"} );
  // DefineVector( sampled_d, "Tneg", q_vector< std::vector<float>, std::vector<double> >(qvector_harmonics), std::vector<std::string>{"trPhi", "trTnegW"} );

  DefineVector(sampled_d, "proton", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );
  DefineVector(sampled_d, "Tpos", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );
  DefineVector(sampled_d, "Tneg", u_vector< std::vector<float> >( harmonics ), std::vector<std::string>{"trPhi"s} );

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

  auto t_axes = CorrelationAxes<float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>>{
    .axes_columns={ "centrality", "trEta", "trPt" },
    .axes=std::vector<Qn::AxisD>{
      Qn::AxisD{ "centrality", 6, 0, 60 },
      Qn::AxisD{ "eta", 20, 2, 4 },
      Qn::AxisD{ "trPt", 5, 0, 2.0 },
    }
  };

  auto qvector_axes = CorrelationAxes<float>{
    .axes_columns={ "centrality" },
    .axes=std::vector<Qn::AxisD>{
      Qn::AxisD{ "centrality", 6, 0, 60 },
    }
  };

  auto proton_weight = Weight<std::vector<double>>{ "trProtonWeight" };
  auto f1_weight = Weight< std::vector<double> >{ "F1w" };
  auto f2_weight = Weight< std::vector<double> >{ "F2w" };
  auto f3_weight = Weight< std::vector<double> >{ "F3w" };
  auto tp_weight = Weight< std::vector<double> >{ "trTposW" };
  auto tn_weight = Weight< std::vector<double> >{ "trTnegW" };
  auto qvector_weight = Weight<double>{ "One" };

  auto start = std::chrono::steady_clock::now();
  auto handler = CorrelationHandler{ sampled_d, 100 };
  handler
    .AddComponentMeans( Vector<uvector_t>{ "proton" }, harmonics, proton_weight, proton_axes)
    .AddComponentCovariance( Vector<uvector_t>{ "proton" }, harmonics, proton_weight, proton_axes)

    .AddComponentMeans( Vector<uvector_t>{ "Tpos" }, harmonics, tp_weight, t_axes )
    .AddComponentCovariance( Vector<uvector_t>{ "Tpos" }, harmonics, tp_weight, t_axes )

    .AddComponentMeans( Vector<uvector_t>{ "Tneg" }, harmonics, tn_weight, t_axes )
    .AddComponentCovariance( Vector<uvector_t>{ "Tneg" }, harmonics, tn_weight, t_axes )

    .AddComponentMeans( Vector<qvector_t>{ "F1" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddComponentCovariance( Vector<qvector_t>{ "F1" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddComponentMeans( Vector<qvector_t>{ "F2" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddComponentCovariance( Vector<qvector_t>{ "F2" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)

    .AddComponentMeans( Vector<qvector_t>{ "F3" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)
    .AddComponentCovariance( Vector<qvector_t>{ "F3" }, qvector_harmonics, qvector_weight, qvector_axes, Qn::Stat::WeightType::REFERENCE)    
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