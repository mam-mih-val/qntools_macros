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

void run8_mc_proton_fill( std::string list, 
                          std::string str_effieciency_file,
                          std::string str_file_calib_recentering,
                          std::string calib_in_file="qa.root" ){

  std::cout << "starting execution" << std::endl;
  
  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

  auto file_calib_recentering = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open(str_file_calib_recentering.c_str(), "READ"), [](auto f){ f->Close(); } };
  
  Qn::DataContainerStatCollect* proton_avg_x{nullptr}; 
  Qn::DataContainerStatCollect* proton_avg_y{nullptr}; 
  
  if(file_calib_recentering) file_calib_recentering->GetObject( "proton_PLAIN.x1centrality", proton_avg_x );
  if(file_calib_recentering) file_calib_recentering->GetObject( "proton_PLAIN.y1centrality", proton_avg_y );
  
	const auto rapidity_generator = []( auto particle_m, auto y_cm ){
    return 
    [particle_m, y_cm]( std::vector<float> vec_pz, std::vector<float> vec_pq ){
      ROOT::VecOps::RVec<float> vec_y{};
      vec_y.reserve( vec_pz.size() );
      for( int i=0; i<vec_pz.size(); ++i ){
        auto pz = vec_pz.at(i);
        auto p = vec_pq.at(i);
        auto E = sqrt( p*p + particle_m*particle_m );
        auto y = 0.5 * log( ( E + pz ) / ( E - pz ) ) - y_cm;
        vec_y.push_back( y );
      }
      return vec_y;
    };
  };

  const auto function_fhcal_x = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<std::vector<float>> vec_param ){
      std::vector<float> vec_x{};
      vec_x.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto x = par.at(0);
        auto z = par.at(2);
        auto tx = par.at(3);
        auto dz = FHCAL_Z - z;
        auto dx = tx * dz;
        vec_x.push_back( x+dx );
      }
      return vec_x;
    };
  const auto function_fhcal_y = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<vector<float>> vec_param ){
      std::vector<float> vec_y{};
      vec_y.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto y = par.at(1);
        auto z = par.at(2);
        auto ty = par.at(4);
        auto dz = FHCAL_Z - z;
        auto dy = ty * dz;
        vec_y.push_back( y+dy );
      }
      return vec_y;
    };

  const auto centrality_function = 
  []
  (double multiplicity){
      float centrality;
      std::vector<float> centrality_percentage{ 0, 10, 20, 30, 40, 50, 60, 70, 100 };
      std::vector<int> multiplicity_edges{ 206, 98, 70, 49, 34, 22, 14, 8, 0  };
      if( multiplicity > multiplicity_edges[0] )
        return -1.0f;
      int idx = 0;
      float bin_edge = multiplicity_edges[idx];
      while( multiplicity < bin_edge &&
        idx < multiplicity_edges.size()-1 ){
        idx++;
        bin_edge = multiplicity_edges[idx];
      }
      centrality = (centrality_percentage[idx-1] + centrality_percentage[idx])/2.0f;
      return centrality;
  };
  
  const auto dca_function = [](std::vector<float> vec_x, std::vector<float> vec_y){
    std::vector<float> vec_r{};
    vec_r.reserve(vec_x.size());
    for (int i=0; i<vec_x.size(); ++i) {
      auto x = vec_x.at(i);
      auto y = vec_y.at(i);
      auto r = std::sqrt( x*x + y*y );
      vec_r.push_back(r);
    }
    return vec_r;
  };

  const auto is_sim_particle = []( int pdg_code ) {
    return [pdg_code]( ROOT::VecOps::RVec<int> vec_pdg, ROOT::VecOps::RVec<int> vec_mother_id ){
      auto vec_is = std::vector<int>( vec_pdg.size(), 0 );
      for( auto i=size_t{0}; i < vec_pdg.size(); ++i ){
        auto pdg = vec_pdg[i];
        auto m_id = vec_mother_id[i];
        if( pdg != pdg_code )
          continue;
        if( m_id != -1 )
          continue;
        vec_is[i] = 1;
      }
      return vec_is;
    };
  };

  const auto tr_is_particle = []( ROOT::VecOps::RVec<int> vec_sim_idx, std::vector<int> vec_is_sim_particle ){
    auto vec_is = std::vector<int>( vec_sim_idx.size(), 0 );
    for( auto i=size_t{0}; i<vec_sim_idx.size(); ++i ){
      auto idx = vec_sim_idx[i];
      if( idx > vec_is_sim_particle.size() )
        continue;
      if( idx < 0 )
        continue;
      vec_is[i] = vec_is_sim_particle[idx];
    }
    return vec_is;
  };

  const auto tr_has_tof_hit = []( ROOT::VecOps::RVec<double> vec_beta ){
    auto vec_has = std::vector<int>( vec_beta.size(), 0 );
    for( auto i = size_t{0}; i < vec_beta.size(); ++i ){
      if( vec_beta[i] < -9. )
              continue;
      vec_has[i] = 1;
    }
    return vec_has;
  };

  const auto tr_has_any_tof_hit = []( std::vector<int> vec_is_400, std::vector<int> vec_is_700 ){
    auto vec_has = std::vector<int>( vec_is_400.size(), 0 );
    for( auto i=size_t{0}; i<vec_is_400.size(); ++i ){
      if( vec_is_400[i] == 0 && vec_is_700[i] == 0 )
        continue;
      vec_has[i] = 1;
    }
    return vec_has;
  };

  const auto weight_generator = []( auto efficiency_map ){
    return [efficiency_map](std::vector<float> vec_px, std::vector<float> vec_py, std::vector<float> vec_pz){
      if( !efficiency_map ){
          return std::vector<float>(vec_px.size(), 1);
        }
      auto vec_weight = std::vector<float>( vec_px.size(), 0.0 );
      for( int i=0; i<vec_px.size(); ++i ){
        auto px = vec_px[i];
        auto py = vec_py[i];
        auto pz = vec_pz[i];

        auto px_bin = efficiency_map->GetXaxis()->FindBin( px );
        auto py_bin = efficiency_map->GetYaxis()->FindBin( py );
        auto pz_bin = efficiency_map->GetZaxis()->FindBin( pz );
        
        auto efficiency = efficiency_map->GetBinContent( px_bin, py_bin, pz_bin );
        if( efficiency < 1e-2 )
          continue;
        auto weight = 1.0 / efficiency;
        vec_weight[i] = static_cast<float>(weight);
      }
      return vec_weight;
    };
  };

  const auto proton_weight = []( 
    std::vector<int> vec_is_proton, 
    std::vector<float> vec_efficiency, 
    std::vector<int> has_any_tof_hit,
    std::vector<float> vec_r,
    ROOT::VecOps::RVec<int> vec_nhits,
    ROOT::VecOps::RVec<float> vec_chi2,
    std::vector<float> vec_fhcal_x,
    std::vector<float> vec_fhcal_y
   ){
    auto weights = std::vector<double>( vec_is_proton.size(), 0.0 );
    for(auto i=size_t{}; i<vec_is_proton.size(); ++i){
      if( vec_is_proton[i] != 1 )
        continue;
      if( has_any_tof_hit[i] != 1 )
        continue;
      if( vec_r[i] < 5.0 )
        continue;
      if( vec_nhits[i] < 5 )
        continue;
      if( vec_chi2[i] > 5 )
        continue;
      if( vec_fhcal_x[i] < -30 || vec_fhcal_x[i] > 160  )
        continue;
      if( vec_fhcal_y[i] < -60 || vec_fhcal_y[i] > 60  )
        continue;

      weights[i] = vec_efficiency[i];
    }
    return weights;
  };


  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH3* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("h3_efficiency_2212_good", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No efficiency for both tof was found in file " << str_effieciency_file << "\n";

  std::vector<int> f1_modules = {
    6,  7,  8,
    11, 12, 13,
    16,     17,
    20, 21, 22, 
    25, 26, 27
  };
  std::vector<int> f2_modules = {
    0,  1,  2,  3,  4,
    5,              9,
    10,             14,
    15,             18,
    19,             23,
    24,             28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_modules = {
    35,                 44,
    37,                 46, 
    39,                 48, 
    41,                 50,
    43,                 52
  };

  std::vector<int> f4_modules = {
    34,                     45,
    36,                     47, 
    38,                     49, 
    40,                     51,
    42,                     53
  };

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  
  auto dd=d
    .Define("track_multiplicity", "return static_cast<double>(trMom.size());")
    .Define("centrality", centrality_function, {"track_multiplicity"} )
    .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
    .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
    .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
    .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
    .Define( "trDcaX", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(0) - vtxX ); } return vec_par; " )
    .Define( "trDcaY", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(1) - vtxY ); } return vec_par; " )
    .Define( "trDcaR", dca_function, {"trDcaX", "trDcaY"} )
    .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
    .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
    .Define( "trChi2Ndf", " std::vector<float> vec_par; for( int i=0; i<trChi2.size(); ++i ){ vec_par.push_back( trChi2.at(i)/trNdf.at(i) ); } return vec_par; " )
    .Define( "trPx", " std::vector<float> px; for( auto mom : trMom ){ px.push_back( mom.Px() ); } return px; " )
    .Define( "trPy", " std::vector<float> py; for( auto mom : trMom ){ py.push_back( mom.Py() ); } return py; " )
    .Define("trPhi","std::vector<float> phi;for(auto& mom : trMom) phi.push_back( mom.phi() ); return phi;")
    .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
    .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P() / trCharge.at(i) ); } return pq;" )
    .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"pz", "pq"} )
    .Define( "trProtonEfficiency", weight_generator(efficiency_histo), {"trPx", "trPy", "pz"} )
    .Define( "trHasTof400Hit", tr_has_tof_hit, { "trBetaTof400" } )
    .Define( "trHasTof700Hit", tr_has_tof_hit, { "trBetaTof700" } )
    .Define( "trHasAnyTofHit", tr_has_any_tof_hit, { "trHasTof400Hit", "trHasTof700Hit" } )


    .Alias("trStsNhits", "stsTrackNhits")
    .Alias("trStsChi2", "stsTrackChi2Ndf")


    .Define( "simP", "std::vector<float> simP; for( auto mom : simMom ){ simP.push_back( mom.P() ); } return simP; " )
    .Define( "simPt", "std::vector<float> simPt; for( auto mom : simMom ){ simPt.push_back( mom.Pt() ); } return simPt; " )
    .Define( "simPz", "std::vector<float> simPz; for( auto mom : simMom ){ simPz.push_back( mom.Pz() ); } return simPz; " )
    .Define( "simPhi", "std::vector<float> simPhi; for( auto mom : simMom ){ simPhi.push_back( mom.Phi() ); } return simPhi; " )
          
    .Define( "simIsProton", is_sim_particle(2212), {"simPdg", "simMotherId"} )
    .Define( "simProtonY", rapidity_generator(PROTON_M, Y_CM), {"simPz", "simP"} )
    
    .Define( "trIsProton", tr_is_particle, {"trSimIndex", "simIsProton"} )
    .Define( "trProtonWeight", proton_weight, {"trIsProton", "trProtonEfficiency", "trHasAnyTofHit", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"} )

    .Define( "x1", u_generator(1, [](double x){ return std::cos(x); }), {"trPhi"} )
    // .Range( 1000 )

    .Filter("vtxNtracks > 2")      
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