#ifndef BMN_ENV_H
#define BMN_ENV_H

#include <algorithm>
#include <vector>
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


const auto function_fhcal_x = []( const float FHCAL_Z ){
  return [FHCAL_Z]( ROOT::VecOps::RVec<std::vector<float>> vec_param ){
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
};
const auto function_fhcal_y = []( const float FHCAL_Z ){
  return [FHCAL_Z]( ROOT::VecOps::RVec<vector<float>> vec_param ){
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

const auto fhcal_weight_generator = []( const auto& layout ){
  return [&layout]( ROOT::VecOps::RVec<int> vec_mod_id, ROOT::VecOps::RVec<float> vec_mod_energy ){
    auto vec_weights = std::vector<double>( vec_mod_energy.size(), 0. );
    for( auto i=size_t{0}; i<vec_mod_id.size(); ++i ){
      auto mod_id = vec_mod_id[i];
      if( std::find( layout.begin(), layout.end(), mod_id ) == layout.end() )
        continue;
      vec_weights[i] = vec_mod_energy[i];
    }
    return vec_weights;
  };
};

const auto proton_weight = []( 
  std::vector<int> vec_is_proton, 
  std::vector<float> vec_efficiency, 
  std::vector<int> has_any_tof_hit,
  std::vector<float> vec_r,
  ROOT::VecOps::RVec<int> vec_nhits,
  ROOT::VecOps::RVec<float> vec_chi2,
  std::vector<float> vec_eta,
  std::vector<float> vec_fhcal_x,
  std::vector<float> vec_fhcal_y
  // "trIsProton", "trProtonEfficiency", "trHasAnyTofHit", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"
  ){
  auto weights = std::vector<double>( vec_is_proton.size(), 0.0 );
  for(auto i=size_t{}; i<vec_is_proton.size(); ++i){
    if( vec_is_proton[i] != 1 )
      continue;
    if( has_any_tof_hit[i] != 1 )
      continue;
    if( vec_r[i] > 5.0 )
      continue;
    if( vec_nhits[i] < 5 )
      continue;
    if( vec_chi2[i] > 5 )
      continue;
    if( vec_eta[i] > 2.5 )
      continue;
    if( -30 <  vec_fhcal_x[i]  && vec_fhcal_x[i] < 160 &&
        -60 < vec_fhcal_y[i] && vec_fhcal_y[i] < 60   )
      continue;
    
    weights[i] = vec_efficiency[i];
  }
  return weights;
};

const auto tpos_weight = []( 
  std::vector<float>       vec_eta, 
  ROOT::VecOps::RVec<float> vec_pT, 
  std::vector<float> vec_pq,
  std::vector<float> vec_r,
  ROOT::VecOps::RVec<int> vec_nhits,
  ROOT::VecOps::RVec<float> vec_chi2,
  std::vector<float> vec_fhcal_x,
  std::vector<float> vec_fhcal_y
  ){
  auto weights = std::vector<double>( vec_eta.size(), 0.0 );
  for(auto i=size_t{}; i<vec_eta.size(); ++i){
    if( vec_eta[i] < 2.0 )
      continue;
    if( vec_eta[i] > 3.0 )
      continue;
    if( vec_pT[i] < 0.2 )
      continue;
    if( vec_pq[i] < 0.0 )
      continue;
    if( vec_r[i] > 5.0 )
      continue;
    if( vec_nhits[i] < 5 )
      continue;
    if( vec_chi2[i] > 5 )
      continue;
    if( -30 <  vec_fhcal_x[i]  && vec_fhcal_x[i] < 160 &&
        -60 < vec_fhcal_y[i] && vec_fhcal_y[i] < 60   )
      continue;
    
    weights[i] = 1.0;
  }
  return weights;
};

const auto tneg_weight = []( 
  std::vector<float>       vec_eta, 
  ROOT::VecOps::RVec<float> vec_pT, 
  std::vector<float> vec_pq,
  std::vector<float> vec_r,
  ROOT::VecOps::RVec<int> vec_nhits,
  ROOT::VecOps::RVec<float> vec_chi2,
  std::vector<float> vec_fhcal_x,
  std::vector<float> vec_fhcal_y
  ){
  auto weights = std::vector<double>( vec_eta.size(), 0.0 );
  for(auto i=size_t{}; i<vec_eta.size(); ++i){
    if( vec_eta[i] < 1.5 )
      continue;
    if( vec_eta[i] > 4.0 )
      continue;
    if( vec_pT[i] < 0.2 )
      continue;
    if( vec_pq[i] > 0.0 )
      continue;
    if( vec_r[i] > 5.0 )
      continue;
    if( vec_nhits[i] < 5 )
      continue;
    if( vec_chi2[i] > 5 )
      continue;
    if( -30 <  vec_fhcal_x[i]  && vec_fhcal_x[i] < 160 &&
        -60 < vec_fhcal_y[i] && vec_fhcal_y[i] < 60   )
      continue;
    
    weights[i] = 1.0;
  }
  return weights;
};

const auto sim_f_weight = []( double eta1, double eta2 ){
  return [eta1, eta2]( std::vector<float> vec_eta, std::vector<float> vec_ekin, ROOT::VecOps::RVec<int> vec_m_id ){
    auto vec_weights = std::vector<double>( vec_eta.size(), 0 );
    for( auto i=size_t{0}; i<vec_eta.size(); ++i ){
      auto eta = vec_eta[i];
      auto m_id = vec_m_id[i];
      if( eta < eta1 )
        continue;
      if( eta > eta2 )
        continue;
      if( m_id != -1 )
        continue;

      vec_weights[i] = static_cast<double>(vec_ekin[i]);
    }
    return vec_weights;
  };
};

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

template<typename DataFrame>
const auto GenerateBmnExtendedTree(DataFrame& d, TH3* efficiency_histo){

  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

  auto dd=d
    .Define("track_multiplicity", "return static_cast<double>(trMom.size());")
    .Define("centrality", centrality_function, {"track_multiplicity"} )
    .Define("fhcalModPhi","std::vector<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
    .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
    .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
    .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
    .Define( "trDcaX", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(0) - vtxX ); } return vec_par; " )
    .Define( "trDcaY", " std::vector<float> vec_par; for( auto par : globalTrackParameters ){ vec_par.push_back( par.at(1) - vtxY ); } return vec_par; " )
    .Define( "trDcaR", dca_function, {"trDcaX", "trDcaY"} )
    .Define( "trFhcalX", function_fhcal_x(FHCAL_Z), {"trParamLast"} )
    .Define( "trFhcalY", function_fhcal_y(FHCAL_Z), {"trParamLast"} )
    .Define( "trChi2Ndf", " std::vector<float> vec_par; for( int i=0; i<trChi2.size(); ++i ){ vec_par.push_back( trChi2.at(i)/trNdf.at(i) ); } return vec_par; " )
    .Define( "trPx", " std::vector<float> px; for( auto mom : trMom ){ px.push_back( mom.Px() ); } return px; " )
    .Define( "trPy", " std::vector<float> py; for( auto mom : trMom ){ py.push_back( mom.Py() ); } return py; " )
    .Define( "trPhi","std::vector<float> phi;for(auto& mom : trMom) phi.push_back( mom.phi() ); return phi;")
    .Define( "trEta","std::vector<float> eta;for(auto& mom : trMom) eta.push_back( mom.Eta() ); return eta;")
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
    .Define( "simEta", "std::vector<float> simEta; for( auto mom : simMom ){ simEta.push_back( mom.Eta() ); } return simEta; " )
    .Define( "simEkin", "std::vector<float> simEkin; for( auto mom : simMom ){ simEkin.push_back( mom.E() - mom.M() ); } return simEkin; " ) 
    .Define( "simPhi", "std::vector<float> simPhi; for( auto mom : simMom ){ simPhi.push_back( mom.Phi() ); } return simPhi; " )
    .Define( "simF1w", sim_f_weight(4.4, 5.5), {"simEta", "simEkin", "simMotherId"} )
    .Define( "simF2w", sim_f_weight(3.9, 4.4), {"simEta", "simEkin", "simMotherId"} )
    .Define( "simF3w", sim_f_weight(3.1, 3.9), {"simEta", "simEkin", "simMotherId"} )

    .Define( "simIsProton", is_sim_particle(2212), {"simPdg", "simMotherId"} )
    .Define( "simProtonY", rapidity_generator(PROTON_M, Y_CM), {"simPz", "simP"} )
    
    .Define( "trIsProton", tr_is_particle, {"trSimIndex", "simIsProton"} )
    .Define( "trProtonWeight", proton_weight, {"trIsProton", "trProtonEfficiency", "trHasAnyTofHit", "trDcaR", "trStsNhits", "trStsChi2", "trEta", "trFhcalX", "trFhcalY"} )
    .Define( "trTposW", tpos_weight, {"trEta", "trPt", "pq", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"} )
    .Define( "trTnegW", tneg_weight, {"trEta", "trPt", "pq", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"} )
    .Define( "One", "return static_cast<double>(1.0)" )
    // .Range( 1000 )

    .Filter("vtxNtracks > 2")
  ;
  return dd;
}


#endif // BMN_ENV_H