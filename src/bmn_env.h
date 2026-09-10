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

auto vtx_correction_generator = 
[]( TGraphErrors* g1_calib ){
  return [g1_calib](double _vtx, UInt_t _runId){return _vtx - g1_calib->Eval( static_cast<double>(_runId) ); };
};
auto ref_mult_generator =
[]( TGraphErrors* g1_calib ){
  return [g1_calib](unsigned long _mult, UInt_t _runId){ return (_mult * g1_calib->Eval( static_cast<double>(_runId) )); };
};

const auto n_sigma_generator = []( auto f1_mean, auto f1_sigma ){
  return 
  [ f1_mean, f1_sigma ]
  ( std::vector<float> vec_pq, ROOT::VecOps::RVec<float> vec_m2 ){
      auto vec_n_sigma = std::vector<float>( vec_pq.size(), 999. );
      for( size_t i=0; i < vec_pq.size(); ++i ){
        auto m2 = vec_m2.at(i);
        auto pq = vec_pq.at(i);
        auto mean = f1_mean->Eval(pq);
        auto sigma = f1_sigma->Eval(pq);
        auto n_sigma = fabs( m2 - mean ) / sigma;
        if( pq < 1.0 )
          continue;
        if( pq > 8.0 )
          continue;
        vec_n_sigma[i] = n_sigma;
      }
    return vec_n_sigma;
  };
};

const auto trWeightFunction = []( auto efficiency_eta_pT_phi_run_id ){
  return [efficiency_eta_pT_phi_run_id]( 
  ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >> vec_mom, UInt_t run_id 
  ){  
    std::vector<float> vec_efficiency(vec_mom.size(), 0.f);
    for( int i=0; i<vec_mom.size(); i++ ){
      auto p = vec_mom.at(i);
      auto pT = p.Pt();
      auto eta = p.Eta();
      auto phi = p.Phi();
      auto coord = std::vector<double>{ eta, pT, phi, static_cast<double>(run_id) };
      auto bin = efficiency_eta_pT_phi_run_id->GetBin( coord.data() );
      if( bin <= 0 )
        continue;
      if( bin > efficiency_eta_pT_phi_run_id->GetNbins() )
        continue;
      auto weight = efficiency_eta_pT_phi_run_id->GetBinContent(bin);
      vec_efficiency[i] = weight;
    }
    return vec_efficiency;
  };
};
const auto n_sigma_particle_function = 
[]
( std::vector<float> n_sigma_400, 
  std::vector<float> n_sigma_700 ){
    std::vector<float> vec_n_simga{};
    vec_n_simga.reserve( n_sigma_400.size() );
    for( int i=0; i<n_sigma_400.size(); ++i ){ 
      vec_n_simga.push_back( std::min( n_sigma_400.at(i), n_sigma_700.at(i) ) ); }
    return vec_n_simga;
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
  return [efficiency_map](std::vector<float> vec_p, std::vector<float> vec_eta, std::vector<float> vec_phi){
    if( !efficiency_map ){
        return std::vector<float>(vec_p.size(), 1);
      }
    auto vec_weight = std::vector<float>( vec_p.size(), 0.0 );
    for( int i=0; i<vec_p.size(); ++i ){
      auto p = vec_p[i];
      auto eta = vec_eta[i];
      auto phi = vec_phi[i];

      auto eta_bin = efficiency_map->GetXaxis()->FindBin( eta );
      auto p_bin = efficiency_map->GetYaxis()->FindBin( p );
      auto phi_bin = efficiency_map->GetZaxis()->FindBin( phi );
      
      auto efficiency = efficiency_map->GetBinContent( eta_bin, p_bin, phi_bin );
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
    // if( vec_eta[i] > 3.0 )
    //   continue;
    // if( -30 <  vec_fhcal_x[i]  && vec_fhcal_x[i] < 160 &&
    //     -60 < vec_fhcal_y[i] && vec_fhcal_y[i] < 60   )
    //   continue;
    
    weights[i] = vec_efficiency[i];
  }
  return weights;
};

const auto proton_weight_data = []( 
  std::vector<float> vec_n_sigma_proton, 
  std::vector<float> vec_efficiency, 
  std::vector<float> vec_r,
  ROOT::VecOps::RVec<int> vec_nhits,
  ROOT::VecOps::RVec<float> vec_chi2,
  std::vector<float> vec_eta,
  std::vector<float> vec_fhcal_x,
  std::vector<float> vec_fhcal_y
  // "trIsProton", "trProtonEfficiency", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"
  ){
  auto weights = std::vector<double>( vec_n_sigma_proton.size(), 0.0 );
  for(auto i=size_t{}; i<vec_n_sigma_proton.size(); ++i){
    if( vec_n_sigma_proton[i] > 3.0 )
      continue;
    if( vec_r[i] > 5.0 )
      continue;
    if( vec_nhits[i] < 5 )
      continue;
    if( vec_chi2[i] > 5 )
      continue;
    // if( vec_eta[i] > 3.0 )
    //   continue;
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

auto charge_function = []( std::vector<float> vec_pq, ROOT::VecOps::RVec<float> vec_dedx ){
  auto bb_body = [](Double_t *x, Double_t *par) {
  // x[0] = p/q (momentum over charge) in GeV/c
    Double_t p = x[0];
    Double_t m = par[0];
    Double_t A = par[1];
    Double_t delta = par[2];
    Double_t norm = par[3];
    Double_t me = 0.000511;

    Double_t bg = p / m;
    Double_t beta = bg / sqrt(1. + bg * bg);
    Double_t gamma_fac = sqrt(1. + bg * bg);
    Double_t ekin_max = 2*pow(m * gamma_fac + me, 2) * me / (m * m);
    Double_t dEdx = A * (1. / (beta * beta)) * (-5.296 + log(bg) + log(ekin_max) - 2 * beta * beta - delta);

    return norm * dEdx;
  };
  auto f1_bethebloch_d = new TF1("fBetheBloch_d", bb_body, 0.1, 10, 4);
  f1_bethebloch_d->SetParameter(0, 2.26);
  f1_bethebloch_d->SetParameter(1, -1.64);
  f1_bethebloch_d->SetParameter(2, 35.73);
  f1_bethebloch_d->SetParameter(3, 2.02);

  auto vec_q = std::vector<float>( vec_pq.size(), -1.f );
  for( auto i = size_t{0}; i<vec_pq.size(); ++i ){
    auto pq = vec_pq.at(i);
    auto dedx = vec_dedx.at(i);
    if( pq < 0 ) 
      continue; 
    auto cut = f1_bethebloch_d->Eval( pq );
    vec_q.at(i) =  dedx < cut ? 1.f : 2.f;
  }
  return vec_q;
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
const auto GenerateBmnExtendedTreeMC(DataFrame& d, TH3* efficiency_histo){

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
    .Define( "trProtonEfficiency", weight_generator(efficiency_histo), {"pq", "trEta", "trPhi"} )
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

struct DataCalibration{
  THn* efficiency_eta_pT_phi_run_id{nullptr}; 
  TH3* efficiency_histo{nullptr};
  TF1* f1_2212_m_400{nullptr};
  TF1* f1_2212_s_400{nullptr};
  TF1* f1_2212_m_700{nullptr};
  TF1* f1_2212_s_700{nullptr};
  TGraphErrors* g1_FitVtxX{nullptr};
  TGraphErrors* g1_FitVtxY{nullptr};
  TGraphErrors* g1_FitVtxZ{nullptr};
  TGraphErrors* g1_FitRunIdFactor{nullptr};
  std::vector<int> selected_runs{};
  std::vector<int> rejected_runs{};
};

template<typename DataFrame>
const auto GenerateBmnExtendedTreeData(DataFrame& d, const DataCalibration& calibration ){

  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

  std::vector<int> physical_runs{ 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199 };
  std::vector<int> bad_runs{7313, 7415, 7417, 7435, 7469, 7517, 7519, 7520, 7537, 7575, 7604, 7630, 7657, 7659, 7679, 7681, 7705, 7735, 7843, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7907, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7954, 7955, 8031, 8032, 8033, 8115, 8121, 8167, 8201, 8204, 8205, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8247, 8265, 8266, 8267, 8281, 8289};

  std::cout << "Preparing the RDF" << std::endl;
  auto dd=d
    .Define( "vtxXcorr", vtx_correction_generator(calibration.g1_FitVtxX), {"vtxX","runId"})
    .Define( "vtxYcorr", vtx_correction_generator(calibration.g1_FitVtxY), {"vtxY","runId"})
    .Define( "vtxZcorr", vtx_correction_generator(calibration.g1_FitVtxY), {"vtxZ","runId"})
    .Define( "vtxRcorr", "return sqrt(vtxXcorr*vtxXcorr + vtxYcorr*vtxYcorr);" )
    .Define("track_multiplicity", "return trMom.size();")
    .Define( "ref_multiplicity", ref_mult_generator( calibration.g1_FitRunIdFactor ), {"track_multiplicity","runId"} )
    .Define("stsNdigits","return stsDigits.size()" )
    .Define("centrality", centrality_function, {"ref_multiplicity"} )
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
    .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
    .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P() / trCharge.at(i) ); } return pq;" )
    .Define( "trQ", charge_function, {"pq", "trEnergyLoss"} )
    .Define( "trNsigmaProton400", n_sigma_generator(calibration.f1_2212_m_400, calibration.f1_2212_s_400), { "pq", "trM2Tof400" } )
    .Define( "trNsigmaProton700", n_sigma_generator(calibration.f1_2212_m_700, calibration.f1_2212_s_700), { "pq", "trM2Tof700"  } )
    .Define( "trNsigmaProton", n_sigma_particle_function, {"trNsigmaProton400", "trNsigmaProton700"} )
    .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"pz", "pq"} )
    .Define("trEta","std::vector<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
    .Define("trPhi","std::vector<float> phi; for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
    .Define( "trWeight", trWeightFunction(calibration.efficiency_eta_pT_phi_run_id), {"trMom", "runId"} )
    .Define( "trProtonEfficiency", weight_generator(calibration.efficiency_histo), {"pq", "trEta", "trPhi"} )
    // .Define( "trProtonEfficiencyTof400", weight_generator(efficiency_tof400), {"trProtonY", "trPt"} )
    // .Define( "trProtonEfficiencyTof700", weight_generator(efficiency_tof700), {"trProtonY", "trPt"} )
    .Alias("trStsNhits", "stsTrackNhits")
    .Alias("trStsChi2", "stsTrackChi2Ndf")
    .Define( "trProtonEffEff", "std::vector<float> weights{}; for( auto i=size_t{0}; i<trWeight.size(); ++i ){ weights.push_back( trWeight[i]*trProtonEfficiency[i] ); } return weights;" )
    .Define( "trProtonWeight", proton_weight_data, {"trNsigmaProton", "trProtonEffEff", "trDcaR", "trStsNhits", "trStsChi2", "trEta", "trFhcalX", "trFhcalY"} )
    .Define( "trTposW", tpos_weight, {"trEta", "trPt", "pq", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"} )
    .Define( "trTnegW", tneg_weight, {"trEta", "trPt", "pq", "trDcaR", "trStsNhits", "trStsChi2", "trFhcalX", "trFhcalY"} )
    .Filter([&calibration]( UInt_t run_id ){ 
      if( std::find( calibration.selected_runs.begin(), calibration.selected_runs.end(), run_id) == calibration.selected_runs.end() )
        return false;
      if( std::find( calibration.rejected_runs.begin(), calibration.rejected_runs.end(), run_id) != calibration.rejected_runs.end() )
        return false;
      return true;
    }, {"runId"} )
    .Filter("runId < 8312")
    .Filter( []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }, {"triggerMapAR"} )
    .Filter("vtxNtracks > 2")
    .Filter("fabs(vtxRcorr)<1.5")
    .Filter("fabs(vtxZcorr)<1.0")
    .Filter("noPileup == 1")
    .Define( "One", "return static_cast<double>(1.0)" )
    .Range( 1000 )

  ;

  return dd;
}


#endif // BMN_ENV_H