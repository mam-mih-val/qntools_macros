//
// Created by Misha on 3/7/2023.
//

#include <cassert>
#include <cmath>
#include <random>
#include <vector>

void run8_mc_proton_correct( std::string list, 
                          std::string str_effieciency_file,
                          std::string calib_in_file="qa.root" ){

  std::cout << "starting execution" << std::endl;
  
  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm
  
	const auto rapidity_generator = []( auto particle_m, auto y_cm ){
    return 
    [particle_m, y_cm]( std::vector<float> vec_pz, std::vector<float> vec_pq ){
      std::vector<float> vec_y{};
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
    return [efficiency_map](std::vector<float> vec_y, ROOT::VecOps::RVec<float> vec_pT){
      if( !efficiency_map ){
          return std::vector<float>(vec_y.size(), 1);
        }
      std::vector<float> vec_weight{};
      vec_weight.reserve(vec_y.size());
      for( int i=0; i<vec_y.size(); ++i ){
        auto pT = vec_pT.at(i);
        auto y = vec_y.at(i);
        auto y_bin = efficiency_map->GetXaxis()->FindBin( y );
        auto pT_bin = efficiency_map->GetYaxis()->FindBin( pT );
        auto efficiency = efficiency_map->GetBinContent( y_bin, pT_bin );
        auto weight = 5e-2 < efficiency && efficiency < 1 ? 1.0 / efficiency : 0.0;
        vec_weight.push_back( weight );
      }
      return vec_weight;
    };
  };

  auto device = std::random_device{};
  auto engine = std::mt19937{ device() };
  auto distribution = std::uniform_real_distribution<double>{ 0.0, 1.0 };

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2D* efficiency_histo{nullptr};
  
  effieciency_file->GetObject("efficiency_2212_tof", efficiency_histo);
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
    .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
    .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P() / trCharge.at(i) ); } return pq;" )
    .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"pz", "pq"} )
    .Define( "trProtonEfficiency", weight_generator(efficiency_histo), {"trProtonY", "trPt"} )
    .Define( "trHasTof400Hit", tr_has_tof_hit, { "trBetaTof400" } )
    .Define( "trHasTof700Hit", tr_has_tof_hit, { "trBetaTof700" } )
    .Define( "trHasAnyTofHit", tr_has_any_tof_hit, { "trHasTof400Hit", "trHasTof700Hit" } )

    .Alias("trStsNhits", "stsTrackNhits")
    .Alias("trStsChi2", "stsTrackChi2Ndf")
    .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
    .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")

    .Define( "simP", "std::vector<float> simP; for( auto mom : simMom ){ simP.push_back( mom.P() ); } return simP; " )
    .Define( "simPt", "std::vector<float> simPt; for( auto mom : simMom ){ simPt.push_back( mom.Pt() ); } return simPt; " )
    .Define( "simPz", "std::vector<float> simPz; for( auto mom : simMom ){ simPz.push_back( mom.Pz() ); } return simPz; " )
    .Define( "simPhi", "std::vector<float> simPhi; for( auto mom : simMom ){ simPhi.push_back( mom.Phi() ); } return simPhi; " )
          
    .Define( "simIsProton", is_sim_particle(2212), {"simPdg", "simMotherId"} )
    .Define( "simProtonY", rapidity_generator(PROTON_M, Y_CM), {"simPz", "simP"} )
    
    .Define( "trIsProton", tr_is_particle, {"trSimIndex", "simIsProton"} )

    .Filter("vtxNtracks > 2")      
  ; 

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("centrality|psiRP"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                      std::regex("tr(Pt|Px|Py|Eta|Phi|Charge|ProtonY|DcaR|ProtonEfficiency|FhcalX|FhcalY|StsNhits|StsChi2|HasAnyTofHit|IsProton)"),
                                      std::regex("sim(Pt|ProtonY|Phi|IsProton)"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 6, 0, 60} );
  // correction_task.AddEventAxis( { "runId", 12, 7100, 8300 } );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( { 1, 2, 3, 4, 5, 6, 7, 8 } );
  f1.SetCorrections( {CORRECTION::PLAIN } );
  f1.AddCut( "fhcalModId", [&f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2, 3, 4, 5, 6, 7, 8 } );
  f2.SetCorrections( {CORRECTION::PLAIN } );
  f2.AddCut( "fhcalModId", [&f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2, 3, 4, 5, 6, 7, 8 } );
  f3.SetCorrections( {CORRECTION::PLAIN } );
  f3.AddCut( "fhcalModId", [&f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  VectorConfig f4( "F4", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f4.SetHarmonicArray( {1, 2, 3, 4, 5, 6, 7, 8 } );
  f4.SetCorrections( {CORRECTION::PLAIN } );
  f4.AddCut( "fhcalModId", [&f4_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f4_modules.begin(), f4_modules.end(), id) != f4_modules.end();
    }, "F3 Cut" );
  f4.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f4);

  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2, 3, 4, 5, 6, 7, 8 } );
  Tneg.SetCorrections( {CORRECTION::PLAIN } );
  Tneg.AddCut( "trCharge", [](double charge){
    return charge < 0.0;
    }, "charge" );
  Tneg.AddCut( "trEta", [](double eta){
    return 1.5 < eta && eta < 4.0;
    }, "eta cut" );
  Tneg.AddCut( "trPt", [](double pT){
    return pT > 0.2;
    }, "pT cut" );
  Tneg.AddCut( "trFhcalX", [](double pos){
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tneg.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tneg);

  VectorConfig Tpos( "Tpos", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tpos.SetHarmonicArray( {1, 2, 3, 4, 5, 6, 7, 8 } );
  Tpos.SetCorrections( {CORRECTION::PLAIN } );
  Tpos.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  Tpos.AddCut( "trEta", [](double eta){
    return 2.0 < eta && eta < 3.0;
  }, "eta cut" );
  Tpos.AddCut( "trPt", [](double pT){
    return pT > 0.2;
  }, "pT cut" );
  Tpos.AddCut( "trFhcalX", [](double pos){
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tpos.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tpos);

  std::vector<Qn::AxisD> proton_axes{
        { "trProtonY", 6, 0.0, 1.2 },
        { "trPt", 5, 0.0, 2.0 },
  };
  
  VectorConfig proton( "proton", "trPhi", "trProtonEfficiency", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( { 1, 2, 3, 4, 5, 6, 7, 8 } );
  proton.SetCorrections( { CORRECTION::PLAIN } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "trIsProton", [](double is){
    return 0.9 < is && is < 1.1;
  }, "cut on if it is proton from MC info" );
  proton.AddCut( "trHasAnyTofHit", [](double has){
    return 0.9 < has && has < 1.1;
  }, "cut on if track has any tof hits" );
  proton.AddCut( "trFhcalX", [](double pos){
    return pos < -30.0 || pos > 160;
  }, "cut on x-pos in fhcal plane" );
  proton.AddCut( "trFhcalY", [](double pos){
    return pos < -60.0 || pos > 60;
  }, "cut on y-pos in fhcal plane" );
  proton.AddCut( "trStsNhits", [](double nhits){
    return nhits > 5.5;
  }, "cut on fake tracks" );
  proton.AddCut( "trDcaR", [](double dca){
    return dca < 5.0;
  }, "DCA cut" );
  proton.AddCut( "trStsChi2", [](double chi2){
    return chi2 < 5.0;
  }, "cut on chi2 in sts" );
  correction_task.AddVector(proton);

  std::vector<Qn::AxisD> sim_proton_axes{
        { "simProtonY", 6, 0.0, 1.2 },
        { "simPt", 5, 0.0, 2.0 },
  };

  VectorConfig sim_proton( "tru_proton", "simPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  sim_proton.SetHarmonicArray( { 1, 2, 3, 4, 5, 6, 7, 8 } );
  sim_proton.SetCorrections( { CORRECTION::PLAIN } );
  sim_proton.SetCorrectionAxes( sim_proton_axes );
  sim_proton.AddCut( "simIsProton", [](double is){
    return 0.9 < is && is < 1.1;
  }, "cut on if it is proton from MC info" );
  correction_task.AddVector(sim_proton);

  VectorConfig psi_rp( "psi_rp", "psiRP", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( { 1, 2 } );
  psi_rp.SetCorrections( { CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);


  std::cout << "Initialized" << std::endl;

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}