//
// Created by Misha on 3/7/2023.
//

TVector3 GetFHCalPos(int iModule)
{
  const int Nmodules = 45;
  int xAxisSwitch = (iModule > Nmodules) ? 1 : -1;
  int module = (iModule < Nmodules) ? iModule : iModule - Nmodules;
  float x, y, z;
  z = (iModule > Nmodules) ? -320. : 320.;
  if (module >= 0 && module <= 4)
  {
    y = 45.;
    x = -1. * (module - 2) * 15.;
  }
  else if ((module >= 5) && (module <= 39))
  {
    y = (3 - (module + 2) / 7) * 15.;
    x = (3 - (module + 2) % 7) * 15.;
  }
  else if ((module >= 40) && (module <= 44))
  {
    y = -45.;
    x = -1. * (module - 42) * 15.;
  }
  TVector3 vec(x * xAxisSwitch, y, z);

  return vec;
}

void mpd_fixed_target_correct(std::string list, std::string collision_energy="2.5", std::string str_nucleus_mass="209"){

  const double T = std::stod( collision_energy );
  const double M = 0.938;
  const double GAMMA = (T + M) / M;
  const double BETA = sqrt(1 - (M * M) / (M + T) / (M + T));
  const double PZ = M * BETA * GAMMA;
  const double E = T + M;
  const double Y_BEAM = 0.25 * log((E + PZ) / (E - PZ));
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  std::vector<int> f1_modules = { 14, 15, 16, 21, 23, 28, 29, 30 };
  std::vector<int> f2_modules = { 6, 7, 8, 9, 10, 13, 17, 20, 24, 27, 31, 34, 35, 36, 37, 38 };
  std::vector<int> f3_modules = { 0, 1, 2, 3, 4, 5, 11, 12, 18, 19, 25, 26, 32, 33, 39, 40, 41, 42, 43, 44 };
  TStopwatch timer;
  timer.Start();
  std::string treename = "picodst";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define( "b_norm", [NUCLEUS_RADIUS](float b){ return b / NUCLEUS_RADIUS; },{"mcevent.fB"})
          .Define( "psi_rp", "mcevent.fPhiRP")
          .Define( "sim_y", [Y_BEAM]( ROOT::VecOps::RVec<float> vec_energy, ROOT::VecOps::RVec<TVector3> vec_momentum ){
            std::vector<float> vec_y;
            for( int i=0; i<vec_energy.size(); ++i ){
              auto E = vec_energy.at(i);
              auto pz = vec_momentum.at(i).Pz();
              auto y = 0.5 * log( (E+pz) / (E-pz) ) - Y_BEAM;
              vec_y.push_back(y);
            }
            return vec_y;
          }, { "mctracks.fEnergy", "mctracks.fMom" })
          .Define( "sim_pT", "std::vector<float> vec_pT; for( auto mom : mctracks.fMom ){ vec_pT.push_back( mom.Pt() ); } return vec_pT;" )
          .Define( "sim_phi", "std::vector<float> vec_phi; for( auto mom : mctracks.fMom ){ vec_phi.push_back( mom.Phi() ); } return vec_phi;" )
          .Define( "sim_pdg", "std::vector<int> vec_pdg; for( auto pdg : mctracks.fPdg ){ vec_pdg.push_back( pdg ); } return vec_pdg;" )
          .Define( "sim_mother_id", "std::vector<int> vec_mother_id; for( auto id : mctracks.fMotherId ){ vec_mother_id.push_back( id ); } return vec_mother_id;" )
          .Define( "tr_pT", "std::vector<float> vec_pT; for( auto mom : recotracks.fMom ){ vec_pT.push_back( mom.Pt() ) ;} return vec_pT;" )
          .Define( "tr_phi", "std::vector<float> vec_phi; for( auto mom : recotracks.fMom ){ vec_phi.push_back( mom.Phi() ); } return vec_phi;" )
          .Define( "tr_eta", "std::vector<float> vec_eta; for( auto mom : recotracks.fMom ){ vec_eta.push_back( mom.Eta() ); } return vec_eta;" )
          .Define( "tr_charge", "std::vector<short> vec_charge{ recotracks.fChargeSign.data(), recotracks.fChargeSign.data() + recotracks.fChargeSign.size() }; return vec_charge;" )
          .Define( "tr_pdg", []( ROOT::VecOps::RVec<Int_t> vec_tr_sim_id, ROOT::VecOps::RVec<Int_t> vec_mc_pdg, ROOT::VecOps::RVec<Int_t> vec_mother_id ){
            std::vector<int> vec_tru_pdg;
            vec_tru_pdg.reserve( vec_tr_sim_id.size() );
            for( auto sim_id : vec_tr_sim_id ){
              if( sim_id >= vec_mc_pdg.size() ){
                vec_tru_pdg.push_back(0);
                continue;
              }
              if( sim_id < 0 ){
                vec_tru_pdg.push_back(0);
                continue;
              }
              auto pdg = vec_mc_pdg.at(sim_id);
              auto mother_id = vec_mother_id.at(sim_id);
              if( mother_id != -1 ){
                vec_tru_pdg.push_back(0);
                continue;
              }
              vec_tru_pdg.push_back(pdg);
            }
            return vec_tru_pdg;
          }, {"recotracks.fId", "mctracks.fPdg", "mctracks.fMotherId"} )
          .Define( "tr_y", [Y_BEAM]( std::vector<int> vec_pdg, ROOT::VecOps::RVec<TVector3> vec_momentum ){
            std::vector<float> vec_y;
            vec_y.reserve( vec_pdg.size() );
            for( int i=0; i<vec_pdg.size(); ++i ){
              auto pdg = vec_pdg.at(i);
              if( pdg == 0 ){
                vec_y.push_back( -999. );
                continue;
              }
              auto pz = vec_momentum.at(i).Pz();
              auto p = vec_momentum.at(i).Mag();
              auto m = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
              auto E = sqrt( p*p + m*m );
              auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
              vec_y.push_back(y);
            }
            return vec_y;
          }, {"tr_pdg", "recotracks.fMom"} )
          .Define( "tr_nhits", "std::vector<int> vec_nhits; for( auto n : recotracks.fNhits ){ vec_nhits.push_back( n ); } return vec_nhits;" )
          .Define( "fhcal_module_id", "std::vector<int> module_id( FHCalModules.fEnergy.size() ); std::iota( module_id.begin(), module_id.end(), 0 ); return module_id;" )
          .Define( "fhcal_module_pos", []( std::vector<int> vec_id ){
            std::vector<TVector3> vec_module_pos{};
            for( auto id : vec_id )
              vec_module_pos.push_back( GetFHCalPos(id) );
            return vec_module_pos;
          }, {"fhcal_module_id"} )
          .Define( "fhcal_module_phi", "std::vector<float> vec_phi{}; for( auto pos : fhcal_module_pos ){ vec_phi.push_back( pos.Phi() ); }; return vec_phi;" )
          .Define( "fhcal_module_energy", "std::vector<float> vec_phi{FHCalModules.fEnergy.data(), FHCalModules.fEnergy.data()+FHCalModules.fEnergy.size()}; return vec_phi;" )
          .Filter("b_norm < 16."); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("b_norm|psi_rp"));
  correction_task.SetChannelVariables({std::regex("fhcal_module_(id|phi|energy)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr_(pT|y|phi|charge|pdg)"),
                                            std::regex("sim_(pT|y|phi|pdg|mother_id)")
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( { "b_norm", 10, 0, 2.0 } );

  VectorConfig f1( "F1", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f1.AddCut( "fhcal_module_id", [f1_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f2.AddCut( "fhcal_module_id", [f2_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcal_module_phi", "fhcal_module_energy", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f3.AddCut( "fhcal_module_id", [f3_modules](double mod_id){
    auto id = std::round(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  correction_task.AddVector(f3);

  std::vector<Qn::AxisD> proton_axes{
          { "tr_y", 15, -1.5, 1.5 },
          { "tr_pT", 15, 0.0, 1.5 },
  };

  VectorConfig proton( "proton", "tr_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  proton.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "q > 0" );
  correction_task.AddVector(proton);

  VectorConfig Tp( "Tp", "tr_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tp.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tp.AddCut( "tr_charge", [](double charge){ return charge > 0; }, "proton cut" );
  Tp.AddCut( "tr_y", [](double ycm){
    return -1.2 < ycm && ycm < 0.8;
    }, "Tp ycm cut" );
  correction_task.AddVector(Tp);

  VectorConfig Tneg( "Tpi", "tr_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "tr_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return abs(pdg_code) == 211;
  }, "pion cut" );
  Tneg.AddCut( "tr_y", [](double eta){
    return -1.5 < eta && eta < 0.2;
    }, "Tpi y cut" );
  correction_task.AddVector(Tneg);

  std::vector<Qn::AxisD> sim_proton_axes{
          { "tr_y", 15, -1.5, 1.5 },
          { "tr_pT", 15, 0.0, 1.5 },
  };

  VectorConfig tru_proton( "tru_proton", "sim_phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( {1, 2} );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( sim_proton_axes );
  tru_proton.AddCut( "sim_pdg", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == 2212;
  }, "tru_proton cut" );
  tru_proton.AddCut( "sim_mother_id", [](double pid){
    auto pdg_code = std::round(pid);
    return pdg_code == -1;
  }, "cut on primary" );
  correction_task.AddVector(tru_proton);

  VectorConfig psi_rp( "psi_rp", "psi_rp", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}