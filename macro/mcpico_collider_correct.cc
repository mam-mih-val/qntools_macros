
//
// Created by Misha on 3/7/2023.
//

void mcpico_collider_correct(std::string list, std::string str_sqrt_snn="2.4", std::string str_nucleus_mass="197"){
  const double ETA_MIN = -1.0;
  const double ETA_MAX = 1.0;

  const double sqrt_snn = std::stod(str_sqrt_snn);
  const double M = 0.938;
  const double T = sqrt_snn * sqrt_snn/ 2 / M - 2*M;
  const double GAMMA = (T + M) / M;
  const double BETA = sqrt(1 - (M * M) / (M + T) / (M + T));
  const double PZ = M * BETA * GAMMA;
  const double E = T + M;
  const double Y_BEAM = 0.5 * log((E + PZ) / (E - PZ)) / 2.0;
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  std::cout << "Energy: " << sqrt_snn << "; Y beam: " << Y_BEAM << std::endl;
  std::cout << "Nucleus: " << nucleus_mass << "; Radius: " << NUCLEUS_RADIUS << std::endl;

  TStopwatch timer;
  timer.Start();
  std::string treename = "mctree";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define( "b_norm", [NUCLEUS_RADIUS](float b){
            return b/NUCLEUS_RADIUS;
          }, {"bimp"} )
          .Define( "reaction_plane", [](){
            std::random_device rnd_device;
            std::mt19937 generator(rnd_device());
            std::uniform_real_distribution<float> distribution(-M_PI,M_PI); // distribution in range [0, 1]
            return distribution(generator);
          }, {} )
          .Define("phi", []( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py, float psi_rp ){
            ROOT::VecOps::RVec<float> vec_phi;
            for( int i=0; i<vec_px.size(); ++i ){
              auto px = vec_px.at(i);
              auto py = vec_py.at(i);
              auto phi = atan2(py, px);
              phi+=psi_rp;
              vec_phi.push_back( phi );
            }
            return vec_phi;
          },{ "momx", "momy", "reaction_plane" })
          .Define("pT", []( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py ){
            ROOT::VecOps::RVec<float> vec_pT;
            for( int i=0; i<vec_px.size(); ++i ){
              auto px = vec_px.at(i);
              auto py = vec_py.at(i);
              vec_pT.push_back( sqrt(px*px + py*py) );
            }
            return vec_pT;
          },{ "momx", "momy", })
          .Define("y", []( ROOT::VecOps::RVec<float> vec_pz, ROOT::VecOps::RVec<float> vec_E ){
            ROOT::VecOps::RVec<float> vec_y;
            for( int i=0; i<vec_pz.size(); ++i ){
              auto pz = vec_pz.at(i);
              auto E = vec_E.at(i);
              auto y = 0.5 * log( (E+pz)/(E-pz) );
              vec_y.push_back( y );
            }
            return vec_y;
          },{ "momz", "ene", })
          .Define("eta_cm", []( ROOT::VecOps::RVec<float> vec_y, ROOT::VecOps::RVec<float> vec_pT, ROOT::VecOps::RVec<int> vec_pid ){
            ROOT::VecOps::RVec<float> vec_eta;
            for( int i=0; i<vec_y.size(); ++i ){
              auto y = vec_y.at(i);
              auto pT = vec_pT.at(i);
              auto pid = vec_pid.at(i);
              auto mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
              auto mT = sqrt( pT*pT + mass*mass );
              auto eta_cm = asinh( mT/pT * sinh( y ) );
              vec_eta.push_back(eta_cm);
            }
            return vec_eta;
          },{ "y", "pT", "pdg", })
          .Filter("for(auto x : phi){ if( std::isnan(x) ) return false; } return true;")
          .Filter("for(auto x : y){ if( std::isnan(x) ) return false; } return true;")
          ;

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables( std::regex("bimp|b_norm|reaction_plane|psi12") );
  correction_task.SetTrackVariables({
    std::regex("phi|pT|y|eta_cm|pdg|charge|")
  });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"b_norm", 20, 0, 2} );

  VectorConfig fp( "Fp", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  fp.SetHarmonicArray( {1, 2} );
  fp.SetCorrections( {CORRECTION::PLAIN } );
  fp.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2112 || pdg_code == 2212;
  }, "proton cut" );
  fp.AddCut( "eta_cm", [](double eta){
    return 2.0 < eta && eta < 5.0;
    }, "fp Cut" );
  correction_task.AddVector(fp);
  
  VectorConfig fn( "Fn", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  fn.SetHarmonicArray( {1, 2} );
  fn.SetCorrections( {CORRECTION::PLAIN } );
  fn.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2112 || pdg_code == 2212;
  }, "proton cut" );
  fn.AddCut( "eta_cm", [](double eta){
    return -5.0 < eta && eta < -2.0;
    }, "fn Cut" );
  correction_task.AddVector(fn);
  
  std::vector<Qn::AxisD> proton_axes{
          { "y", 20, -1.0, 1.0 },
          { "pT", 20, 0.0, 2.0 },
  };

  std::vector<Qn::AxisD> pion_axes{
          { "y", 20, -1.0, 1.0 },
          { "pT", 20, 0.0, 1.0 },
  };

  VectorConfig proton( "proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  proton.AddCut( "eta_cm", [ETA_MIN, ETA_MAX](double eta){
    return ETA_MIN < eta && eta < ETA_MAX;
    }, "acceptance cut" );
  correction_task.AddVector(proton);

  VectorConfig Tp( "Tp", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1} );
  Tp.SetCorrections( {CORRECTION::PLAIN } );
  Tp.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tp.AddCut( "eta_cm", [](double eta){
    return 0.4 < eta && eta < 1.0;
    }, "Tp ycm cut" );
  Tp.AddCut( "pT", [](double pT){
    return 0.4 < pT && pT < 2.0;
    }, "Tp pT cut" );
  correction_task.AddVector(Tp);

  VectorConfig Tn( "Tn", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tn.SetHarmonicArray( {1} );
  Tn.SetCorrections( {CORRECTION::PLAIN } );
  Tn.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tn.AddCut( "eta_cm", [](double eta){
    return -1.0 < eta && eta < -0.4;
    }, "Tn ycm cut" );
  Tn.AddCut( "pT", [](double pT){
    return 0.4 < pT && pT < 2.0;
    }, "Tn pT cut" );
  correction_task.AddVector(Tn);

  VectorConfig Tp2( "Tp2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp2.SetHarmonicArray( {2} );
  Tp2.SetCorrections( {CORRECTION::PLAIN } );
  Tp2.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tp2.AddCut( "y", [](double y){
    return 0.2 < y && y < 0.5;
    }, "Tp2 ycm cut" );
  Tp2.AddCut( "pT", [](double pT){
    return 0.4 < pT && pT < 2.0;
    }, "Tp2 pT cut" );
  correction_task.AddVector(Tp2);

  VectorConfig Tn2( "Tn2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tn2.SetHarmonicArray( { 2} );
  Tn2.SetCorrections( {CORRECTION::PLAIN } );
  Tn2.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tn2.AddCut( "y", [](double y){
    return -0.5 < y && y < -0.2;
    }, "Tn2 ycm cut" );
  Tn2.AddCut( "pT", [](double pT){
    return 0.4 < pT && pT < 2.0;
    }, "Tn2 pT cut" );
  correction_task.AddVector(Tn2);

  VectorConfig psi_rp( "psi_rp", "reaction_plane", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}