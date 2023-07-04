
//
// Created by Misha on 3/7/2023.
//

void mcpico_correct(std::string list, std::string str_sqrt_snn="2.4"){
  const double sqrt_snn = std::stod(str_sqrt_snn);
  const double M = 0.938;
  const double T = sqrt_snn * sqrt_snn/ 2 / M - 2*M;
  const double GAMMA = (T + M) / M;
  const double BETA = sqrt(1 - (M * M) / (M + T) / (M + T));
  const double PZ = M * BETA * GAMMA;
  const double E = T + M;
  const double Y_BEAM = 0.5 * log((E + PZ) / (E - PZ)) / 2.0;
  TStopwatch timer;
  timer.Start();
  std::string treename = "mctree";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
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
          .Define("eta_lab",[Y_BEAM]( ROOT::VecOps::RVec<float> vec_y, ROOT::VecOps::RVec<float> vec_pT, ROOT::VecOps::RVec<int> vec_pid ){
            ROOT::VecOps::RVec<float> vec_eta;
            for( int i=0; i<vec_y.size(); ++i ){
              auto y = vec_y.at(i);
              auto pT = vec_pT.at(i);
              auto pid = vec_pid.at(i);
              auto mass = TDatabasePDG::Instance()->GetParticle(pid)->Mass();
              auto mT = sqrt( pT*pT + mass*mass );
              auto eta_lab = asinh( mT/pT * sinh( y + Y_BEAM ) );
              vec_eta.push_back(eta_lab);
            }
            return vec_eta;
          },{ "y", "pT", "pdg", })
          .Filter("bimp<14.0")
          ;

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables( std::regex("bimp|reaction_plane") );
  correction_task.SetTrackVariables({
    std::regex("phi|pT|y|eta_lab|pdg|rnd_sub|charge")
  });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"bimp", 14, 0, 14} );

  VectorConfig f1( "F1", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN } );
//  f1.AddCut( "pdg", [](double pid){
//    auto pdg_code = static_cast<int>(pid);
//    return pdg_code == 2112;
//  }, "proton cut" );
  f1.AddCut( "eta_lab", [](double eta){
    return 3.8 < eta && eta < 5.4;
    }, "F1 Cut" );
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN } );
//  f2.AddCut( "pdg", [](double pid){
//    auto pdg_code = static_cast<int>(pid);
//    return pdg_code == 2112;
//  }, "proton cut" );
  f2.AddCut( "eta_lab", [](double eta){
    return 3.3 < eta && eta < 3.8;
    }, "F2 Cut" );
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN } );
//  f3.AddCut( "pdg", [](double pid){
//    auto pdg_code = static_cast<int>(pid);
//    return pdg_code == 2112;
//  }, "proton cut" );
  f3.AddCut( "eta_lab", [](double eta){
    return 2.7 < eta && eta < 3.3;
    }, "F3 Cut" );
  correction_task.AddVector(f3);

  VectorConfig psi_rp( "psi_rp", "reaction_plane", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}