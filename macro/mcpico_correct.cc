
//
// Created by Misha on 3/7/2023.
//

void mcpico_correct(std::string list, std::string str_sqrt_snn="2.4", std::string str_nucleus_mass="197"){
  const double ETA_MIN = 0.5;
  const double ETA_MAX = 2.5;

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
          .Define("y_norm",[Y_BEAM]( ROOT::VecOps::RVec<float> vec_y ){
            ROOT::VecOps::RVec<float> y_norm;
            for( int i=0; i<vec_y.size(); ++i ){
              auto y = vec_y.at(i);
              y_norm.push_back( y/Y_BEAM );
            }
            return y_norm;
          },{ "y" })
          .Define( "is_accepted", [ETA_MIN, ETA_MAX](ROOT::VecOps::RVec<float> vec_eta, ROOT::VecOps::RVec<float> vec_pT){
            ROOT::VecOps::RVec<int> is_accepted;
            for( int i=0; i<vec_eta.size(); ++i ){
              auto eta = vec_eta.at(i);
              auto pT = vec_pT.at(i);
              if( eta < ETA_MIN || eta > ETA_MAX ){
                is_accepted.push_back( 0 );
                continue;
              }
              if( pT < 0.3 ){
                is_accepted.push_back( 0 );
                continue;
              }
              is_accepted.push_back(1);
            }
            return is_accepted;
          }, { "eta_lab", "pT" } )
          .Define( "rnd_weights", [](
                  ROOT::VecOps::RVec<int> is_accepted,
                  ROOT::VecOps::RVec<int> vec_pid,
                  ROOT::VecOps::RVec<float> vec_y
                  ){
            ROOT::VecOps::RVec<float> vec_weights;
            for( int i=0; i<is_accepted.size(); ++i ){
              auto acceptance = is_accepted.at(i);
              auto pid = vec_pid.at(i);
              auto y = vec_y.at(i);
              if( acceptance == 0 ){
                vec_weights.push_back(0);
                continue;
              }
              if( pid != 2212 ){
                vec_weights.push_back(0);
                continue;
              }
              auto weight = fabs(y) < 0.8 ? y / 0.8 : y/fabs(y);
              vec_weights.push_back(weight);
            }
            return vec_weights;
          }, {"is_accepted", "pdg", "y_norm"} )
          .Define("dphi",[](
                  ROOT::VecOps::RVec<float> vec_phi,
                  ROOT::VecOps::RVec<float> rnd_weights
                  ){
            ROOT::VecOps::RVec<float> vec_dphi;
            float sum_wx{};
            float sum_wy{};
            // Filling the global Q-vector
            for(int i=0; i<vec_phi.size(); ++i){
              auto phi = vec_phi.at(i);
              auto weight = rnd_weights.at(i);
              sum_wx += weight * cos(phi);
              sum_wy += weight * sin(phi);
            }
            // Calculating the dphi for each particle
            for(int i=0; i<vec_phi.size(); ++i){
              auto phi = vec_phi.at(i);
              auto weight = rnd_weights.at(i);
              auto Qx = sum_wx - weight * cos(phi);
              auto Qy = sum_wy - weight * sin(phi);
              auto psi = atan2( Qy, Qx );
              vec_dphi.push_back( phi - psi );
            }
            return vec_dphi;
          },{ "phi", "rnd_weights" })
          .Define("psi12",[ETA_MIN, ETA_MAX](
                  ROOT::VecOps::RVec<float> vec_phi,
                  ROOT::VecOps::RVec<float> rnd_weights
                  ){
            std::random_device rnd_device;
            std::mt19937 generator(rnd_device());
            std::uniform_int_distribution<int> distribution(0,1);
            std::array<float, 2> sum_wx{};
            std::array<float, 2> sum_wy{};
            for(int i=0; i<vec_phi.size(); ++i){
              auto phi = vec_phi.at(i);
              auto weight = rnd_weights.at(i);
              auto idx = distribution(generator);
              sum_wx.at(idx) += weight * cos(phi);
              sum_wy.at(idx) += weight * sin(phi);
            }
            auto psi1 = atan2( sum_wy.at(0), sum_wx.at(0) );
            auto psi2 = atan2( sum_wy.at(1), sum_wx.at(1) );
            return psi1 - psi2;
          },{ "phi", "rnd_weights" })
          .Filter("for(auto x : phi){ if( std::isnan(x) ) return false; } return true;")
          .Filter("for(auto x : y){ if( std::isnan(x) ) return false; } return true;")
          ;

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables( std::regex("bimp|b_norm|reaction_plane|psi12") );
  correction_task.SetTrackVariables({
    std::regex("phi|pT|y|eta_lab|pdg|charge|dphi|is_accepted")
  });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"b_norm", 20, 0, 2} );

  VectorConfig f1( "F1", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN } );
  f1.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2112 || pdg_code == 2212;
  }, "proton cut" );
  f1.AddCut( "eta_lab", [](double eta){
    return 4.4 < eta && eta < 5.5;
    }, "F1 Cut" );
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN } );
  f2.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2112 || pdg_code == 2212;
  }, "proton cut" );
  f2.AddCut( "eta_lab", [](double eta){
    return 3.9 < eta && eta < 4.4;
    }, "F2 Cut" );
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN } );
  f3.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2112 || pdg_code == 2212;
  }, "proton cut" );
  f3.AddCut( "eta_lab", [](double eta){
    return 3.1 < eta && eta < 3.9;
    }, "F3 Cut" );
  correction_task.AddVector(f3);


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
  proton.AddCut( "eta_lab", [ETA_MIN, ETA_MAX](double eta){
    return ETA_MIN < eta && eta < ETA_MAX;
    }, "acceptance cut" );
  correction_task.AddVector(proton);

  VectorConfig pion( "pi_neg", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  pion.SetHarmonicArray( {1, 2} );
  pion.SetCorrections( {CORRECTION::PLAIN } );
  pion.SetCorrectionAxes( pion_axes );
  pion.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == -211;
    }, "pion cut" );
  pion.AddCut( "eta_lab", [ETA_MIN, ETA_MAX](double eta){
    return ETA_MIN < eta && eta < ETA_MAX;
  }, "acceptance cut" );
  correction_task.AddVector(pion);

  VectorConfig Tp( "Tp", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN } );
  Tp.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  Tp.AddCut( "y", [](double ycm){
    return 0.4 < ycm && ycm < 0.6;
    }, "Tp ycm cut" );
  Tp.AddCut( "pT", [](double pT){
    return 0.4 < pT && pT < 2.0;
    }, "Tp pT cut" );
  correction_task.AddVector(Tp);

  VectorConfig Tneg( "Tpi", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "pdg", [](double double_pid){
    auto pid = static_cast<int>(double_pid);
    return pid == -211 || pid == 211;
    }, "pion cut" );
  Tneg.AddCut( "y", [](double y){
    return 0.2 < y && y < 0.8;
    }, "Tneg y cut" );
  Tneg.AddCut( "pT", [](double pT){
    return 0.1 < pT && pT < 0.5;
    }, "Tneg pT cut" );
  correction_task.AddVector(Tneg);

  VectorConfig Tneg2( "Tpi2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg2.SetHarmonicArray( {1, 2} );
  Tneg2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg2.AddCut( "pdg", [](double double_pid){
    auto pid = static_cast<int>(double_pid);
    return pid == -211 || pid == 211;
    }, "pion cut" );
  Tneg2.AddCut( "is_accepted", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 1;
  }, "acceptance cut" );
  correction_task.AddVector(Tneg2);

  VectorConfig rs_proton( "rnd_proton", "dphi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  rs_proton.SetHarmonicArray( {1, 2} );
  rs_proton.SetCorrections( {CORRECTION::PLAIN } );
  rs_proton.SetCorrectionAxes( proton_axes );
  rs_proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
  }, "proton cut" );
  rs_proton.AddCut( "eta_lab", [ETA_MIN, ETA_MAX](double eta){
    return ETA_MIN < eta && eta < ETA_MAX;
  }, "acceptance cut" );
  correction_task.AddVector(rs_proton);

  VectorConfig rnd_sub( "rnd_sub", "psi12", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  rnd_sub.SetHarmonicArray( {1, 2} );
  rnd_sub.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(rnd_sub);

  VectorConfig psi_rp( "psi_rp", "reaction_plane", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}