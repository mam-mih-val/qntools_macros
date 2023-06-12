
//
// Created by Misha on 3/7/2023.
//

void mcpico_correct(std::string list){
  const float Y_BEAM = 0.731032;
  TStopwatch timer;
  timer.Start();
  std::string treename = "mctree";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define("phi", []( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py ){
            ROOT::VecOps::RVec<float> vec_phi;
            for( int i=0; i<vec_px.size(); ++i ){
              auto px = vec_px.at(i);
              auto py = vec_py.at(i);
              vec_phi.push_back( atan2(py, px) );
            }
            return vec_phi;
          },{ "momx", "momy", })
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
          .Define("rnd_sub",[]( ROOT::VecOps::RVec<float> vec_y ){
            ROOT::VecOps::RVec<int> rnd_sub;
            std::random_device rnd_device;
            std::mt19937 generator(rnd_device());
            std::uniform_int_distribution<int> distribution(0,1); // distribution in range [0, 1]
            for( int i=0; i<vec_y.size(); ++i ){
              rnd_sub.push_back( distribution(generator) );
            }
            return rnd_sub;
          },{ "y" })
          .Define("y_norm",[Y_BEAM]( ROOT::VecOps::RVec<float> vec_y ){
            ROOT::VecOps::RVec<float> y_norm;
            for( int i=0; i<vec_y.size(); ++i ){
              auto y = vec_y.at(i);
              y_norm.push_back( y/Y_BEAM );
            }
            return y_norm;
          },{ "y" })
          .Define("dphi",[](
                  ROOT::VecOps::RVec<float> vec_phi,
                  ROOT::VecOps::RVec<float> vec_y,
                  ROOT::VecOps::RVec<int> vec_pid,
                  ROOT::VecOps::RVec<float> vec_eta
                  ){
            ROOT::VecOps::RVec<float> vec_dphi;
            for(int i=0; i<vec_phi.size(); ++i){
              auto phi = vec_phi.at(i);
              float sum_wx{};
              float sum_wy{};
              for( int j=0; j<vec_phi.size(); ++j ){
                if( i == j )
                  continue;
                auto pid = vec_pid.at(j);
                if( pid != 2212 )
                  continue;
                auto eta = vec_eta.at(j);
                if( eta < 0 || eta > 2 )
                  continue;
                auto y = vec_y.at(j);
                auto phi2 = vec_phi.at(j);
                auto weight = y < 0.8 ? y / 0.8 : 1.0;
                sum_wx += weight * cos(phi2);
                sum_wy += weight * sin(phi2);
              }
              auto psi = atan2( sum_wy, sum_wx );
              vec_dphi.push_back( phi - psi );
            }
            return vec_dphi;
          },{ "phi", "y_norm", "pdg", "eta_lab" })
          .Define("psi12",[](
                  ROOT::VecOps::RVec<float> vec_phi,
                  ROOT::VecOps::RVec<float> vec_y,
                  ROOT::VecOps::RVec<int> vec_pid,
                  ROOT::VecOps::RVec<int> rnd_sub,
                  ROOT::VecOps::RVec<int> vec_eta
                  ){
            float vec_dphi;
            // First RS
            std::array<float, 2> sum_wx{};
            std::array<float, 2> sum_wy{};
            for(int i=0; i<vec_phi.size(); ++i){
              auto pid = vec_pid.at(i);
              if( pid != 2212 )
                continue;
              auto eta = vec_eta.at(i);
              if( eta < 0 || eta > 2 )
                continue;
              auto phi = vec_phi.at(i);
              auto y = vec_y.at(i);
              auto weight = y < 0.8 ? y / 0.8 : 1.0;
              auto idx = rnd_sub.at(i);
              sum_wx.at(idx) += weight * cos(phi);
              sum_wy.at(idx) += weight * sin(phi);
              for( int j=0; j<vec_phi.size(); ++j ){
                if( i == j )
                  continue;
              }
            }
            auto psi1 = atan2( sum_wy.at(0), sum_wx.at(0) );
            auto psi2 = atan2( sum_wy.at(1), sum_wx.at(1) );
            return psi1 - psi2;
          },{ "phi", "y_norm", "pdg", "rnd_sub", "eta_lab" })
          .Filter("bimp<14.0")

          ;

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables( std::regex("bimp|phi2|psi12") );
  correction_task.SetTrackVariables({
    std::regex("phi|pT|y|eta_lab|pdg|rnd_sub|charge|dphi")
  });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"bimp", 14, 0, 14} );

  VectorConfig f1( "F1", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN } );
  f1.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212 || pdg_code == 2112;
  }, "proton cut" );f1.AddCut( "eta_lab", [](double eta){
    return 3.8 < eta && eta < 5.4;
    }, "F1 Cut" );
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN } );
  f2.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212 || pdg_code == 2112;
  }, "proton cut" );f2.AddCut( "eta_lab", [](double eta){
    return 3.3 < eta && eta < 3.8;
    }, "F2 Cut" );
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN } );
  f3.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212 || pdg_code == 2112;
  }, "proton cut" );
  f3.AddCut( "eta_lab", [](double eta){
    return 2.7 < eta && eta < 3.3;
    }, "F3 Cut" );
  correction_task.AddVector(f3);


  std::vector<Qn::AxisD> proton_axes{
          { "y", 20, -1.0, 1.0 },
          { "pT", 20, 0.0, 2.0 },
  };

  VectorConfig proton( "proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "eta_lab", [](double eta){
    return 0. < eta && eta < 2.;
  }, "acceptance cut" );
  proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  correction_task.AddVector(proton);

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
  Tneg.AddCut( "pdg", [](double pid){
    return static_cast<int>(pid) == -211;
    }, "negative pion cut" );
  Tneg.AddCut( "y", [](double y){
    return 0.2 < y && y < 0.8;
    }, "Tneg y cut" );
  Tneg.AddCut( "pT", [](double pT){
    return 0.1 < pT && pT < 0.5;
    }, "Tneg pT cut" );
  correction_task.AddVector(Tneg);

  VectorConfig rs_proton( "rnd_proton", "dphi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  rs_proton.SetHarmonicArray( {1, 2} );
  rs_proton.SetCorrections( {CORRECTION::PLAIN } );
  rs_proton.SetCorrectionAxes( proton_axes );
  rs_proton.AddCut( "eta_lab", [](double eta){
    return 0. < eta && eta < 2.;
  }, "acceptance cut" );
  rs_proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
  }, "proton cut" );
  correction_task.AddVector(rs_proton);

  VectorConfig rnd_sub( "rnd_sub", "psi12", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  rnd_sub.SetHarmonicArray( {1, 2} );
  rnd_sub.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(rnd_sub);

  VectorConfig psi_rp( "psi_rp", "phi2", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  correction_task.Run();
}