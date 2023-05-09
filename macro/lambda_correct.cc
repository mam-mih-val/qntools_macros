//
// Created by Misha on 3/7/2023.
//

void lambda_correct(std::string list, std::string str_efficiency_file){
  std::vector<int> f1_modules = {11, 12, 13, 16, 17, 20, 21, 22};
  std::vector<int> f2_modules = {5, 6, 7, 8, 9, 10, 14, 15, 18, 19, 23, 24, 25, 26, 27, 28};
  std::vector<int> f3_modules = {0, 1, 2, 3, 4, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
  TFilePtr efficiency_file{ str_efficiency_file };
  TH2D* efficiency_histo{nullptr};
  efficiency_file->GetObject( "h2_efficiency", efficiency_histo );
  if( !efficiency_histo ){ std::cout << "Efficiency histogram cannot be retrieved from file" << std::endl; }
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("scwallModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
          .Define("simPt","ROOT::VecOps::RVec<float> pt;for(auto& mom:simMom) pt.push_back(mom.pt()); return pt;")
          .Define("simEta","ROOT::VecOps::RVec<float> eta;for(auto& mom:simMom) eta.push_back(mom.eta()); return eta;")
          .Define("simPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:simMom) phi.push_back(mom.phi()); return phi;")
          .Define("simY","ROOT::VecOps::RVec<float> y;for(auto& mom:simMom) y.push_back( mom.Rapidity() - 1.0 ); return y;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
          .Define("trY",
                  "ROOT::VecOps::RVec<float> ycm{};\n"
                  "for(int i=0; i<trMom.size(); ++i){\n"
                  " auto matching_id = trSimIndex.at(i);\n"
                  " if( matching_id < 0 || matching_id > simMom.size() ){\n"
                  "   ycm.push_back(-999.);\n"
                  "   continue;\n"
                  " }\n"
                  " auto m = simMom.at(matching_id).mass();\n"
                  " auto pz = trMom.at(i).pz();\n"
                  " auto p = trMom.at(i).P();\n"
                  " auto E = sqrt(m*m + p*p);\n"
                  " auto y = 0.5*log( (E + pz) / (E - pz) );\n"
                  " ycm.push_back( y - 1.0 );\n"
                  "}\n"
                  "return ycm;\n")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom:trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:trMom) phi.push_back(mom.phi()); return phi;")
          .Define("trPid",
                  "ROOT::VecOps::RVec<int> pid{};\n"
                  "for( auto& id : trSimIndex){\n"
                  " if(id > simPdg.size()){\n"
                  "   pid.push_back(-1);\n"
                  "   continue;\n"
                  " }\n"
                  " if(id < 0){\n"
                  "   pid.push_back(-1);\n"
                  "   continue;\n"
                  " }\n"
                  " pid.push_back(simPdg.at(id));\n"
                  "}\n"
                  "return pid;\n")
          .Define("trMotherId",
                  "ROOT::VecOps::RVec<int> mother_id{};\n"
                  "for( auto& id : trSimIndex){\n"
                  " if(id > simMotherId.size()){\n"
                  "   mother_id.push_back(-999);\n"
                  "   continue;\n"
                  " }\n"
                  " if(id < 0){\n"
                  "   mother_id.push_back(-999);\n"
                  "   continue;\n"
                  "  }\n"
                  " mother_id.push_back(simMotherId.at(id));\n"
                  "}\n"
                  "return mother_id;\n")
          .Define( "trIsProton", "trPid == 2212" )
          .Define( "trIsPiNeg", "trPid == -211" )
          .Define( "trIsPiPos", "trPid == 211" )
          .Define("candidate_pT", "std::vector<float> pT; for( auto mom : candidate_momenta ){ pT.push_back( mom.Pt() ); } return pT;")
          .Define("candidate_phi", "std::vector<float> phi; for( auto mom : candidate_momenta ){ phi.push_back( mom.Phi() ); } return phi;")
          .Define("candidate_rapidity", "std::vector<float> rapidity; for( auto mom : candidate_momenta ){ rapidity.push_back( mom.Rapidity() - 1.0 ); } return rapidity;")
          .Define("candidate_weight",
                  [efficiency_histo]( ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >> momenta ){
            std::vector<float> weights;
            for( auto mom : momenta ){
              if( !efficiency_histo ){
                weights.push_back(1.f);
              }
              auto y = mom.Rapidity();
              auto pT = mom.Pt();
              auto y_bin = efficiency_histo->GetXaxis()->FindBin( y );
              auto pT_bin = efficiency_histo->GetYaxis()->FindBin( pT );
              auto efficiency = efficiency_histo->GetBinContent( y_bin, pT_bin );
              auto weight = efficiency > 0.0001 ? 1.0/efficiency : 0.0;
              weights.push_back( weight );
            }
            return weights;
          },
                  {"candidate_momenta"})
          .Define("m_err", "std::vector<float> err; for( auto mom : candidate_momentum_errors ){ err.push_back( mom.at(3) ); } return err;")
          .Define("daughter1_cos", "std::vector<float> cosine; for( int i=0; i<daughter_cosines.at(0).size(); ++i ){ cosine.push_back( daughter_cosines.at(0).at(i) ); } return cosine;")
          .Define("daughter2_cos", "std::vector<float> cosine; for( int i=0; i<daughter_cosines.at(1).size(); ++i ){ cosine.push_back( daughter_cosines.at(1).at(i) ); } return cosine;")
          .Define("daughter1_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(0).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(0).at(i) ); } return chi2;")
          .Define("daughter2_chi2_prim", "std::vector<float> chi2; for( int i=0; i<daughter_chi2_prim.at(1).size(); ++i ){ chi2.push_back( daughter_chi2_prim.at(1).at(i) ); } return chi2;")
          .Define("candidate_signal", "candidate_true_pid == 3122")
          .Define("candidate_background", "candidate_true_pid != 3122")
          .Define("candidate_good",
                  "std::vector<int> good_candidate;\n"
                  "for(int i=0; i<daughter_cosines.at(0).size(); ++i){\n"
                  " if( m_err.at(i) > 0.001 ){\n"
                  "  good_candidate.push_back(0);\n"
                  "  continue;\n"
                  " }\n"
                  " if( candidate_chi2_topo.at(i) > 51. ){\n"
                  "  good_candidate.push_back(0);\n"
                  "  continue;\n"
                  " }\n"
                  " good_candidate.push_back( 1 );\n"
                  "}\n"
                  "return good_candidate;\n"
          )
          .Define("candidate_good_signal",
                  "std::vector<int> selected_signal;"
                  "for(int i=0; i<candidate_good.size(); ++i ){\n"
                  "  if( candidate_signal.at(i) == 0 ){\n"
                  "    selected_signal.push_back(0);\n"
                  "    continue;\n"
                  "  }\n"
                  "  if( candidate_good.at(i) == 0 ){\n"
                  "    selected_signal.push_back(0);\n"
                  "    continue;\n"
                  "  }\n"
                  "  selected_signal.push_back( 1 );\n"
                  "}\n"
                  "return selected_signal;")
          .Define("candidate_good_background",
                  "std::vector<int> selected_background;\n"
                  "for(int i=0; i<candidate_good.size(); ++i ){\n"
                  " if( candidate_signal.at(i) == 1 ){\n"
                  "   selected_background.push_back(0);\n"
                  "   continue;\n"
                  " }\n"
                  " if( candidate_good.at(i) == 0 ){\n"
                  "   selected_background.push_back(0);\n"
                  "   continue;\n"
                  " }\n"
                  " selected_background.push_back( 1 );\n"
                  "}\n"
                  "return selected_background;")
          .Filter("vtxChi2>0.0001"); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("centrality|psiRP"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex|Y|Pid|IsProton|MotherId|Charge)"),
                                            std::regex("candidate_(pT|rapidity|phi|weight|mass|signal|background|good|good_signal|good_background)"),
                                            std::regex("sim(Pt|Eta|Phi|Pdg|MotherId|Y)")
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 1, 10, 40} );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f1.AddCut( "fhcalModId", [f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f2.AddCut( "fhcalModId", [f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f3.AddCut( "fhcalModId", [f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  VectorConfig Tp( "Tp", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tp.SetHarmonicArray( {1, 2} );
  Tp.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tp.AddCut( "trPid", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
  }, "proton cut" );
  Tp.AddCut( "trY", [](double ycm){
    return 0.4 < ycm && ycm < 0.6;
  }, "Tp ycm cut" );
  Tp.AddCut( "trPt", [](double pT){
    return 0.4 < pT && pT < 2.0;
  }, "Tp pT cut" );
  Tp.AddCut( "trMotherId", [](double mother_id){
    auto int_mother_id = static_cast<int>(mother_id);
    return int_mother_id == -1;
  }, "cut on primary" );
  correction_task.AddVector(Tp);

  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "trCharge", [](double charge){
    return charge < 0.0;
  }, "charge" );
  Tneg.AddCut( "trEta", [](double eta){
    return 1.0 < eta && eta < 2.0;
  }, "Tneg eta cut" );
  Tneg.AddCut( "trPt", [](double pT){
    return 0.1 < pT && pT < 0.5;
  }, "Tneg pT cut" );
  correction_task.AddVector(Tneg);

  VectorConfig psi_rp( "psi_rp", "psiRP", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  psi_rp.AddHisto1D({"psi_rp", 100, -3.5, 3.5}, "psiRP");
  correction_task.AddVector(psi_rp);

  std::vector<Qn::AxisD> rec_lamda_axes{
          { "candidate_rapidity", 6, -0.2, 1.0 },
          { "candidate_pT", 7, 0.0, 1.4 },
          { "candidate_mass", 10, 1.10, 1.13 },
  };
  std::vector<Qn::AxisD> tru_lamda_axes{
          { "simY", 12, -0.2, 1.0 },
          { "simPt", 14, 0.0, 1.4 },
  };

  VectorConfig lambda_signal( "lambda_signal", "candidate_phi", "candidate_weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_signal.SetHarmonicArray( {1, 2} );
  lambda_signal.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, } );
  lambda_signal.SetCorrectionAxes( rec_lamda_axes );
  lambda_signal.AddCut( "candidate_good_signal", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is signal" );
  correction_task.AddVector(lambda_signal);

  VectorConfig lambda_background( "lambda_background", "candidate_phi", "candidate_weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_background.SetHarmonicArray( {1, 2} );
  lambda_background.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, } );
  lambda_background.SetCorrectionAxes( rec_lamda_axes );
  lambda_background.AddCut( "candidate_good_background", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is signal" );
  correction_task.AddVector(lambda_background);

  VectorConfig lambda_good( "lambda_good", "candidate_phi", "candidate_weight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_good.SetHarmonicArray( {1, 2} );
  lambda_good.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, } );
  lambda_good.SetCorrectionAxes( rec_lamda_axes );
  lambda_good.AddCut( "candidate_good", [](double is_signal){
    auto int_is_signal = static_cast<int>(is_signal);
    return int_is_signal == 1;
  }, "cut on if is good candidate" );
  correction_task.AddVector(lambda_good);

  VectorConfig lambda_tru( "lambda_true", "simPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  lambda_tru.SetHarmonicArray( {1, 2} );
  lambda_tru.SetCorrections( {CORRECTION::PLAIN } );
  lambda_tru.SetCorrectionAxes( tru_lamda_axes );
  lambda_tru.AddCut( "simPdg", [](double pid_code){
    auto int_pid_code = static_cast<int>(pid_code);
    return int_pid_code == 3122;
  }, "cut on pid" );
  correction_task.AddVector(lambda_tru);


  correction_task.Run();
}