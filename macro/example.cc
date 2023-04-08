//
// Created by Misha on 3/7/2023.
//

void example(){
  std::vector<int> f1_modules = {11, 12, 13, 16, 17, 20, 21, 22};
  std::vector<int> f2_modules = {5, 6, 7, 8, 9, 10, 14, 15, 18, 19, 23, 24, 25, 26, 27, 28};
  std::vector<int> f3_modules = {0, 1, 2, 3, 4, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  ROOT::RDataFrame d("t", "/home/mikhail/bmn_plain_tree/330.tree.root");
  auto dd=d
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("scwallModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
          .Define("trY","ROOT::VecOps::RVec<float> ycm;"
                        "for(int i=0; i<trMom.size(); ++i){ "
                          "auto matching_id = trSimIndex.at(i);"
                          "auto m = simMom.at(i).mass();"
                          "auto pz = trMom.at(i).pz();"
                          "auto p = trMom.at(i).P();"
                          "auto E = sqrt(m*m + p*p);"
                          "auto y = 0.5*log( (E+pz)/(E-pz) );"
                          "ycm.push_back(y);"
                        "}"
                        " return ycm;")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom:trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:trMom) phi.push_back(mom.phi()); return phi;")
          .Define("simPt","ROOT::VecOps::RVec<float> pt;for(auto& mom:simMom) pt.push_back(mom.pt()); return pt;")
          .Define("simEta","ROOT::VecOps::RVec<float> eta;for(auto& mom:simMom) eta.push_back(mom.eta()); return eta;")
          .Define("simPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:simMom) phi.push_back(mom.phi()); return phi;")
          .Define("simY","ROOT::VecOps::RVec<float> y;for(auto& mom:simMom) y.push_back(mom.Rapidity()); return y;")
          .Filter("vtxChi2>0.0001"); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "correction_in.root" );
  correction_task.SetEventVariables(std::regex("b"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex|Y)"),
                                            std::regex("sim(Pt|Eta|Phi|Pdg|MotherId|Y)")
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"b", 4,0,10} );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  f1.AddCut( "fhcalModId", [f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "Bruh" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  //
//  std::vector<Qn::AxisD> corrAxesParticle=
//          {
//                  {"trPt",4,0,2},
//                  {"trEta",4,0,2},
//          };
//
//  std::vector<Qn::AxisD> corrAxesSim=
//          {
//                  {"simPt",4,0,2},
//                  {"simEta",4,0,2},
//          };
//
//
//  for (auto &axis:corrAxesEvent)
//    man->AddCorrectionAxis(axis);
//
//  Qn::Recentering recentering;
//  recentering.SetApplyWidthEqualization(false);
//  Qn::TwistAndRescale twistRescale;
//  twistRescale.SetApplyRescale(true);
//  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
//
//  auto sumW=Qn::QVector::Normalization::M;
//  auto track=Qn::DetectorType::TRACK;
//  auto channel=Qn::DetectorType::CHANNEL;
//  auto plain=Qn::QVector::CorrectionStep::PLAIN;
//  auto recentered=Qn::QVector::CorrectionStep::RECENTERED;
//  auto twisted=Qn::QVector::CorrectionStep::TWIST;
//  auto rescaled=Qn::QVector::CorrectionStep::RESCALED;
//
//  man->AddDetector("fhcal1", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
//  man->AddCorrectionOnQnVector("fhcal1", recentering);
//  man->AddCorrectionOnQnVector("fhcal1", twistRescale);
//  man->SetOutputQVectors("fhcal1", {plain, recentered});
//  man->AddCutOnDetector("fhcal1", {"fhcalModInSub1"}, [](double val){ return static_cast<bool>(val == 1); }, "fhcal1");
//  man->AddHisto2D("fhcal1", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}});
//  man->AddHisto2D("fhcal1", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
//
//  man->AddDetector("sim", track, "simPhi", "Ones", corrAxesSim, {1,2}, sumW);
//  man->AddCutOnDetector("sim", {"track_type"}, [](double val){ return static_cast<bool>(fabs(val - 1)<0.1); }, "sim");
////  correction_manager_.AddCutOnDetector("tr", {"simPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
////  man->AddCorrectionOnQnVector("sim", recentering);
////  man->AddCorrectionOnQnVector("sim", twistRescale);
//  man->SetOutputQVectors("sim", {plain});
//  man->AddHisto1D("sim", {"simPhi", 100, -3.15, 3.15}, "Ones");
//  man->AddHisto2D("sim", {{"simY", 100, 0.0, 4.0},
//                         {"simPt", 100, 0.0, 2.0}}, "Ones");
////  correction_manager_.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");

  correction_task.Run();
}