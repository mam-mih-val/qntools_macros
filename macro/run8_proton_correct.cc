//
// Created by Misha on 3/7/2023.
//

void run8_proton_correct(std::string list){
  std::vector<int> f1_modules = {11, 12, 13, 16, 17, 20, 21, 22};
  std::vector<int> f2_modules = {5, 6, 7, 8, 9, 10, 14, 15, 18, 19, 23, 24, 25, 26, 27, 28};
  std::vector<int> f3_modules = {0, 1, 2, 3, 4, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53};
  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  auto dd=d
          .Define("centrality",
                  "float centrality;"
                  "std::vector<float> centrality_percentage{ 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100 };\n"
                  "std::vector<int> multiplicity_edges{ 269, 170, 143, 121, 103, 86, 72, 60, 49, 33, 21, 12, 7, 4, 0 };\n"
                  "auto multiplicity = trMom.size();\n"
                  "int idx = 0;\n"
                  "float bin_edge = multiplicity_edges[idx];\n"
                  "while( multiplicity < bin_edge &&\n"
                    "idx < multiplicity_edges.size()-1 ){\n"
                    "idx++;\n"
                    "bin_edge = multiplicity_edges[idx];\n"
                  "}\n"
                  "centrality = (centrality_percentage[idx-1] + centrality_percentage[idx])/2.0f;\n"
                  "return centrality;")
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
          .Filter("1e4 < bc1Integral.at(0) && bc1Integral.at(0) < 4e4" )
          .Filter("vtxChi2/vtxNdf > 0.1")
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("centrality"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex|Y|Pid|IsProton|MotherId|Charge)"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 8, 0, 40} );

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

  std::vector<Qn::AxisD> negative_axes{
          { "trEta", 5, 0.5, 5.5 },
          { "trPt", 5, 0.0, 1.0 },
  };

  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tneg.AddCut( "trCharge", [](double charge){
    return charge < 0.0;
    }, "charge" );
  Tneg.AddCut( "trEta", [](double eta){
    return 1.5 < eta && eta < 4.0;
    }, "eta cut" );
  Tneg.AddCut( "trPt", [](double pT){
    return pT > 0.2;
    }, "pT cut" );
  correction_task.AddVector(Tneg);

  VectorConfig Tpos( "Tpos", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tpos.SetHarmonicArray( {1, 2} );
  Tpos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  Tpos.SetCorrectionAxes( negative_axes );
  Tpos.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  correction_task.AddVector(Tpos);

  correction_task.Run();
}