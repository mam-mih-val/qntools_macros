//
// Created by Misha on 3/7/2023.
//

#include <cmath>
void fhcal65_correct(std::string list){
  std::vector<int> f1_modules = {
    22, 23, 24, 
    32,     33, 
    41, 42, 43 
  };
  std::vector<int> f2_modules = {
    11, 12, 13, 14, 15, 16,
    21,             25, 26,
    31,             34, 35, 
    40,             44, 45, 
    50, 51, 52, 53, 54, 55, 
  };
  std::vector<int> f3_modules = {
         1,  2,  3,  4,  5,  6,  7,  8,  
     9, 10,                         17, 18,
    19, 20,                         27, 28,
    29, 30,                         36, 37,
    38, 39,                         46, 47,
    48, 49,                         56, 57,
        58, 59, 60, 61, 62, 63, 64, 65,    
  };
  const double R0 = 1.25; // fm
  const double XE_A = 131.0;
  const double XE_RADIUS = R0 * std::pow( XE_A, 1.0 / 3.0 );
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
          .Define( "bNorm", [XE_RADIUS]( float b ){ return b / XE_RADIUS; }, {"b"}  )
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("scwallModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
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
          .Define("simPt","ROOT::VecOps::RVec<float> pt;for(auto& mom:simMom) pt.push_back(mom.pt()); return pt;")
          .Define("simEta","ROOT::VecOps::RVec<float> eta;for(auto& mom:simMom) eta.push_back(mom.eta()); return eta;")
          .Define("simPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:simMom) phi.push_back(mom.phi()); return phi;")
          .Define("simY","ROOT::VecOps::RVec<float> y;for(auto& mom:simMom) y.push_back(mom.Rapidity()); return y;")
          .Filter("vtxChi2>0.0001"); // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("centrality|bNorm"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex|Y|Pid|IsProton|MotherId|Charge)"),
                                            std::regex("sim(Pt|Eta|Phi|Pdg|MotherId|Y)")
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"bNorm", 10, 0, 2} );

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

  std::vector<Qn::AxisD> proton_axes{
          { "trY", 12, -0.2, 1.0 },
          { "trPt", 12, 0.0, 1.2 },
  };

  VectorConfig proton( "proton", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "trPid", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  proton.AddCut( "trMotherId", [](double mother_id){
    auto int_mother_id = static_cast<int>(mother_id);
    return int_mother_id == -1;
    }, "cut on primary" );
  proton.AddHisto2D({{"trY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}}, "trIsProton");
  correction_task.AddVector(proton);

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


  correction_task.Run();
}