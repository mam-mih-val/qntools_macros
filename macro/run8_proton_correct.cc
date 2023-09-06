//
// Created by Misha on 3/7/2023.
//

void run8_proton_correct(std::string list){
  
  const float PROTON_M = 0.938;
  const float Y_CM = 1.15141;

  auto f1_m2_400 = new TF1( "m2_p_400", "pol2", 0, 10 );
  f1_m2_400->SetParameter( 0, 0.965499 );
  f1_m2_400->SetParameter( 1, -0.0625193 );
  f1_m2_400->SetParameter( 2, -0.0217673 );

	auto f1_s_400 = new TF1( "s_p_400", "pol3", 0, 10 );
  f1_s_400->SetParameter( 0, 0.220837 );
  f1_s_400->SetParameter( 1, -0.214113 );
  f1_s_400->SetParameter( 2, 0.161722 );
  f1_s_400->SetParameter( 3, -0.0251886 );
  
  auto f1_m2_700 = new TF1( "m2_p_700", "pol4", 0, 10 );
  f1_m2_700->SetParameter( 0, 1.0847 );
  f1_m2_700->SetParameter( 1, -0.330513 );
  f1_m2_700->SetParameter( 2, 0.220286 );
  f1_m2_700->SetParameter( 3, -0.064973 );
  f1_m2_700->SetParameter( 4, 0.00705849 );

  auto f1_s_700 = new TF1( "s_p_700", "pol3", 0, 10 );
  f1_s_700->SetParameter( 0,  0.102933 );
  f1_s_700->SetParameter( 1,  -0.115384 );
  f1_s_700->SetParameter( 2,  0.088186 );
  f1_s_700->SetParameter( 3,  -0.0115386 );
    
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
          .Define( "trDcaX", " std::vector<float> vec_par; for( auto par : trParamFirst ){ vec_par.push_back( par.at(0) - vtxX ); } return vec_par; " )
		      .Define( "trDcaY", " std::vector<float> vec_par; for( auto par : trParamFirst ){ vec_par.push_back( par.at(1) - vtxY ); } return vec_par; " )
          .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
          .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P()/trCharge.at(i) ); } return pq;" )
          .Define( "trM2Tof400","std::vector<float> vec_m2;\n"
                    "for( int i=0; i<trMom.size(); i++ ){\n"
                    " auto p = trMom.at(i).P();\n"
                    " auto p2 = p*p;\n"
                    " auto beta = trBetaTof400.at(i);\n"
                    " auto beta2 = beta*beta;\n"
                    " auto gamma2 = 1 - beta2;\n"
                    " auto m2 = beta > -990. ? p2 / beta2 * gamma2 : -999.0;\n"
                    " vec_m2.push_back( m2 );\n"
                    "}\n"
                    "return vec_m2;" )
          .Define( "trM2Tof700","std::vector<float> vec_m2;\n"
                    "for( int i=0; i<trMom.size(); i++ ){\n"
                    " auto p = trMom.at(i).P();\n"
                    " auto p2 = p*p;\n"
                    " auto beta = trBetaTof700.at(i);\n"
                    " auto beta2 = beta*beta;\n"
                    " auto gamma2 = 1 - beta2;\n"
                    " auto m2 = beta > -990. ? p2 / beta2 * gamma2 : -999.0;\n"
                    " vec_m2.push_back( m2 );\n"
                    "}\n"
                    "return vec_m2;" )
          .Define( "trIsProton400", [ f1_m2_400, f1_s_400 ]( 
                std::vector<float> vec_m2, 
                std::vector<float> vec_pq 
                ){
                  std::vector<int> vec_is{};
                  vec_is.reserve( vec_pq.size() );
                  for( int i=0; i < vec_pq.size(); ++i ){
                    auto pq = vec_pq.at(i);
                    auto m2 = vec_m2.at(i);
                    if( pq < 0 ){ vec_is.push_back(0); continue; }
                    auto mean = f1_m2_400->Eval(pq);
                    auto sigma = f1_s_400->Eval(pq);
                    auto lo = mean - 2*sigma;
                    auto hi = mean + 2*sigma;
                    vec_is.push_back( lo < m2 && m2 < hi ? 1 : 0 );
                  }
                  return vec_is;
                }
                ,{ "trM2Tof400", "pq" } )
          .Define( "trIsProton700", [ f1_m2_700, f1_s_700 ]( 
                std::vector<float> vec_m2, 
                std::vector<float> vec_pq 
                ){
                  std::vector<int> vec_is{};
                  vec_is.reserve( vec_pq.size() );
                  for( int i=0; i < vec_pq.size(); ++i ){
                    auto pq = vec_pq.at(i);
                    auto m2 = vec_m2.at(i);
                    if( pq < 0 ){ vec_is.push_back(0); continue; }
                    auto mean = f1_m2_700->Eval(pq);
                    auto sigma = f1_s_700->Eval(pq);
                    auto lo = mean - 2*sigma;
                    auto hi = mean + 2*sigma;
                    vec_is.push_back( lo < m2 && m2 < hi ? 1 : 0 );
                  }
                  return vec_is;
                }
                ,{ "trM2Tof700", "pq" } )
          .Define( "trIsProton", []( std::vector<int> is_400, std::vector<int> is_700 ){
                  std::vector<int> vec_is{};
                  vec_is.reserve( is_400.size() );
                  for( int i=0; i<is_400.size(); ++i ){ vec_is.push_back( is_400.at(i) == 1 || is_700.at(i) == 1 ? 1 : 0 ); }
                  return vec_is;
                }, {"trIsProton400", "trIsProton700"} )
          .Define( "trProtonY", [PROTON_M, Y_CM]( std::vector<float> vec_pz, std::vector<float> vec_pq ){
                  std::vector<float> vec_y{};
                  vec_y.reserve( vec_pz.size() );
                  for( int i=0; i<vec_pz.size(); ++i ){
                    auto pz = vec_pz.at(i);
                    auto p = vec_pq.at(i);
                    auto E = sqrt( p*p + PROTON_M*PROTON_M );
                    auto y = 0.5 * log( ( E + pz )/( E - pz ) ) - Y_CM;
                    vec_y.push_back( y );
                  }
                  return vec_y;
                }, {"pz", "pq"} )
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
          .Filter("1e4 < bc1Integral.at(0) && bc1Integral.at(0) < 4e4" )
          .Filter("vtxChi2/vtxNdf > 0.1")
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", "qa.root" );
  correction_task.SetEventVariables(std::regex("centrality"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Eta|Phi|IsProton|Charge|ProtonY|IsProton|DcaX|DcaY)"),
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
  Tpos.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  Tpos.AddCut( "trEta", [](double eta){
    return 2.0 < eta && eta < 3.0;
  }, "eta cut" );
  Tpos.AddCut( "trPt", [](double pT){
    return pT > 0.2;
  }, "pT cut" );
  correction_task.AddVector(Tpos);

  std::vector<Qn::AxisD> proton_axes{
        { "trProtonY", 12, -0.2, 1.0 },
        { "trPt", 15, 0.0, 1.5 },
  };

  VectorConfig proton( "proton", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton.SetHarmonicArray( {1, 2} );
  proton.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::RESCALING } );
  proton.SetCorrectionAxes( proton_axes );
  proton.AddCut( "trIsProton", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 1;
    }, "proton cut" );
  proton.AddCut( "trDcaX", [](double dca){
    return fabs(dca) < 3.0;
    }, "cut on dca x" );
  proton.AddCut( "trDcaY", [](double dca){
    return fabs(dca) < 5.0;
    }, "cut on dca y" );
  proton.AddHisto2D({{"trProtonY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}}, "trIsProton");
  correction_task.AddVector(proton);

  correction_task.Run();
}