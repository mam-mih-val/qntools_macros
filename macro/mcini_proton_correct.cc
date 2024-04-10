//
// Created by Misha on 3/7/2023.
//

#include <cmath>
#include <math.h>
#include <vector>

struct PidFunctions{
	TF1* mean_400;
	TF1* sigma_400;
	TF1* mean_700;
	TF1* sigma_700;
};


void mcini_proton_correct(  std::string list, 
                          std::string calib_in_file="qa.root"){
  
  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

  std::random_device dev{}; 
  std::mt19937 engine(dev()); 
  std::uniform_real_distribution<float> distribution{M_PI*-1, M_PI};
  auto psi_rp_function = [ &distribution, &engine ] ( float psi_rp ){
    return distribution(engine);
  };
  auto phi_function = []( float psi_rp, ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py ){
    ROOT::VecOps::RVec<float> vec_phi{};
    vec_phi.reserve( vec_px.size() );
    for( size_t i=0; i<vec_px.size(); ++i ){
      auto px = vec_px.at(i);
      auto py = vec_py.at(i);
      auto phi = atan2( py, px );
      vec_phi.push_back( phi+psi_rp );
    }
    return vec_phi;
  };
  auto pT_function = []( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py ){
    ROOT::VecOps::RVec<float> vec_pT{};
    vec_pT.reserve( vec_py.size() );
    for( size_t i=0; i<vec_py.size(); ++i ){
      auto px = vec_px.at(i);
      auto py = vec_py.at(i);
      auto pT = sqrt(px*px + py*py);
      vec_pT.push_back( pT );
    }
    return vec_pT;
  };
  auto ycm_function = [Y_CM]( ROOT::VecOps::RVec<float> vec_e, ROOT::VecOps::RVec<float> vec_pz ){
    ROOT::VecOps::RVec<float> vec_ycm{};
    vec_ycm.reserve( vec_e.size() );
    for( size_t i=0; i<vec_e.size(); ++i ){
      auto pz = vec_pz.at(i);
      auto E = vec_e.at(i);
      auto ycm = log( E + pz ) - log( E - pz ) - Y_CM;
      vec_ycm.push_back( ycm );
    }
    return vec_ycm;
  };
  auto eta_function = [Y_CM]( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py, ROOT::VecOps::RVec<float> vec_pz ){
    ROOT::VecOps::RVec<float> vec_eta{};
    vec_eta.reserve( vec_px.size() );
    for( size_t i=0; i<vec_px.size(); ++i ){
      auto px = vec_px.at(i);
      auto py = vec_py.at(i);
      auto pT = sqrt(px*px + py*py);
      auto pz = vec_pz.at(i);
      auto theta = atan2(pT, pz);
      auto eta = - log( tan( theta /2 ) );
      vec_eta.push_back( eta );
    }
    return vec_eta;
  };
  auto ekin_function = [Y_CM]( ROOT::VecOps::RVec<float> vec_px, ROOT::VecOps::RVec<float> vec_py, ROOT::VecOps::RVec<float> vec_pz, ROOT::VecOps::RVec<float> vec_e ){
    ROOT::VecOps::RVec<float> vec_ekin{};
    vec_ekin.reserve( vec_px.size() );
    for( size_t i=0; i<vec_px.size(); ++i ){
      auto E = vec_e.at(i);
      auto px = vec_px.at(i);
      auto py = vec_py.at(i);
      auto pz = vec_pz.at(i);
      auto M = sqrt( E*E - px*px - py*py - pz*pz );
      auto Ekin = E - M;
      vec_ekin.push_back( Ekin );
    }
    return vec_ekin;
  };
  
  TStopwatch timer;
  timer.Start();
  std::string treename = "events";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  
  auto dd=d
          .Alias( "b", "fB" )
          .Define( "psi_rp", psi_rp_function, {"fPhi"} )
          .Define( "phi", phi_function, { "psi_rp", "event.fParticles.fPx", "event.fParticles.fPy" } )
          .Define( "pT", pT_function, { "event.fParticles.fPx", "event.fParticles.fPy" } )
          .Define( "y", ycm_function, { "event.fParticles.fE", "event.fParticles.fPz" } )
          .Define( "eta", eta_function, { "event.fParticles.fPx", "event.fParticles.fPy", "event.fParticles.fPz" } )
          .Define( "Ekin", ekin_function, {  "event.fParticles.fPx", "event.fParticles.fPy", "event.fParticles.fPz", "event.fParticles.fE" } )
          // .Alias( "pdg", "event.fParticles.fPdg"  )
          .Filter( "b > 0." )
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("fB|psi_rp"));
  correction_task.SetTrackVariables({
                                      std::regex("phi|pT|y|Ekin|pdg"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"b", 7, 0, 14} );

  VectorConfig psi_rp( "psi_rp", "psi_rp", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( {1, 2, 3} );
  psi_rp.SetCorrections( {CORRECTION::PLAIN } );
  // psi_rp.AddHisto1D({"psiRp", 100, -3.5, 3.5}, "psiRP");
  correction_task.AddVector(psi_rp);

  std::vector<Qn::AxisD> tru_proton_axes{
        { "y", 10, -0.6, 1.4 },
        { "pT", 10, 0.0, 2.0 },
  };

  VectorConfig tru_proton( "tru_proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( {1, 2, 3} );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( tru_proton_axes );
  tru_proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  tru_proton.AddCut( "eta", [](double eta){
    return eta < 2.7;
    }, "rapidity cut" );
  tru_proton.AddHisto2D({{"simProtonY", 100, -0.5, 1.5}, {"simPt", 100, 0.0, 2.0}}, "simIsProton");
  correction_task.AddVector(tru_proton);

  VectorConfig S1( "S1", "phi", "Ekin", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  S1.SetHarmonicArray( {1, 2, 3} );
  S1.SetCorrections( {CORRECTION::PLAIN } );
  S1.AddCut( "eta", [](double eta){
    return 3.8 < eta && eta < 5.4;
    }, "rapidity cut" );
  S1.AddHisto2D({{"simProtonY", 100, -0.5, 1.5}, {"simPt", 100, 0.0, 2.0}}, "simIsProton");
  correction_task.AddVector(S1);

  VectorConfig S2( "S2", "phi", "Ekin", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  S2.SetHarmonicArray( {1, 2, 3} );
  S2.SetCorrections( {CORRECTION::PLAIN } );
  S2.AddCut( "eta", [](double eta){
    return 3.3 < eta && eta < 3.8;
    }, "rapidity cut" );
  S2.AddHisto2D({{"simProtonY", 100, -0.5, 1.5}, {"simPt", 100, 0.0, 2.0}}, "simIsProton");
  correction_task.AddVector(S2);

  VectorConfig S3( "S3", "phi", "Ekin", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  S3.SetHarmonicArray( {1, 2, 3} );
  S3.SetCorrections( {CORRECTION::PLAIN } );
  S3.AddCut( "eta", [](double eta){
    return 2.7 < eta && eta < 3.3;
    }, "rapidity cut" );
  S3.AddHisto2D({{"simProtonY", 100, -0.5, 1.5}, {"simPt", 100, 0.0, 2.0}}, "simIsProton");
  correction_task.AddVector(S3);

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}