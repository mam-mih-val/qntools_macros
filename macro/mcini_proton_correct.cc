auto//
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
  auto psi_rp_function = [ &distribution, &engine ] ( double psi_rp ){
    return distribution(engine);
  };
  auto phi_function = []( float psi_rp, ROOT::VecOps::RVec<double> vec_px, ROOT::VecOps::RVec<double> vec_py ){
    ROOT::VecOps::RVec<float> vec_phi{};
    vec_phi.reserve( vec_px.size() );
    for( size_t i=0; i<vec_px.size(); ++i ){
      auto px = vec_px.at(i);
      auto py = vec_py.at(i);
      auto px_new = px * cos( psi_rp ) - py * sin( psi_rp );
      auto py_new = px * sin( psi_rp ) + py * cos( psi_rp );
      auto phi = atan2( py_new, px_new );
      vec_phi.push_back( phi );
    }
    return vec_phi;
  };
  auto pT_function = []( ROOT::VecOps::RVec<double> vec_px, ROOT::VecOps::RVec<double> vec_py ){
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
  auto ycm_function = []( ROOT::VecOps::RVec<double> vec_e, ROOT::VecOps::RVec<double> vec_pz ){
    ROOT::VecOps::RVec<float> vec_ycm{};
    vec_ycm.reserve( vec_e.size() );
    for( size_t i=0; i<vec_e.size(); ++i ){
      auto pz = vec_pz.at(i);
      auto E = vec_e.at(i);
      auto ycm = 0.5 * log( E + pz ) - 0.5 * log( E - pz );
      vec_ycm.push_back( ycm );
    }
    return vec_ycm;
  };
  auto eta_function = [Y_CM]( ROOT::VecOps::RVec<double> vec_px, ROOT::VecOps::RVec<double> vec_py, ROOT::VecOps::RVec<double> vec_pz ){
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
  
  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  
  auto dd=d
          .Filter( "b > 0." )
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("b|psi_rp"));
  correction_task.SetTrackVariables({
                                      std::regex("phi|pT|y|pdg|mate"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( { "b", 16, 0, 16 } );

  VectorConfig psi_rp( "psi_rp", "psi_rp", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( { 1, 2 } );
  psi_rp.SetCorrections( { CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  std::vector<Qn::AxisD> tru_proton_axes{
        { "y", 20, -0.5, 1.5 },
        { "pT", 10, 0.0, 2.0 },
  };

  VectorConfig tru_proton( "proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( { 1, 2 } );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( tru_proton_axes );
  tru_proton.AddCut( "pdg", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  correction_task.AddVector(tru_proton);

  tru_proton.AddCut( "mate", [](double mate){
    auto mate_code = static_cast<int>(mate);
    return mate_code != 1 && mate_code != -1;
    }, "mate cut" );
  correction_task.AddVector(tru_proton);

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}