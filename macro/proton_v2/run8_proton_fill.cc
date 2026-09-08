//
// Created by Misha on 3/7/2023.
//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <random>
#include <vector>
#include <functional>

#include <DataContainer.hpp>
#include <QnDataFrame.hpp>
#include "correlation_helper.h"
#include "bmn_env.h"
#include "vector_generators.h"

void run8_proton_fill( std::string file_list, 
                          std::string str_run_id_efficiency_file,
                          std::string str_effieciency_file,
                          std::string centrality_calib_file,
                          std::string str_pid_tof400_file,
                          std::string str_pid_tof700_file ){

  std::cout << "starting execution" << std::endl;

  auto proton_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "runId", 12, 7100, 8300 },
    Qn::AxisD{ "y", 12, 0.0, 1.2 },
    Qn::AxisD{ "pT", 10, 0.0, 2.0 },
  };

  auto qvector_axes = std::vector<Qn::AxisD>{
    Qn::AxisD{ "centrality", 6, 0, 60 },
    Qn::AxisD{ "runId", 12, 7100, 8300 },
  };

  auto calibration = DataCalibration{};

  calibration.selected_runs = std::vector<int>{ 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199 };
  calibration.rejected_runs = std::vector<int>{7313, 7415, 7417, 7435, 7469, 7517, 7519, 7520, 7537, 7575, 7604, 7630, 7657, 7659, 7679, 7681, 7705, 7735, 7843, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7907, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7954, 7955, 8031, 8032, 8033, 8115, 8121, 8167, 8201, 8204, 8205, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8247, 8265, 8266, 8267, 8281, 8289};

  auto file_pid400 = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( str_pid_tof400_file.c_str(), "READ" ), [](auto f){f->Close(); } };
  auto file_pid700 = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( str_pid_tof700_file.c_str(), "READ" ), [](auto f){f->Close(); } };

  assert(file_pid400);
  assert(file_pid700);

  file_pid400->GetObject( "fit_2212_x0", calibration.f1_2212_m_400 );
  file_pid700->GetObject( "fit_2212_x0", calibration.f1_2212_m_700 );
  file_pid400->GetObject( "fit_2212_sigma", calibration.f1_2212_s_400 );
  file_pid700->GetObject( "fit_2212_sigma", calibration.f1_2212_s_700 );

  auto file_fit = TFile::Open( centrality_calib_file.c_str(), "READ" );
	file_fit->cd();
	calibration.g1_FitVtxX = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_x");
	calibration.g1_FitVtxY = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_y");
	calibration.g1_FitVtxZ = file_fit->Get<TGraphErrors>("grNew_def_h2_RunId_vtx_z");

	calibration.g1_FitRunIdFactor = file_fit->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_nTracks_7400_7450");

  auto file_run_id_eff = std::unique_ptr< TFile, std::function< void(TFile*) > >{ TFile::Open( str_run_id_efficiency_file.c_str(), "READ" ), [](auto f){f->Close(); } };
  assert(file_run_id_eff);
	file_run_id_eff->GetObject("hn_efficiency", calibration.efficiency_eta_pT_phi_run_id);

  auto effieciency_file = std::unique_ptr<TFile>{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  effieciency_file->GetObject("h3_efficiency_2212_good", calibration.efficiency_histo);

  if( !calibration.efficiency_histo )
    std::cerr << "Warning: No efficiency for both tof was found in file " << str_effieciency_file << "\n";

  std::vector<int> f1_mod = {
    6,  7,  8,
    11, 12, 13,
    16,     17,
    20, 21, 22, 
    25, 26, 27
  };
  std::vector<int> f2_mod = {
    0,  1,  2,  3,  4,
    5,              9,
    10,             14,
    15,             18,
    19,             23,
    24,             28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_mod = {
    35,                 44,
    37,                 46, 
    39,                 48, 
    41,                 50,
    43,                 52
  };

  std::vector<int> f4_mod = {
    34,                     45,
    36,                     47, 
    38,                     49, 
    40,                     51,
    42,                     53
  };

  // std::for_each( f1_mod.begin(), f1_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f2_mod.begin(), f2_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f3_mod.begin(), f3_mod.end(), [](auto& m){ m += 1; } );
  // std::for_each( f4_mod.begin(), f4_mod.end(), [](auto& m){ m += 1; } );

  auto harmonics = std::vector<size_t>( 5 );
  std::iota( harmonics.begin(), harmonics.end(), 1 );

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", file_list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  
  auto dd = GenerateBmnExtendedTreeData(d, calibration);
  auto sampled_d = Qn::Correlation::Resample(dd, 100);

  sampled_d = sampled_d.Define( "F1w", fhcal_weight_generator(f1_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F2w", fhcal_weight_generator(f2_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F3w", fhcal_weight_generator(f3_mod), { "fhcalModId", "fhcalModE" } );
  sampled_d = sampled_d.Define( "F4w", fhcal_weight_generator(f4_mod), { "fhcalModId", "fhcalModE" } );

  auto p_components_names = AddUVectorComponents(sampled_d, "proton", harmonics, "trPhi" );
  auto p_cov_names = AddUVectorCovariance(sampled_d, "proton", harmonics, "trPhi" );

  auto f1_means_str = AddUVectorComponents( sampled_d, "F1", harmonics, "fhcalModPhi" );
  auto f2_means_str = AddUVectorComponents( sampled_d, "F2", harmonics, "fhcalModPhi" );
  auto f3_means_str = AddUVectorComponents( sampled_d, "F3", harmonics, "fhcalModPhi" );
  auto f4_means_str = AddUVectorComponents( sampled_d, "F4", harmonics, "fhcalModPhi" );
  auto tp_means_str = AddUVectorComponents( sampled_d, "Tpos", harmonics, "trPhi" );
  auto tn_means_str = AddUVectorComponents( sampled_d, "Tneg", harmonics, "trPhi" );

  auto f1_cov_str = AddUVectorCovariance( sampled_d, "F1", harmonics, "fhcalModPhi" );
  auto f2_cov_str = AddUVectorCovariance( sampled_d, "F2", harmonics, "fhcalModPhi" );
  auto f3_cov_str = AddUVectorCovariance( sampled_d, "F3", harmonics, "fhcalModPhi" );
  auto f4_cov_str = AddUVectorCovariance( sampled_d, "F4", harmonics, "fhcalModPhi" );
  auto tp_cov_str = AddUVectorCovariance( sampled_d, "Tpos", harmonics, "trPhi" );
  auto tn_cov_str = AddUVectorCovariance( sampled_d, "Tneg", harmonics, "trPhi" );

  auto p_components_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  p_components_ptr.reserve( p_components_names.size() );
  auto p_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  p_cov_ptr.reserve( p_cov_names.size() );

  auto f1_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f1_means_ptr.reserve( f1_means_str.size() );
  auto f1_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f1_cov_ptr.reserve( f1_cov_str.size() );

  auto f2_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f2_means_ptr.reserve( f2_means_str.size() );
  auto f2_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f2_cov_ptr.reserve( f2_cov_str.size() );

  auto f3_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f3_means_ptr.reserve( f3_means_str.size() );
  auto f3_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f3_cov_ptr.reserve( f3_cov_str.size() );

  auto f4_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f4_means_ptr.reserve( f4_means_str.size() );
  auto f4_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  f4_cov_ptr.reserve( f4_cov_str.size() );

  auto tp_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tp_means_ptr.reserve( tp_means_str.size() );
  auto tp_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tp_cov_ptr.reserve( tp_cov_str.size() );

  auto tn_means_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tn_means_ptr.reserve( tn_means_str.size() );
  auto tn_cov_ptr = std::vector< ROOT::RDF::RResultPtr< Qn::DataContainerStatCollect > >{};
  tn_cov_ptr.reserve( tn_cov_str.size() );

  for( const auto& name : p_components_names ){
    p_components_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), std::vector<std::string>{name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : p_cov_names ){
    p_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> >( CorrelationHelper(proton_axes), std::vector<std::string>{name, "trProtonWeight", "samples", "centrality", "trProtonY", "trPt" } )
    ); 
  }

  for( const auto& name : f1_means_str ){
    f1_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F1w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_means_str ){
    f2_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_means_str ){
    f3_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F3w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f4_means_str ){
    f4_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F3w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tp_means_str ){
    tp_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTposW", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tn_means_str ){
    tn_means_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTnegW", "samples", "centrality" } )
    ); 
  }


  for( const auto& name : f1_cov_str ){
    f1_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F1w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f2_cov_str ){
    f2_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f3_cov_str ){
    f3_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : f4_cov_str ){
    f4_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "F2w", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tp_cov_str ){
    tp_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTposW", "samples", "centrality" } )
    ); 
  }

  for( const auto& name : tn_cov_str ){
    tn_cov_ptr.emplace_back(
      sampled_d.Book< std::vector<double>, std::vector<double>,  ROOT::VecOps::RVec<ULong64_t>, float>( CorrelationHelper(qvector_axes), std::vector<std::string>{name, "trTnegW", "samples", "centrality" } )
    ); 
  }

  auto file_out = std::unique_ptr<TFile, std::function<void(TFile*)> >{ TFile::Open( "corr.root", "RECREATE"), [](auto f){ f->Close(); } };
  file_out->cd();
  std::for_each( p_components_ptr.begin(), p_components_ptr.end(), [i=0, &p_components_names]( auto& p ) mutable { p->Write( p_components_names.at(i).c_str() ); ++i; } );
  std::for_each( p_cov_ptr.begin(), p_cov_ptr.end(), [i=0, &p_cov_names]( auto& p ) mutable { p->Write( p_cov_names.at(i).c_str() ); ++i; } );
  
  std::for_each( f1_means_ptr.begin(), f1_means_ptr.end(), [i=0, &f1_means_str]( auto& p ) mutable { p->Write( f1_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f1_cov_ptr.begin(), f1_cov_ptr.end(), [i=0, &f1_cov_str]( auto& p ) mutable { p->Write( f1_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( f2_means_ptr.begin(), f2_means_ptr.end(), [i=0, &f2_means_str]( auto& p ) mutable { p->Write( f2_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f2_cov_ptr.begin(), f2_cov_ptr.end(), [i=0, &f2_cov_str]( auto& p ) mutable { p->Write( f2_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( f3_means_ptr.begin(), f3_means_ptr.end(), [i=0, &f3_means_str]( auto& p ) mutable { p->Write( f3_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f3_cov_ptr.begin(), f3_cov_ptr.end(), [i=0, &f3_cov_str]( auto& p ) mutable { p->Write( f3_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( f4_means_ptr.begin(), f4_means_ptr.end(), [i=0, &f4_means_str]( auto& p ) mutable { p->Write( f4_means_str.at(i).c_str() ); ++i; } );
  std::for_each( f4_cov_ptr.begin(), f4_cov_ptr.end(), [i=0, &f4_cov_str]( auto& p ) mutable { p->Write( f4_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( tp_means_ptr.begin(), tp_means_ptr.end(), [i=0, &tp_means_str]( auto& p ) mutable { p->Write( tp_means_str.at(i).c_str() ); ++i; } );
  std::for_each( tp_cov_ptr.begin(), tp_cov_ptr.end(), [i=0, &tp_cov_str]( auto& p ) mutable { p->Write( tp_cov_str.at(i).c_str() ); ++i; } );

  std::for_each( tn_means_ptr.begin(), tn_means_ptr.end(), [i=0, &tn_means_str]( auto& p ) mutable { p->Write( tn_means_str.at(i).c_str() ); ++i; } );
  std::for_each( tn_cov_ptr.begin(), tn_cov_ptr.end(), [i=0, &tn_cov_str]( auto& p ) mutable { p->Write( tn_cov_str.at(i).c_str() ); ++i; } );

  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}