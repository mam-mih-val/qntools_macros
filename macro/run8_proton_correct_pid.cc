//
// Created by Misha on 3/7/2023.
//

#include <cmath>
#include <vector>

void run8_proton_correct_pid( std::string list, 
                          std::string str_effieciency_file,
                          std::string calib_in_file="qa.root" ){

  std::cout << "starting execution" << std::endl;
  
  const float PROTON_M = 0.938; // GeV/c2
  const float PI_POS_M = 0.134;
  const float DEUTERON_M = 1.875;  
  const float Y_CM = 1.15141;
  const float FHCAL_Z = 980; // cm

	auto f1_2212_m_400 = new TF1("2212_mean_400", "pol1");
  f1_2212_m_400->SetParameter(0, 0.94612766);
  f1_2212_m_400->SetParameter(1, -0.026582296);

  auto f1_2212_s_400 = new TF1("2212_mean_400", "pol2");
  f1_2212_s_400->SetParameter(0, 0.087635220);
  f1_2212_s_400->SetParameter(1, -0.028676680);
  f1_2212_s_400->SetParameter(2, 0.013336097);

  auto f1_2212_m_700 = new TF1("2212_mean_700", "pol1");
  f1_2212_m_700->SetParameter(0, 0.91900468);
  f1_2212_m_700->SetParameter(1, 0.00070073927);

  auto f1_2212_s_700 = new TF1("2212_mean_700", "pol2");
  f1_2212_s_700->SetParameter(0, 0.017564553);
  f1_2212_s_700->SetParameter(1, 0.022982477);
  f1_2212_s_700->SetParameter(2, 0.011826542);

  auto m2_function = 
  []
  ( ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >> vec_mom, 
    ROOT::VecOps::RVec<double> vec_beta){
      std::vector<float> vec_m2;
      vec_m2.reserve( vec_beta.size() );
      for( size_t i=0; i<vec_mom.size(); i++ ){
        auto p = vec_mom.at(i).P();
        auto p2 = p*p;
        auto beta = vec_beta.at(i);
        auto beta2 = beta*beta;
        auto gamma2 = 1 - beta2;
        auto m2 = beta > -990. ? p2 / beta2 * gamma2 : -999.0;
        vec_m2.push_back( m2 );
      }
      return vec_m2;
    };
  const auto n_sigma_generator = []( auto f1_mean, auto f1_sigma ){
    return 
    [ f1_mean, f1_sigma ]
    ( std::vector<float> vec_pq, std::vector<float> vec_m2 ){
        auto vec_n_sigma = std::vector<float>{};
        vec_n_sigma.reserve( vec_pq.size() );
        for( size_t i=0; i < vec_pq.size(); ++i ){
          auto m2 = vec_m2.at(i);
          auto pq = vec_pq.at(i);
          auto mean = f1_mean->Eval(pq);
          auto sigma = f1_sigma->Eval(pq);
          auto n_sigma = fabs( m2 - mean ) / sigma;
          vec_n_sigma.push_back( pq > 0 ? n_sigma : 999. );
        }
      return vec_n_sigma;
    };
  };
   
	const auto n_sigma_particle_function = 
  []
  ( std::vector<float> n_sigma_400, 
    std::vector<float> n_sigma_700 ){
      std::vector<int> vec_n_simga{};
      vec_n_simga.reserve( n_sigma_400.size() );
      for( int i=0; i<n_sigma_400.size(); ++i ){ 
        vec_n_simga.push_back( std::min( n_sigma_400.at(i), n_sigma_700.at(i) ) ); }
      return vec_n_simga;
  };
  
	const auto rapidity_generator = []( auto particle_m, auto y_cm ){
    return 
    [particle_m, y_cm]( std::vector<float> vec_pz, std::vector<float> vec_pq ){
      std::vector<float> vec_y{};
      vec_y.reserve( vec_pz.size() );
      for( int i=0; i<vec_pz.size(); ++i ){
        auto pz = vec_pz.at(i);
        auto p = vec_pq.at(i);
        auto E = sqrt( p*p + particle_m*particle_m );
        auto y = 0.5 * log( ( E + pz ) / ( E - pz ) ) - y_cm;
        vec_y.push_back( y );
      }
      return vec_y;
    };
  };
  const auto function_fhcal_x = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<std::vector<float>> vec_param ){
      std::vector<float> vec_x{};
      vec_x.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto x = par.at(0);
        auto z = par.at(2);
        auto tx = par.at(3);
        auto dz = FHCAL_Z - z;
        auto dx = tx * dz;
        vec_x.push_back( x+dx );
      }
      return vec_x;
    };
  const auto function_fhcal_y = 
  [FHCAL_Z]
  ( ROOT::VecOps::RVec<vector<float>> vec_param ){
      std::vector<float> vec_y{};
      vec_y.reserve( vec_param.size() );
      for( auto par : vec_param ){
        auto y = par.at(1);
        auto z = par.at(2);
        auto ty = par.at(4);
        auto dz = FHCAL_Z - z;
        auto dy = ty * dz;
        vec_y.push_back( y+dy );
      }
      return vec_y;
    };

  const auto centrality_function = 
  []
  (ROOT::VecOps::RVec<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >> vec_mom){
      float centrality;
      std::vector<float> centrality_percentage{ 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100 };
      std::vector<int> multiplicity_edges{ 248, 150, 122, 100, 82, 67, 54, 43, 34, 21, 13, 8, 4, 3, 1 };
      auto multiplicity = vec_mom.size();
      if( multiplicity > multiplicity_edges[0] )
        return -1.0f;
      int idx = 0;
      float bin_edge = multiplicity_edges[idx];
      while( multiplicity < bin_edge &&
        idx < multiplicity_edges.size()-1 ){
        idx++;
        bin_edge = multiplicity_edges[idx];
      }
      centrality = (centrality_percentage[idx-1] + centrality_percentage[idx])/2.0f;
      return centrality;
    };
  
  const auto dca_function = [](std::vector<float> vec_x, std::vector<float> vec_y){
    std::vector<float> vec_r{};
    vec_r.reserve(vec_x.size());
    for (int i=0; i<vec_x.size(); ++i) {
      auto x = vec_x.at(i);
      auto y = vec_y.at(i);
      auto r = std::sqrt( x*x + y*y );
      vec_r.push_back(r);
    }
    return vec_r;
  };

  const auto weight_generator = []( auto efficiency_map ){
    return [efficiency_map](std::vector<float> vec_y, ROOT::VecOps::RVec<float> vec_pT){
      if( !efficiency_map ){
          return std::vector<float>(vec_y.size(), 1);
        }
      std::vector<float> vec_weight{};
      vec_weight.reserve(vec_y.size());
      for( int i=0; i<vec_y.size(); ++i ){
        auto pT = vec_pT.at(i);
        auto y = vec_y.at(i);
        auto y_bin = efficiency_map->GetXaxis()->FindBin( y );
        auto pT_bin = efficiency_map->GetYaxis()->FindBin( pT );
        auto efficiency = efficiency_map->GetBinContent( y_bin, pT_bin );
        auto weight = efficiency > 1e-2 ? 1.0 / efficiency : 0.0;
        vec_weight.push_back( weight );
      }
      return vec_weight;
    };
  };

  std::unique_ptr<TFile> effieciency_file{TFile::Open( str_effieciency_file.c_str(), "READ" )};
  TH2D* efficiency_histo{nullptr};
  TH2D* efficiency_tof400{nullptr};
  TH2D* efficiency_tof700{nullptr};

  effieciency_file->GetObject("efficiency_2212_tof", efficiency_histo);
  if( !efficiency_histo )
    std::cerr << "Warning: No efficiency for both tof was found in file " << str_effieciency_file << "\n";
  effieciency_file->GetObject("efficiency_2212_tof400", efficiency_tof400);
  if( !efficiency_tof400 )
    std::cerr << "Warning: No efficiency for tof-400 was found in file " << str_effieciency_file << "\n";
  effieciency_file->GetObject("efficiency_2212_tof700", efficiency_tof700);
  if( !efficiency_tof700 )
    std::cerr << "Warning: No efficiency for tof-700 was found in file " << str_effieciency_file << "\n";
  
  std::vector<int> physical_runs{6667, 6668, 6669, 6670, 6671, 6672, 6673, 6674, 6675, 6676, 6677, 6678, 6679, 6680, 6681, 6683, 6684, 6685, 666, 6687, 6689, 6690, 6691, 6692, 6694, 6695, 6696, 6698, 6699, 6732, 6733, 6734, 6737, 6738, 6739, 6740, 6745, 6752, 6753, 6760, 6761, 6765, 6766, 6767, 6768, 6769, 6771, 6772, 6773, 6774, 6779, 6780, 6782, 6783, 6785, 6786, 6788, 6794, 6795, 6797, 6799, 6800, 6803, 6815, 6816, 6817, 6818, 6819, 6820, 6821, 6822, 6879, 6882, 6883, 6884, 6886, 6887, 6889, 6891, 6900, 6901, 6902, 6903, 6904, 6905, 6906, 6907, 6908, 6909, 6910, 6911, 6915, 6916, 6918, 6919, 6920, 6921, 6923, 6924, 6926, 6927, 6928, 6929, 6930, 6931, 6932, 6933, 6934, 6935, 6936, 6937, 6939, 6940, 6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6990, 6991, 6992, 6993, 6994, 6995, 6997, 6998, 6999, 7000, 7002, 7003, 7004, 7005, 7006, 7008, 7009, 7010, 7011, 7012, 7030, 7031, 7032, 7033, 7034, 7035, 7037, 7038, 7040, 7041, 7042, 7043, 7044, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7075, 7076, 7077, 7078, 7081, 7082, 7083, 7084, 7086, 7087, 7091, 7092, 7093, 7094, 7096, 7097, 7098, 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7842, 7843, 7845, 7846, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7941, 7942, 7944, 7948, 7949, 7950, 7952, 7954, 7955, 7957, 7958, 7960, 7961, 7962, 7963, 7965, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8026, 8027, 8028, 8029, 8030, 8031, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8217, 8219, 8220, 8221, 8228, 8229, 8230, 8231, 8235, 8236, 8238, 8239, 8240, 8242, 8244, 8245, 8246, 8247, 8248, 8250, 8251, 8253, 8254, 8255, 8256, 8257, 8258, 8265, 8266, 8267, 8268, 8270, 8271, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8281, 8284, 8286, 8287, 8288, 8289, 8290, 8292, 8293, 8294, 8295, 8297, 8298, 8299, 8300, 8305, 8306};

  std::vector<int> f1_modules = {
    34, 35, 
    36, 37, 
    38, 39, 
    40, 41, 
    42, 43, 
  };
  std::vector<int> f2_modules = {
     0,  1,  2,  3,  4,
     5,  6,  7,  8,  9,
    10, 11, 12, 13, 14,
    15, 16,     17, 18, 
    19, 20, 21, 22, 23,
    24, 25, 26, 27, 28,
    29, 30, 31, 32, 33,
  };
  std::vector<int> f3_modules = {
    44, 45,
    46, 47,
    48, 49,
    50, 51, 
    52, 53,
  };

  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  TFileCollection collection( "collection", "", list.c_str() );
  auto* chain = new TChain( treename.c_str() );
  chain->AddFileInfoList( collection.GetList() );
  ROOT::RDataFrame d( *chain );
  std::cout << "Preparing the RDF" << std::endl;
  auto dd=d
          .Define("track_multiplicity", "return trMom.size();")
          .Define( "stsNdigits","return stsDigits.size()" )
          .Define("centrality", centrality_function, {"trMom"} )
          .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
          .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
          .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
          .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
          .Define( "trDcaX", " std::vector<float> vec_par; for( auto par : trParamFirst ){ vec_par.push_back( par.at(0) - vtxX ); } return vec_par; " )
		      .Define( "trDcaY", " std::vector<float> vec_par; for( auto par : trParamFirst ){ vec_par.push_back( par.at(1) - vtxY ); } return vec_par; " )
          .Define( "trDcaR", dca_function, {"trDcaX", "trDcaY"} )
          .Define( "trFhcalX", function_fhcal_x, {"trParamLast"} )
          .Define( "trFhcalY", function_fhcal_y, {"trParamLast"} )
          .Define( "trChi2Ndf", " std::vector<float> vec_par; for( int i=0; i<trChi2.size(); ++i ){ vec_par.push_back( trChi2.at(i)/trNdf.at(i) ); } return vec_par; " )
          .Define( "trPx", " std::vector<float> px; for( auto mom : trMom ){ px.push_back( mom.Px() ); } return px; " )
          .Define( "trPy", " std::vector<float> py; for( auto mom : trMom ){ py.push_back( mom.Py() ); } return py; " )
          .Define( "pz", " std::vector<float> pz; for( auto mom : trMom ){ pz.push_back( mom.Pz() ); } return pz; " )
          .Define( "pq", " std::vector<float> pq; for( int i=0; i<trMom.size(); i++ ){ pq.push_back( trMom.at(i).P() / trCharge.at(i) ); } return pq;" )
          .Define( "trM2Tof700", m2_function, { "trMom", "trBetaTof700" } )
          .Define( "trM2Tof400", m2_function, { "trMom", "trBetaTof400" } )
          .Define( "trNsigmaProton400", n_sigma_generator(f1_2212_m_400, f1_2212_s_400), { "pq", "trM2Tof400" } )
          .Define( "trNsigmaProton700", n_sigma_generator(f1_2212_m_700, f1_2212_s_700), { "pq", "trM2Tof700"  } )
          .Define( "trNsigmaProton", n_sigma_particle_function, {"trNsigmaProton400", "trNsigmaProton700"} )
          .Define( "trProtonY", rapidity_generator(PROTON_M, Y_CM), {"pz", "pq"} )
          .Define( "trWeight", weight_generator(efficiency_histo), {"trProtonY", "trPt"} )
          .Define( "trWeightTof400", weight_generator(efficiency_tof400), {"trProtonY", "trPt"} )
          .Define( "trWeightTof700", weight_generator(efficiency_tof700), {"trProtonY", "trPt"} )
          .Alias("trStsNhits", "stsTrackNhits")
          .Alias("trStsChi2", "stsTrackChi2Ndf")
          .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom : trMom) eta.push_back(mom.eta()); return eta;")
          .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom : trMom) phi.push_back(mom.phi()); return phi;")
          // .Filter([]( ROOT::VecOps::RVec<float> bc1_int, ROOT::VecOps::RVec<float> fd_int ){
          //   float trigger = fd_int[0]-bc1_int[0]*0.602;
          //   return -25900 < trigger && trigger < -6030;
          // }, {"bc1sIntegral", "fdIntegral"} )
          .Filter([&physical_runs]( UInt_t run_id ){ 
            return std::find( physical_runs.begin(), physical_runs.end(), run_id) != physical_runs.end();
          }, {"runId"} )
          .Filter("runId < 8312")
          .Filter( []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }, {"triggerMapAR"} )
          .Filter([]( unsigned long sts_digits, unsigned long n_tracks ){ 
            double sts_min = sts_digits-n_tracks*(4.81632+0.0332792*n_tracks-9.62078e-05*n_tracks*n_tracks);
            double sts_max = sts_digits-n_tracks*(19.4203-0.0518774*n_tracks+4.56033e-05*n_tracks*n_tracks);
            return -74.0087 < sts_min && sts_max < 188.248; 
          }, {"stsNdigits", "track_multiplicity"} )
          .Filter("vtxChi2/vtxNdf > 0.1")
          .Filter("fabs(vtxX)<1.0")
          .Filter("fabs(vtxY)<1.0")
          .Filter("fabs(vtxZ)<0.5")
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("centrality"));
  correction_task.SetChannelVariables({std::regex("fhcalMod(X|Y|Phi|E|Id)")});
  correction_task.SetTrackVariables({
                                            std::regex("tr(Pt|Px|Py|Eta|Phi|NsigmaProton|NsigmaProton400|NsigmaProton700|Charge|ProtonY|DcaR|Chi2Ndf|Nhits|Weight|WeightTof400|WeightTof700|FhcalX|FhcalY|StsNhits|StsChi2)"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( {"centrality", 4, 0, 40} );

  VectorConfig f1( "F1", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f1.SetHarmonicArray( {1, 2} );
  f1.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f1.AddCut( "fhcalModId", [&f1_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f1_modules.begin(), f1_modules.end(), id) != f1_modules.end();
    }, "F1 Cut" );
  f1.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f1);

  VectorConfig f2( "F2", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f2.SetHarmonicArray( {1, 2} );
  f2.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f2.AddCut( "fhcalModId", [&f2_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f2_modules.begin(), f2_modules.end(), id) != f2_modules.end();
    }, "F2 Cut" );
  f2.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f2);

  VectorConfig f3( "F3", "fhcalModPhi", "fhcalModE", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  f3.SetHarmonicArray( {1, 2} );
  f3.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  f3.AddCut( "fhcalModId", [&f3_modules](double mod_id){
    auto id = static_cast<int>(mod_id);
    return std::find( f3_modules.begin(), f3_modules.end(), id) != f3_modules.end();
    }, "F3 Cut" );
  f3.AddHisto2D({{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}});
  correction_task.AddVector(f3);

  VectorConfig Tneg( "Tneg", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tneg.SetHarmonicArray( {1, 2} );
  Tneg.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  Tneg.AddCut( "trCharge", [](double charge){
    return charge < 0.0;
    }, "charge" );
  Tneg.AddCut( "trEta", [](double eta){
    return 1.5 < eta && eta < 4.0;
    }, "eta cut" );
  Tneg.AddCut( "trPt", [](double pT){
    return pT > 0.2;
    }, "pT cut" );
  Tneg.AddCut( "trFhcalX", [](double pos){
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tneg.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tneg);

  VectorConfig Tpos( "Tpos", "trPhi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  Tpos.SetHarmonicArray( {1, 2, 3} );
  Tpos.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  Tpos.AddCut( "trCharge", [](double charge){
    return charge >= 0.0;
    }, "charge" );
  Tpos.AddCut( "trEta", [](double eta){
    return 2.0 < eta && eta < 3.0;
  }, "eta cut" );
  Tpos.AddCut( "trPt", [](double pT){
    return pT > 0.2;
  }, "pT cut" );
  Tpos.AddCut( "trFhcalX", [](double pos){
    return pos < -40.0 || pos > 170;
    }, "cut on x-pos in fhcal plane" );
  Tpos.AddCut( "trFhcalY", [](double pos){
    return pos < -100.0 || pos > 100;
    }, "cut on y-pos in fhcal plane" );
  correction_task.AddVector(Tpos);

  std::vector<Qn::AxisD> proton_axes{
        { "trProtonY", 12, -0.2, 1.0 },
        { "trPt", 10, 0.0, 2.0 },
  };
  
  VectorConfig proton_1sigma( "proton_1sigma", "trPhi", "trWeight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton_1sigma.SetHarmonicArray( {1, 2, 3} );
  proton_1sigma.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  proton_1sigma.SetCorrectionAxes( proton_axes );
  proton_1sigma.AddCut( "trNsigmaProton", [](double n_sigma){
    return n_sigma < 1;
  }, "proton_1sigma cut" );
  proton_1sigma.AddCut( "trFhcalX", [](double pos){
    return pos < -30.0 || pos > 160;
  }, "cut on x-pos in fhcal plane" );
  proton_1sigma.AddCut( "trFhcalY", [](double pos){
    return pos < -60.0 || pos > 60;
  }, "cut on y-pos in fhcal plane" );
  proton_1sigma.AddCut( "trStsNhits", [](double nhits){
    return nhits > 5.5;
  }, "cut on fake tracks" );
  proton_1sigma.AddCut( "trDcaR", [](double dca){
    return dca < 5.0;
  }, "DCA cut" );
  proton_1sigma.AddCut( "trStsChi2", [](double chi2){
    return chi2 < 5.0;
  }, "cut on chi2 in sts" );
  proton_1sigma.AddHisto2D({{"trProtonY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}});
  correction_task.AddVector(proton_1sigma);

  VectorConfig proton_2sigma( "proton_2sigma", "trPhi", "trWeight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton_2sigma.SetHarmonicArray( {1, 2, 3} );
  proton_2sigma.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  proton_2sigma.SetCorrectionAxes( proton_axes );
  proton_2sigma.AddCut( "trNsigmaProton", [](double n_sigma){
    return n_sigma < 2;
  }, "proton_2sigma cut" );
  proton_2sigma.AddCut( "trFhcalX", [](double pos){
    return pos < -30.0 || pos > 160;
  }, "cut on x-pos in fhcal plane" );
  proton_2sigma.AddCut( "trFhcalY", [](double pos){
    return pos < -60.0 || pos > 60;
  }, "cut on y-pos in fhcal plane" );
  proton_2sigma.AddCut( "trStsNhits", [](double nhits){
    return nhits > 5.5;
  }, "cut on fake tracks" );
  proton_2sigma.AddCut( "trDcaR", [](double dca){
    return dca < 5.0;
  }, "DCA cut" );
  proton_2sigma.AddCut( "trStsChi2", [](double chi2){
    return chi2 < 5.0;
  }, "cut on chi2 in sts" );
  proton_2sigma.AddHisto2D({{"trProtonY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}});
  correction_task.AddVector(proton_2sigma);

  VectorConfig proton_3sigma( "proton_3sigma", "trPhi", "trWeight", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  proton_3sigma.SetHarmonicArray( {1, 2, 3} );
  proton_3sigma.SetCorrections( {CORRECTION::PLAIN, CORRECTION::RECENTERING, CORRECTION::TWIST_RESCALING } );
  proton_3sigma.SetCorrectionAxes( proton_axes );
  proton_3sigma.AddCut( "trNsigmaProton", [](double n_sigma){
    return n_sigma < 3;
  }, "proton_3sigma cut" );
  proton_3sigma.AddCut( "trFhcalX", [](double pos){
    return pos < -30.0 || pos > 160;
  }, "cut on x-pos in fhcal plane" );
  proton_3sigma.AddCut( "trFhcalY", [](double pos){
    return pos < -60.0 || pos > 60;
  }, "cut on y-pos in fhcal plane" );
  proton_3sigma.AddCut( "trStsNhits", [](double nhits){
    return nhits > 5.5;
  }, "cut on fake tracks" );
  proton_3sigma.AddCut( "trDcaR", [](double dca){
    return dca < 5.0;
  }, "DCA cut" );
  proton_3sigma.AddCut( "trStsChi2", [](double chi2){
    return chi2 < 5.0;
  }, "cut on chi2 in sts" );
  proton_3sigma.AddHisto2D({{"trProtonY", 100, -0.5, 1.5}, {"trPt", 100, 0.0, 2.0}});
  correction_task.AddVector(proton_3sigma);

  std::cout << "Initialized" << std::endl;

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}