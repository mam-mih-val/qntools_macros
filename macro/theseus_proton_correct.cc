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


void mcini_proton_correct(  std::string list, std::string str_sqrt_snn="2.4", std::string str_nucleus_mass="197",
                            std::string calib_in_file="qa.root"){
  
  const double sqrt_snn = std::stod(str_sqrt_snn);
  const double M = 0.938;
  const double T = sqrt_snn * sqrt_snn/ 2 / M - 2*M;
  const double GAMMA = (T + M) / M;
  const double BETA = sqrt(1 - (M * M) / (M + T) / (M + T));
  const double PZ = M * BETA * GAMMA;
  const double E = T + M;
  const double Y_BEAM = 0.5 * log((E + PZ) / (E - PZ)) / 2.0;
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );


  std::random_device dev{}; 
  std::mt19937 engine(dev()); 
  std::uniform_real_distribution<float> distribution{M_PI*-1, M_PI};
  auto psi_rp_function = [ &distribution, &engine ] (){
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
    // Автоматически извлекаем b из имени файла и нормируем на радиус ядра
          .Define( "b_norm", [chain, NUCLEUS_RADIUS](){
            // Получаем имя текущего открытого файла из TChain
            std::string filepath = chain->GetFile()->GetName();
            
            // Ищем подстроку вида "_Xfm_" (например, _2fm_, _6fm_, _10fm_)
            size_t pos_fm = filepath.find("fm_");
            double b_value = 0.0;
            
            if (pos_fm != std::string::npos) {
              size_t start_pos = filepath.rfind('_', pos_fm) + 1;
              std::string b_str = filepath.substr(start_pos, pos_fm - start_pos);
              b_value = std::stod(b_str); // Превращаем "2" или "6" в double
            } else {
              // На случай, если имя файла изменится, ставим дефолт 2.0, чтобы код не падал
              b_value = 2.0; 
            }
            return b_value / NUCLEUS_RADIUS;
          }, {} )
          
          // Истинная плоскость реакции в THESEUS всегда равна 0
          .Define( "reaction_plane", [](){ return 0.0; }, {} )
          // Rotated reaction plane:
          // .Define( "reaction_plane", psi_rp_function, {} )
         
          // Вычисляем phi напрямую (без искусственных сдвигов)
          .Define( "phi", []( ROOT::VecOps::RVec<double> vec_px, ROOT::VecOps::RVec<double> vec_py ){
            ROOT::VecOps::RVec<double> phi;
            for(size_t i=0; i<vec_px.size(); ++i) {
              phi.push_back(atan2(vec_py.at(i), vec_px.at(i)));
            }
            return phi;
          }, {"px", "py"} )
          // Rotated phi w.r.r reaction plane angle
          // .Define( "phi", phi_function, {"reaction_plane", "px", "py"} )
           // calculate pT
           .Define("pT", []( ROOT::VecOps::RVec<double> vec_px, ROOT::VecOps::RVec<double> vec_py ){
            ROOT::VecOps::RVec<double> vec_pT;
            for( size_t i=0; i<vec_px.size(); ++i ){
              auto px = vec_px.at(i);
              auto py = vec_py.at(i);
              vec_pT.push_back( sqrt(px*px + py*py) );
            }
            return vec_pT;
          },{ "px", "py" })
          
          // y_cm  =  rapidity!
          .Define("y_cm", []( ROOT::VecOps::RVec<double> vec_pz, ROOT::VecOps::RVec<double> vec_E ){
            ROOT::VecOps::RVec<double> vec_y;
            for( size_t i=0; i<vec_pz.size(); ++i ){
              auto pz = vec_pz.at(i);
              auto E = vec_E.at(i);
              auto y_cm = 0.5 * log( (E+pz)/(E-pz) );
              vec_y.push_back( y_cm );
            }
            return vec_y;
          },{ "pz", "E" })
          .Filter( "b_norm > 0." )
  ; // at least one filter is mandatory!!!

  auto correction_task = CorrectionTask( dd, "correction_out.root", calib_in_file );
  correction_task.SetEventVariables(std::regex("b_norm|reaction_plane"));
  correction_task.SetTrackVariables({
                                      std::regex("phi|pT|y_cm|id"),
                                    });

  correction_task.InitVariables();
  correction_task.AddEventAxis( { "b_norm", 16, 0, 16 } );

  VectorConfig psi_rp( "psi_rp", "reaction_plane", "Ones", VECTOR_TYPE::CHANNEL, NORMALIZATION::M );
  psi_rp.SetHarmonicArray( { 1, 2 } );
  psi_rp.SetCorrections( { CORRECTION::PLAIN } );
  correction_task.AddVector(psi_rp);

  std::vector<Qn::AxisD> tru_proton_axes{
        { "y_cm", 20, -0.5, 1.5 },
        { "pT", 10, 0.0, 2.0 },
  };

  VectorConfig tru_proton( "proton", "phi", "Ones", VECTOR_TYPE::TRACK, NORMALIZATION::M );
  tru_proton.SetHarmonicArray( { 1, 2 } );
  tru_proton.SetCorrections( {CORRECTION::PLAIN } );
  tru_proton.SetCorrectionAxes( tru_proton_axes );
  tru_proton.AddCut( "id", [](double pid){
    auto pdg_code = static_cast<int>(pid);
    return pdg_code == 2212;
    }, "proton cut" );
  correction_task.AddVector(tru_proton);

  correction_task.Run();
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
}