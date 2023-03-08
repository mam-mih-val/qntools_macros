#include "utils.h"
#include "makeQvectors.h"

filteredDF defineVariables(definedDF &d);
void setupQvectors();

void makeQvectors(std::string inputFiles="/home/ogolosov/desktop/bman/data/run8/sim/dcm_4gev.root", std::string calibFilePath="qa.root", std::string outFilePath="qn.root")
{
  TStopwatch timer;
  timer.Start();
  std::string treename = "t";
  ROOT::RDataFrame d(*makeChain(inputFiles, treename.c_str()));
  auto dd=defineVariables(d);
  init(dd, outFilePath, calibFilePath);
  setupQvectors(); 
  timer.Stop();
  timer.Print();
  timer.Start();
  run(dd);
  std::cout << "Done!\n";
  timer.Stop();
  timer.Print();
}

filteredDF defineVariables(definedDF &d)
{
  auto dd=d
    .Define("fhcalModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:fhcalModPos) phi.push_back(pos.phi()); return phi;")
    .Define("fhcalModX","ROOT::VecOps::RVec<float> x; for(auto& pos:fhcalModPos) x.push_back(pos.x()); return x;")
    .Define("fhcalModY","ROOT::VecOps::RVec<float> y; for(auto& pos:fhcalModPos) y.push_back(pos.y()); return y;")
    .Define("fhcalModInSub1","fhcalModId>=0&&fhcalModId<=16")
    .Define("fhcalModInSub2","fhcalModId>=17&&fhcalModId<=34")
    .Define("fhcalModInSub3","fhcalModId>=35&&fhcalModId<=54")
    .Define("scwallModPhi","ROOT::VecOps::RVec<float> phi; for(auto& pos:scwallModPos) phi.push_back(pos.phi()); return phi;")
    .Define("trPt","ROOT::VecOps::RVec<float> pt; for(auto& mom:trMom) pt.push_back(mom.pt()); return pt;")
    .Define("trEta","ROOT::VecOps::RVec<float> eta; for(auto& mom:trMom) eta.push_back(mom.eta()); return eta;")
    .Define("trPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:trMom) phi.push_back(mom.phi()); return phi;")
    .Define("simPt","ROOT::VecOps::RVec<float> pt;for(auto& mom:simMom) pt.push_back(mom.pt()); return pt;")
    .Define("simEta","ROOT::VecOps::RVec<float> eta;for(auto& mom:simMom) eta.push_back(mom.eta()); return eta;")
    .Define("simPhi","ROOT::VecOps::RVec<float> phi;for(auto& mom:simMom) phi.push_back(mom.phi()); return phi;")
    .Filter("vtxChi2>0.0001") // at least one filter is mandatory!!!
  ;
  
  varPatterns=
  {
    "b",                                             // kEvent
    "(fhcal|scwall)Mod(X|Y|Phi|E|Id|InSub.)",        // kChannel 
    "tr(Pt|Eta|Phi|BetaTof400|BetaTof700|SimIndex)", // kRecParticle  
    "sim(Pt|Eta|Phi|Pdg|MotherId)"                   // kSimParticle  
  };

  return dd; 
}

void setupQvectors()
{
  std::vector<Qn::AxisD> corrAxesEvent=
  {
    {"b", 4,0,10},
  };

  std::vector<Qn::AxisD> corrAxesParticle=
  {
    {"trPt",4,0,2},
    {"trEta",4,0,2},
  };

  for (auto &axis:corrAxesEvent)
    man.AddCorrectionAxis(axis);

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(false);
  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(true);
  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
  
  auto sumW=Qn::QVector::Normalization::M;
  auto track=Qn::DetectorType::TRACK;
  auto channel=Qn::DetectorType::CHANNEL;
  auto plain=Qn::QVector::CorrectionStep::PLAIN;
  auto recentered=Qn::QVector::CorrectionStep::RECENTERED;
  auto twisted=Qn::QVector::CorrectionStep::TWIST;
  auto rescaled=Qn::QVector::CorrectionStep::RESCALED;
  
  man.AddDetector("fhcal1", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal1", recentering);
  man.AddCorrectionOnQnVector("fhcal1", twistRescale);
  man.SetOutputQVectors("fhcal1", {plain, recentered});
  man.AddCutOnDetector("fhcal1", {"fhcalModInSub1"}, equal(1), "fhcal1");
  man.AddHisto2D("fhcal1", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub1");
  man.AddHisto2D("fhcal1", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub1");
  
  man.AddDetector("fhcal2", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal2", recentering);
  man.AddCorrectionOnQnVector("fhcal2", twistRescale);
  man.SetOutputQVectors("fhcal2", {plain, recentered});
  man.AddCutOnDetector("fhcal2", {"fhcalModInSub2"}, equal(1), "fhcal2");
  man.AddHisto2D("fhcal2", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub2");
  man.AddHisto2D("fhcal2", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub2");
  
  man.AddDetector("fhcal3", channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCorrectionOnQnVector("fhcal3", recentering);
  man.AddCorrectionOnQnVector("fhcal3", twistRescale);
  man.SetOutputQVectors("fhcal3", {plain, recentered});
  man.AddCutOnDetector("fhcal3", {"fhcalModInSub3"}, equal(1), "fhcal3");
  man.AddHisto2D("fhcal3", {{"fhcalModId", 100, 0., 100}, {"fhcalModE", 100, 0., 10}}, "fhcalModInSub3");
  man.AddHisto2D("fhcal3", {{"fhcalModX", 100, -100, 100}, {"fhcalModY", 100, -100, 100}}, "fhcalModInSub3");

  man.AddDetector("tr", track, "trPhi", "Ones", corrAxesParticle, {1,2}, sumW);
  man.AddCutOnDetector("tr", {"particleType"}, equal(kRecParticle), "recParticle");
//  correction_manager_.AddCutOnDetector("tr", {"simPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCorrectionOnQnVector("tr", recentering);
  man.AddCorrectionOnQnVector("tr", twistRescale);
  man.SetOutputQVectors("tr", {plain, recentered, twisted, rescaled});
  man.AddHisto1D("tr", {"trPhi", 100, -3.15, 3.15}, "Ones");
//  correction_manager_.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");
}

int main( int n_args, char** args ){
  if( n_args < 2 )
    throw std::runtime_error(std::string( "Too few arguments were provided. At least 1 is expected." ));
  std::string in_file = args[1];

  makeQvectors(in_file);
}