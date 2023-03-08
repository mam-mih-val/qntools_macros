#ifndef MAKEQVECTORS_H
#define MAKEQVECTORS_H

#include <iostream>
#include <string>
#include <vector>
#include <regex>

#include <TChain.h>
#include <ROOT/RDataFrame.hxx>
#include <QnTools/CorrectionManager.hpp>

using filteredDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
using definedDF=ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void>;
using ROOT::VecOps::RVec;

enum varType
{
  kEvent=0,
  kChannel,
  kRecParticle,
  kSimParticle,
  kNVarTypes
};

std::string qaFilePath;
Qn::CorrectionManager man; 
TFile *outFile;
TTree *outTree;
short particleTypePosition;

std::vector <std::string> varPatterns(kNVarTypes);
std::vector <std::vector<std::string>> varNames(kNVarTypes);
std::vector<short> vcInitialPositions(kNVarTypes, 0);

void InitVariables(filteredDF &d, std::vector <std::string> &patterns, std::vector<std::vector<std::string>> &varNames)
{
  // Add all needed variables
  short ivar=0;
  for (int type=0;type<kNVarTypes;type++)
  {
    std::vector<short> varSizes;
    for (auto &name : d.GetColumnNames()) 
    {
      short size{1};
      auto columnType=d.GetColumnType(name);
      if (type==kEvent && columnType.find("ROOT::VecOps::RVec")<columnType.size())
        continue;
      if (!regex_match(name, std::regex(patterns.at(type))))
        continue;
      if (type==kChannel)
        size=*d.Range(0,0).Define("n", name+".size()").Mean("n");
      varNames.at(type).push_back(name);
      varSizes.push_back(size);
    }

    vcInitialPositions.at(type)=ivar;
    std::cout << "\nPattern: " << patterns.at(type) << std::endl;
    for (int i=0;i<varNames.at(type).size();i++)
    {
      auto name=varNames.at(type).at(i);
      auto size=varSizes.at(i);
      man.AddVariable(name, ivar, size);
      printf("%s: %i+%i\n", name.c_str(), ivar, size);
      ivar+=size;
      if (type==kEvent)
        man.AddEventVariable(name);
    }
  }
  particleTypePosition=ivar;
  printf("%s: %i+%i\n", "particleType", particleTypePosition, 1);
  man.AddVariable("particleType", particleTypePosition, 1);
}

void DefineVariableFilling(filteredDF &d, std::vector <std::vector<std::string>> &varNames)
{
  std::string vcEventFillExpr="vector<float> vcEvent; ";
  for (auto& var:varNames.at(kEvent))
    vcEventFillExpr+=Form("vcEvent.push_back(%s); ",var.c_str());
  vcEventFillExpr+="return vcEvent;";
  std::cout << std::endl << vcEventFillExpr << std::endl;

  std::string vcModuleFillExpr="vector<float> vcModule; ";
  for (auto& var:varNames.at(kChannel))
    vcModuleFillExpr+=Form("for(auto& val:%s) vcModule.push_back(val); ",var.c_str());
  vcModuleFillExpr+="return vcModule;";
  std::cout << std::endl << vcModuleFillExpr << std::endl;

  std::string nRecPart="0";
  if (varNames.at(kRecParticle).size()>0)
    nRecPart=varNames.at(kRecParticle).front() + ".size()";
  std::string vcRecPartFillExpr=Form("vector<vector<float>> vcRecParticle(%s); for(int i=0;i<vcRecParticle.size();i++) { ", nRecPart.c_str());
  for (auto& var:varNames.at(kRecParticle))
    vcRecPartFillExpr+=Form("vcRecParticle.at(i).push_back(%s.at(i)); ", var.c_str());
  vcRecPartFillExpr+=" } return vcRecParticle;";
  std::cout << std::endl << vcRecPartFillExpr << std::endl;

  std::string nSimPart="0";
  if (varNames.at(kSimParticle).size()>0)
    nSimPart=varNames.at(kSimParticle).front() + ".size()";
  std::string vcSimPartFillExpr=Form("vector<vector<float>> vcSimParticle(%s); for(int i=0;i<vcSimParticle.size();i++) { ", nSimPart.c_str());
  for (auto& var:varNames.at(kSimParticle))
    vcSimPartFillExpr+=Form("vcSimParticle.at(i).push_back(%s.at(i)); ", var.c_str());
  vcSimPartFillExpr+=" } return vcSimParticle;";
  std::cout << std::endl << vcSimPartFillExpr << std::endl;
  
  d=d.Define("vcEvent", vcEventFillExpr)
     .Define("vcModule", vcModuleFillExpr)
     .Define("vcRecParticle", vcRecPartFillExpr)
     .Define("vcSimParticle", vcSimPartFillExpr)
  ;
//  d.Display({/*"vcEvent","vcModule","vcRecParticle",*/"vcSimParticle"},2)->Print();
}

void processEvent(const ULong64_t eventId, const std::vector<float> vcEvent, const std::vector<float> vcModule, const std::vector<std::vector<float>> vcRecParticle, const std::vector<std::vector<float>> vcSimParticle)
{
  std::cout << "\r" << eventId;
  man.Reset();
  double *vc = man.GetVariableContainer();
  for(int ivar=0; ivar<vcEvent.size();ivar++)
    vc[vcInitialPositions.at(kEvent)+ivar]=vcEvent.at(ivar); 
  for(int ivar=0; ivar<vcModule.size();ivar++)
    vc[vcInitialPositions.at(kChannel)+ivar]=vcModule.at(ivar); 
  man.ProcessEvent();
  man.FillChannelDetectors();
  
  vc[particleTypePosition]=kRecParticle; 
  for(int ipart=0;ipart<vcRecParticle.size();ipart++)
  {
    for(int ivar=0;ivar<vcRecParticle.at(0).size();ivar++)
      vc[vcInitialPositions.at(kRecParticle)+ivar]=vcRecParticle.at(ipart).at(ivar);
    man.FillTrackingDetectors();
  }
  
  vc[particleTypePosition]=kSimParticle; 
  for(int ipart=0;ipart<vcSimParticle.size();ipart++)
  {
    for(int ivar=0;ivar<vcSimParticle.at(0).size();ivar++)
      vc[vcInitialPositions.at(kSimParticle)+ivar]=vcSimParticle.at(ipart).at(ivar);
    man.FillTrackingDetectors();
  }
  man.ProcessCorrections();
}

void init(filteredDF &d, std::string outFilePath, std::string calibFilePath)
{
  outFile = TFile::Open(outFilePath.c_str(), "recreate");
  outFile->cd();
  outTree = new TTree("tree", "tree");
  //qa_file_path_=regex_replace(outFilePath, regex("qn\\.root"), "qa.root");
  qaFilePath=regex_replace(outFilePath, std::regex("qn\\.root"), "qa.root");
  man.SetCalibrationInputFileName(calibFilePath);
  man.SetFillOutputTree(true);
  man.SetFillCalibrationQA(true);
  man.SetFillValidationQA(true);
  man.ConnectOutputTree(outTree);
  
  InitVariables(d, varPatterns, varNames);
}

void run(filteredDF &d) {
  man.InitializeOnNode();
  man.SetCurrentRunName("test");
  DefineVariableFilling(d, varNames);
  d.Foreach(processEvent, {"rdfentry_", "vcEvent", "vcModule", "vcRecParticle", "vcSimParticle"});
  std::cout << std::endl;
  
  man.Finalize();
  outFile->cd();
  outTree->Write("tree");
  outFile->Close();

  TFile qaFile(qaFilePath.c_str(), "RECREATE");
  man.GetCorrectionQAList()->Write("CorrectionQAHistograms", TObject::kSingleKey);
  man.GetCorrectionList()->Write("CorrectionHistograms", TObject::kSingleKey);
  qaFile.Close();
}

#endif // MAKEQVECTORS_H