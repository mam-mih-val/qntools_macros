//
// Created by Misha on 3/7/2023.
//

#define QNTOOLS_INCLUDE_PATH "/home/mikhail/QnTools/install/lib/cmake/QnTools/../../../include/QnTools"
#define THIS_SRC_PATH "src"

#include <TROOT.h>
#include <stdexcept>
#include <TSystem.h>
#include <TSystem.h>

int main(int n_args, char** args){
  if( n_args < 2 )
    throw std::runtime_error( "No argumets provided" );
  gInterpreter->AddIncludePath(QNTOOLS_INCLUDE_PATH);
  gInterpreter->AddIncludePath(THIS_SRC_PATH);
  std::string macro{args[1]};
  std::string macro_full = macro+"(";
  for( int i=2; i<n_args; ++i ){
    macro_full+="\""+ std::string{args[i]}+"\",";
  }
  macro_full.back()=')';
  gROOT->Macro(macro_full.c_str());
  return 0;
}
