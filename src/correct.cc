//
// Created by Misha on 3/7/2023.
//

#include <TROOT.h>
#include <stdexcept>

int main(int n_args, char** args){
  if( n_args < 2 )
    throw std::runtime_error( "No argumets provided" );
  std::string macro{args[1]};
  std::string list{args[2]};
  std::string macro_full = macro+"(\""+list+"\")";
  gROOT->Macro(macro_full.c_str());
  return 0;
}