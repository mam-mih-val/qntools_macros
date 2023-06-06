//
// Created by Misha on 3/7/2023.
//

#include <TROOT.h>
#include <stdexcept>
#include <iostream>

int main(int n_args, char** args){
  if( n_args < 2 )
    throw std::runtime_error( "No argumets provided" );
  std::string macro{args[1]};
  std::string macro_full = macro+"(";
  for( int i=2; i<n_args; ++i ){
    macro_full+="\""+ std::string{args[i]}+"\",";
  }
  macro_full.back()=')';
  std::cout << "Executing "+macro_full << std::endl;
  gROOT->Macro(macro_full.c_str());
  return 0;
}