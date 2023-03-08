//
// Created by Misha on 3/7/2023.
//

#include <TROOT.h>
#include <stdexcept>

int main(int n_args, char** args){
  if( n_args < 2 )
    throw std::runtime_error( "No argumets provided" );
  gROOT->Macro(args[1]);
  return 0;
}