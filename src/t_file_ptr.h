//
// Created by Misha on 5/9/2023.
//

#ifndef QNTOOLSINTERFACE_T_FILE_PTR_H
#define QNTOOLSINTERFACE_T_FILE_PTR_H


#include <TFile.h>

class TFilePtr {
public:
  TFilePtr() : file_(nullptr){}
  TFilePtr(const TFilePtr&) = delete;
  TFilePtr& operator=( const TFilePtr& ) = delete;
  TFilePtr(TFilePtr&& other)  noexcept {
    if( file_ )
      file_->Close();
    file_ = other.file_;
    other.file_ = nullptr;
  };
  TFilePtr& operator=( TFilePtr&& other ) noexcept {
    if( file_ )
      file_->Close();
    file_ = other.file_;
    other.file_ = nullptr;
    return *this;
  };
  explicit TFilePtr( const std::string& file_name ) : file_( TFile::Open( file_name.c_str() ) ){  }
  explicit TFilePtr( const std::string& file_name, const std::string& option ) : file_( TFile::Open( file_name.c_str(), option.c_str() ) ){  }
  ~TFilePtr(){ if(file_) file_->Close(); }
  explicit operator bool(){ return file_; }
  auto operator->(){ return file_; }
  auto Get(){ return file_; }
private:
  TFile* file_{nullptr};
};


#endif //QNTOOLSINTERFACE_T_FILE_PTR_H
