#ifndef GENOMICSDB_JNI_EXCEPTION_H
#define GENOMICSDB_JNI_EXCEPTION_H

#include "headers.h"

class GenomicsDBJNIException : public std::exception {
  public:
    GenomicsDBJNIException(const std::string m="") : msg_("GenomicsDBJNIException : "+m) { ; } 
    ~GenomicsDBJNIException() { ; } 
    // ACCESSORS
    const char* what() const noexcept { return msg_.c_str(); }
  private:
    std::string msg_;
};

#endif
