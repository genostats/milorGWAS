#include <Rcpp.h>
#include <iostream>
#include <string>
#include "token.h"
#include "stringstream_lite.h"

using namespace Rcpp;

template<>
double sto<double>(const std::string&x) {
  return std::stod(x);
}

template<>
float sto<float>(const std::string&x) {
  return std::stof(x);
}

template<>
int sto<int>(const std::string&x) {
  return std::stoi(x);
}

template<>
std::string sto<std::string>(const std::string&x) {
  return x;
}

int token_position(std::string s, std::string token) {
  std::istringstream ss(s);
  std::string tok;
  int k = 0;
  while( std::getline(ss, tok, ':') ) {
    if(tok == token) return k;
    k++;
  }
  return -1;
}

