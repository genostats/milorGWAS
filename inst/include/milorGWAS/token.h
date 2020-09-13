#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#ifndef TOKEN
#define TOKEN
using namespace Rcpp;

int token_position(std::string s, std::string token);

template<typename T>
T sto(const std::string & x);

template<typename T>
T token_at_position(std::string & s, int pos) {
  std::istringstream ss(s);
  std::string token;
  for(int i = 0; i < pos && std::getline(ss, token, ':'); i++) {}
  std::getline(ss, token, ':');
  T r = sto<T>(token);
  return r;
}

#endif
