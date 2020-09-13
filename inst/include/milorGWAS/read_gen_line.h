#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "token.h"
#include <algorithm>
#ifndef GASTONread_gen_line
#define GASTONread_gen_line
using namespace Rcpp;

template<typename scalar_t>
inline void parse_gen_line(std::string line, std::vector<scalar_t> & dosage, std::string & snp_id, 
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string tirets;
  if(!(li >> tirets && li >> snp_id && li >> snp_pos && li >> A1 && li >> A2)) 
    stop("gen file format error");
  scalar_t p0, p1, p2;
  while(li >> p0 && li >> p1 && li >> p2) {
    dosage.push_back(p1 + 2*p2);
  }
}

inline void parse_gen_line(std::string line, std::string & snp_id, int & snp_pos, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string tirets;
  if(!(li >> tirets && li >> snp_id && li >> snp_pos && li >> A1 && li >> A2)) 
    stop("gen file format error");
}

template<typename scalar_t>
inline bool read_gen_line(igzstream & in, std::vector<scalar_t> & dosage, std::string & snp_id, 
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::string line;
  if(!std::getline(in, line)) {
    return false;
  }
  parse_gen_line(line, dosage, snp_id, snp_pos, A1, A2);
  return true;
}



template<typename scalar_t>
inline void parse_gen_line_pes(std::string line, std::vector<scalar_t> & dosage, std::string & snp_id, std::string & snp_chr,
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  if(!(li >> snp_id >> snp_chr >> snp_pos >> A1 >> A2))
    stop("gen file format error");
  scalar_t dose;
  /* version simple, qui ne prend pas en compte la possibilité d'avoir des NA 
  while(li >> dose)
    dosage.push_back(dose);
  */
  scalar_t somme(0);
  int n(0);
  std::string tok;
  while(li >> tok) {
    if(tok == "NA") {
      dosage.push_back(-1.0);
    } else {
      dose = sto<scalar_t>(tok);
      dosage.push_back(dose);
      somme += dose;
      n++;
    }
  }
  std::replace(dosage.begin(), dosage.end(), (scalar_t) -1.0, somme/n);
}

inline void parse_gen_line_pes(std::string line, std::string & snp_id, std::string & snp_chr, 
                               int & snp_pos, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  if(!(li >> snp_id >> snp_chr >> snp_pos >> A1 >> A2))
    stop("gen file format error");
}

template<typename scalar_t>
inline bool read_gen_line_pes(igzstream & in, std::vector<scalar_t> & dosage, std::string & snp_id, std::string & snp_chr,
                          int & snp_pos, std::string & A1, std::string & A2) {
  std::string line;
  if(!std::getline(in, line)) {
    return false;
  }
  parse_gen_line_pes(line, dosage, snp_id, snp_chr, snp_pos, A1, A2);
  return true;
}
#endif
