#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "token.h"
#include "read_gen_line.h"
#include "read_vcf_header.h"
#include "read_vcf_line.h"

#ifndef DOSAGES
#define DOSAGES

enum dosage_type {VCF, Impute2, PES};

using namespace Rcpp;

class dosages {
public:
  std::string filename;
  igzstream in;
  std::string line;
  dosage_type type;
  bool good;
  std::vector<std::string> samples; //sera rempli pour les VCF

  dosages(std::string file);
  dosages(const char * file);
  dosages(const CharacterVector Filename);
  ~dosages();


private:  
  void start();

public:
  bool read_line(std::vector<double> & dosage, std::string & snp_id, 
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

  bool read_line(std::vector<float> & dosage, std::string & snp_id, 
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

  bool read_line(std::string & snp_id, int & snp_pos, std::string & chr, std::string & A1, std::string & A2);

};

#endif


