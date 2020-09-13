#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "dosage_files.h"
#include "gaston/gzstream.h"
#include "token.h"
#include "read_gen_line.h"
#include "read_vcf_header.h"
#include "read_vcf_line.h"

using namespace Rcpp;

dosages::dosages(std::string file) : filename(file), in( (char *) &filename[0u] ) {
  start();
}

dosages::dosages(const char * file) : filename(file), in( file ) {
  start();
}

dosages::dosages(const CharacterVector Filename) : filename(Filename[0]), in( (const char *) &filename[0u] ) {
  start();
}

dosages::~dosages() {
  in.close();
}

// se met en début de fichier 
// initialise le vecteur samples quand c'est un VCF
void dosages::start() {
  if(!in.good()) 
    stop("Can't open file");
  if(!std::getline(in, line)) 
    stop("File is empty");

  if(line.substr(0,1) == "#") { // VCF
    // Rcout << "VCF\n";
    type = VCF;
    // skip description informations
    while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
      if(line.substr(0,2) != "##") {
        read_vcf_samples(line, samples);
        break; // fin
      }
    }
    // read first data line
    if(std::getline(in, line)) 
      good = true;
    else
      good = false;
    return;  // we are done
  }
  if(line.substr(0,3) == "---") { // Impute2 (sauce Agnès ? Ou vrai fichier de sortie de Impute2 ?)
    // Rcout << "Impute2\n";
    type = Impute2;
    good = true;
    return;
  }
  std::istringstream li(line);
  std::vector<std::string> splitted;
  std::string str;
  while(li >> str) splitted.push_back(str);
  if(splitted.size() < 5 || (splitted[0] != "id"  && splitted[0] != "ID" )
                         || (splitted[1] != "chr" && splitted[1] != "CHR") 
                         || (splitted[2] != "pos" && splitted[2] != "POS")
                         || (splitted[3] != "A1"  && splitted[3] != "a1" )
                         || (splitted[4] != "A2"  && splitted[4] != "a2" ) ) {
    in.close();
    stop("Unknown file format");
  }
  type = PES;
  // set sample names
  for(int i = 5; i < splitted.size(); i++)
    samples.push_back( splitted[i] );
  // read first data line
  if(std::getline(in, line)) 
    good = true;
  else
    good = false;
}

bool dosages::read_line(std::vector<double> & dosage, std::string & snp_id, 
                   int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  if(!good) return false;
  if(type == Impute2) {
    chr = "NA";
    parse_gen_line<double>(line, dosage, snp_id, snp_pos, A1, A2);
  }
  if(type == VCF)
    parse_vcf_line_dosages<double>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  if(type == PES)
    parse_gen_line_pes<double>(line, dosage, snp_id, chr, snp_pos, A1, A2);
  // read next line now...
  if(std::getline(in, line))
    good = true;
  else
    good = false;

  return true;
}

bool dosages::read_line(std::vector<float> & dosage, std::string & snp_id, 
                   int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  if(!good) return false;
  if(type == Impute2) {
    chr = "NA";
    parse_gen_line<float>(line, dosage, snp_id, snp_pos, A1, A2);
  }
  if(type == VCF)
    parse_vcf_line_dosages<float>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  if(type == PES)
    parse_gen_line_pes<float>(line, dosage, snp_id, chr, snp_pos, A1, A2);
  // read next line now...a
  if(std::getline(in, line))
    good = true;
  else
    good = false;

  return true;
}


bool dosages::read_line(std::string & snp_id, int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  if(!good) return false;
  if(type == Impute2) {
    chr = "NA";
    parse_gen_line(line, snp_id, snp_pos, A1, A2);
  }
  if(type == VCF)
    parse_vcf_line(line, snp_id, snp_pos, chr, A1, A2);
  if(type == PES)
    parse_gen_line_pes(line, snp_id, chr, snp_pos, A1, A2);
  // read next line now...a
  if(std::getline(in, line))
    good = true;
  else
    good = false;

  return true;
}


