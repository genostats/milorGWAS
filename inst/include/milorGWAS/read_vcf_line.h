#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "gaston/gzstream.h"
#include "token.h"
#include "stringstream_lite.h"
#include "chr_convert.h"
#include "snp_filter.h"
#include "gaston/flip_strand.h"
#include "genotype_probas.h"

#ifndef GASTONread_vcf_line
#define GASTONread_vcf_line
using namespace Rcpp;

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT

inline void parse_vcf_line(std::string line, std::string & snp_id, int & snp_pos, 
                            std::string & chr, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string qual, filter, info, format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format))
    stop("VCF file format error");
}


template<typename scalar_t> 
void parse_vcf_line_dosages(std::string line, std::vector<scalar_t> & dosage, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  std::istringstream li(line);
  std::string qual, filter, info, format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format))
    stop("VCF file format error");
  
  int pos = token_position(format, "DS");
  if(pos != -1) {
    std::string G;
    while(li >> G) {
      dosage.push_back( token_at_position<scalar_t>(G, pos) );
    }
    return;
  }

  pos = token_position(format, "GP");
  if(pos != -1) {
    std::string G;
    while(li >> G) {
      std::string GP( token_at_position<std::string>(G, pos) );
      dosage.push_back( genotype_probas_to_dosage(GP) );
    } 
    return;
  }
  stop("No DS / GP field");
}

template<typename scalar_t> 
bool read_vcf_line_dosages(igzstream & in, std::vector<scalar_t> & dosage, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2) {
  std::string line;
  if(!std::getline(in, line))
    return false;
  parse_vcf_line_dosages<scalar_t>(line, dosage, snp_id, snp_pos, chr, A1, A2);
  return true;
}




template<typename scalar_t>
inline scalar_t geno_conv(char * s, int le) {
  scalar_t g = 0;
  if(le == 3) { // cas diploide
    if(s[0] == '1') g++;
    if(s[2] == '1') g++;
    if(s[0] == '.' || s[2] == '.') g = 3; // missing value : NA
  } else if(le == 1) { // cas haploide
    if(s[0] == '1') g++;
    if(s[0] == '.') g = 3; // missing value : NA
  } else {
    g = 3;
  }
  return g;
}

// avec chr = string
// conservé au cas où 
template<typename scalar_t>
void parse_vcf_line_genotypes(std::string line, std::vector<scalar_t> & genotypes, std::string & snp_id,
                     int & snp_pos, std::string & chr, std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format;
  if(!(li >> chr >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  
  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    // conversion du token t1 en génotype
    scalar_t g = geno_conv<scalar_t>(tok.token, tok.token_length);
    genotypes.push_back(g);
  }
}

// idem avec chr = int
template<typename scalar_t>
void parse_vcf_line_genotypes(std::string line, std::vector<scalar_t> & genotypes, std::string & snp_id,
                     int & snp_pos, int & chr, std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    // conversion du token t1 en génotype
    scalar_t g = geno_conv<scalar_t>(tok.token, tok.token_length);
    genotypes.push_back(g);
  }
}

// version 'filtered'
template<typename scalar_t>
bool parse_vcf_line_genotypes_filtered(std::string line, std::vector<scalar_t> & genotypes, 
                     std::string & snp_id, int & snp_pos, int & chr, 
                     std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info, snp_filter & FILTER) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  bool swap = false, flip = false;
  if(!FILTER(snp_id, chr, snp_pos, A1, A2, flip, swap)) 
    return false;

  if(swap) {
    std::string tmp(A1);
    A1 = A2;
    A2 = tmp;
  }
  if(flip) {
    A1 = flip_strand(A1);
    A2 = flip_strand(A2);
  }

  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    // conversion du token t1 en génotype
    scalar_t g = geno_conv<scalar_t>(tok.token, tok.token_length);
    if(swap)
      genotypes.push_back(2-g);
    else
      genotypes.push_back(g);
  }
  return true;
}


// version 'filtered' + vecteur de booleens pour déterminer quels individus on prend
template<typename scalar_t>
bool parse_vcf_line_genotypes_filtered(std::string line, std::vector<scalar_t> & genotypes, 
                     std::string & snp_id, int & snp_pos, int & chr, 
                     std::string & A1, std::string & A2, double & qual,
                     std::string & filter, std::string & info, snp_filter & FILTER, 
                     std::vector<bool> & which_samples) {

  stringstream_lite li(line, 9); // 9 = tab separated
  std::string format, chr_;
  if(!(li >> chr_ >> snp_pos >> snp_id >> A1 >> A2 >> qual >> filter >> info >> format)) {
    stop("VCF file format error");
  }
  chr = chr_to_int(chr_);

  int pos = token_position(format, "GT");
  if(pos < 0) stop("VCF error (No 'GT' found)");

  bool swap = false, flip = false;
  if(!FILTER(snp_id, chr, snp_pos, A1, A2, flip, swap)) 
    return false;

  if(swap) {
    std::string tmp(A1);
    A1 = A2;
    A2 = tmp;
  }
  if(flip) {
    A1 = flip_strand(A1);
    A2 = flip_strand(A2);
  }

  int k = 0;
  while( li.next_token() > 0 ) { // li.token pointe sur une chaîne avec le génotype en position pos
    stringstream_lite tok(li.token, ':'); // les champs sont séparés par des ':'
    for(int i = 0; i <= pos; i++) { // <= pos car même si pos = 0 il faut lire un token... 
      if(tok.next_token() == 0) 
        stop("VCF file format error");
    }

    if(which_samples[k++]) {
      // conversion du token t1 en génotype
      scalar_t g = geno_conv<scalar_t>(tok.token, tok.token_length);
      if(swap)
        genotypes.push_back(2-g);
      else
        genotypes.push_back(g);
    }
  }
  return true;
}

#endif
