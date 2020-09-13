#include <Rcpp.h>
#include "gaston/snp_hash.h"

using namespace Rcpp;

#ifndef SNP_FILTER
#define SNP_FILTER
enum filter_type {nofilter, chr_filter, range_bp, range_cm, hash};

class snp_filter {
  SNPhash H;
  int chr_, low_bp, high_bp;
  double low_cm, high_cm;

  filter_type t;

  public:

  snp_filter() :
    t(nofilter) {}

  snp_filter(int chr) : chr_(chr), t(chr_filter) {}

  snp_filter(int chr, int low, int high) : 
    chr_(chr), low_bp(low), high_bp(high), t(range_bp) {}

  snp_filter(int chr, double low, double high) :
    chr_(chr), low_cm(low), high_cm(high), t(range_cm) {}

  snp_filter(IntegerVector CHR, IntegerVector POS) : 
    H(CHR, POS), t(hash) {}

  snp_filter(IntegerVector CHR, IntegerVector POS, CharacterVector A1, CharacterVector A2) : 
    H(CHR, POS, A1, A2), t(hash) {}

  snp_filter(CharacterVector ID) : 
    H(ID), t(hash) {}

  snp_filter(CharacterVector ID, IntegerVector CHR, IntegerVector POS) : 
    H(ID, CHR, POS), t(hash) {}

  snp_filter(CharacterVector ID, IntegerVector CHR, IntegerVector POS, CharacterVector A1, CharacterVector A2) : 
    H(ID, CHR, POS, A1, A2), t(hash) {}

  bool operator()(const std::string & snp, int chr, int bp);
  bool operator()(const std::string & snp, int chr, int bp, const std::string & A1, const std::string & A2, bool & flip, bool & swap);
  bool operator()(const std::string & snp, int chr, int bp, double cm);
  bool operator()(int chr, int bp);
  bool operator()(const std::string & snp);
};

#endif
