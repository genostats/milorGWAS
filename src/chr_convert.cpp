#include <Rcpp.h>
#include "chr_convert.h"

using namespace Rcpp;

static List chr_ids;

int chr_to_int(std::string & chr) {
  int r = std::atoi(chr.c_str());
  int c = chr_ids.size();
  if( r == 0 && chr_ids.containsElementNamed(chr.c_str()) )
    r = chr_ids[chr];
  return r;
}

int chr_to_int(char * chr) {
  int r = std::atoi(chr);
  if( r == 0 && chr_ids.containsElementNamed(chr) )
    r = chr_ids[chr];
  return r;
}

void set_chr_ids(List L) {
  chr_ids = L;
}

RcppExport SEXP gg_set_chr_ids(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    set_chr_ids(L);
    return R_NilValue;
END_RCPP
}


RcppExport SEXP gg_get_chr_ids(SEXP LSEXP) {
BEGIN_RCPP
  return wrap(chr_ids);
END_RCPP
}
