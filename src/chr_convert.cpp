#include <Rcpp.h>
#include "chr_convert.h"

using namespace Rcpp;

static List chr_ids;

int chr_to_int(std::string & chr) {
  int r = std::atoi(chr.c_str());
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

