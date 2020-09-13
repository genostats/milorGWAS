#include <Rcpp.h>

#ifndef CHR_CONVERT
#define CHR_CONVERT
using namespace Rcpp;

void set_chr_ids(List L);

int chr_to_int(std::string & chr);

inline int chr_to_int(char * chr);

#endif
