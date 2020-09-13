#include <Rcpp.h>
#include "gaston/matrix4.h"
#include "snp_filler.h"

#ifndef GASTON_SNP_FILLER_011_BED
#define GASTON_SNP_FILLER_011_BED

using namespace Rcpp;

template<typename scalar_t>
class snp_filler_011_bed : public snp_filler<scalar_t> {
  public:
  XPtr<matrix4> pA;  
  int ncol, true_ncol;
  NumericVector p; 
  int beg, end;
  int i;
  snp_filler_011_bed(XPtr<matrix4> pA_, NumericVector p_, int beg_, int end_)
    : snp_filler<scalar_t>(), pA(pA_), ncol(pA->ncol), true_ncol(pA->true_ncol), 
      p(p_), beg(beg_), end(end_), i(beg) { };

  bool snp_fill(scalar_t * SNP) {
    if(i > end) {
      this->monomorphic = true; 
      return false; 
    }
    if( std::isnan(p(i)) || p(i) == 0 || p(i) == 1 ) {
      this->monomorphic = true;
      i++;
      return true;
    }
    uint8_t * snp = pA-> data[i];
    scalar_t gg[4] = { 0.0, 1.0, 1.0, (scalar_t) p(i)*(2 - (scalar_t) p(i)) };
    // p(i)*(2-p(i)) est E(X) si X codé 011
    for(int ii = 0; ii < true_ncol-1; ii++) {
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4; ss++) {
        SNP[4*ii+ss] = gg[ x&3 ];
        x >>= 2;
      }
    }
    { int ii = true_ncol-1;
      uint8_t x = snp[ii];
      for(int ss = 0; ss < 4 && 4*ii+ss < ncol; ss++) {
        SNP[4*ii+ss] = gg[ x&3 ];
        x >>= 2;
      }
    }
    i++;
    this->monomorphic = false;
    return true;
  }
};

#endif
