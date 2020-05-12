#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "snp_filler_001_bed.h"
#include "snp_filler_011_bed.h"
#include "gwas_logit_offset.h"

template<typename filler> 
inline List GWAS_logit_offset_bed(NumericVector Y, NumericVector Offset, NumericMatrix Q, 
                                  double tol, int max_iter, filler S) {
  gwas_logit_offset<typename filler::type> x(Y, Offset, Q, (typename filler::type) tol, max_iter, S);
  x.run_tests();
  return S.L;
}


//[[Rcpp::export]]
List GWAS_logit_offset_bed(XPtr<matrix4> pA, NumericVector p, NumericVector Y, NumericVector Offset,
                           NumericMatrix Q, int beg, int end, double tol, int max_iter, std::string coding) {
  if(coding == "012") { 
    snp_filler_additive_bed<float> S(pA, p, beg, end);
    return GWAS_logit_offset_bed(Y, Offset, Q, tol, max_iter, S);
  } else if(coding == "011") {
    snp_filler_011_bed<float> S(pA, p, beg, end);
    return GWAS_logit_offset_bed(Y, Offset, Q, tol, max_iter, S);
  } else if(coding == "001") {
    snp_filler_001_bed<float> S(pA, p, beg, end);
    return GWAS_logit_offset_bed(Y, Offset, Q, tol, max_iter, S);
  } else {
    stop("Unknown coding value");
  }
}

