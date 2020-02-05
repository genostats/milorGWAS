#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_logit_offset.h"

//[[Rcpp::export]]
List GWAS_logit_offset_bed(XPtr<matrix4> pA, NumericVector p, NumericVector Y, NumericVector Offset,
                           NumericMatrix Q, int beg, int end, double tol, int max_iter) {
  snp_filler_additive_bed<float> S(pA, p, beg, end);
  gwas_logit_offset<float> x(Y, Offset, Q, (float) tol, max_iter, S);
  x.run_tests();
  return S.L;
}

