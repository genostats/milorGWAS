#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "gwas_approx_pql.h"

//[[Rcpp::export]]
List GWAS_approx_pql_bed(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end) {
  snp_filler_additive_bed<double> S(pA, p, beg, end);
  gwas_approx_pql<double> x(PY, P, S);
  x.run_tests();
  return S.L;
}

