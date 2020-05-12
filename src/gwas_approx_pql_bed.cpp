#include <Rcpp.h>
#include "snp_filler_bed.h"
#include "snp_filler_001_bed.h"
#include "snp_filler_011_bed.h"
#include "gwas_approx_pql.h"


template<typename filler>
inline List GWAS_approx_pql_bed(NumericVector PY, NumericMatrix P, filler S) {
  gwas_approx_pql<typename filler::type> x(PY, P, S);
  x.run_tests();
  return S.L;
}

//[[Rcpp::export]]
List GWAS_approx_pql_bed(XPtr<matrix4> pA, NumericVector PY, NumericMatrix P, NumericVector p, int beg, int end, std::string coding) {
  if(coding == "012") {
    snp_filler_additive_bed<float> S(pA, p, beg, end);
    return GWAS_approx_pql_bed(PY, P, S);
  } else if(coding == "011") {
    snp_filler_011_bed<float> S(pA, p, beg, end);
    return GWAS_approx_pql_bed(PY, P, S);
  } else if(coding == "001") {
    snp_filler_001_bed<float> S(pA, p, beg, end);
    return GWAS_approx_pql_bed(PY, P, S);
  } else {
    stop("Unknown coding value");
  }
}

