#include <Rcpp.h>
#include "snp_filler_dosages.h"
#include "gwas_approx_pql.h"

// [[Rcpp::export]]
List GWAS_approx_pql_dosages(CharacterVector filename, NumericVector PY, NumericMatrix P, int beg, int end, double tol) {
  snp_filler_dosages<float> S(filename, beg, end, PY.size());
  gwas_approx_pql<float> x(PY, P, S);
  x.run_tests();

  List R;
  R["id"] = wrap(S.SNP_ID);
  R["chr"] = wrap(S.CHR);
  R["pos"] = wrap(S.POS);
  R["A1"] = wrap(S.AL1);
  R["A2"] = wrap(S.AL2);
  R["freq1"] = wrap(S.F1);
  R["freq2"] = wrap(S.F2);
  R["beta"] = S.L["beta"];
  R["sd"] = S.L["sd"];
  
  return R;
}



