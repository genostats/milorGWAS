#include <Rcpp.h>
#include "snp_filler_dosages.h"
#include "gwas_logit_offset.h"

// [[Rcpp::export]]
List GWAS_logit_offset_dosages(CharacterVector filename, NumericVector Y, NumericVector Offset, NumericMatrix Q, int beg, int end, double tol, int max_iter) {
  snp_filler_dosages<float> S(filename, beg, end, Y.size());
  gwas_logit_offset<float> x(Y, Offset, Q, (float) tol, max_iter, S);
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



