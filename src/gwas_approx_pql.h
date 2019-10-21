#include <Rcpp.h>
#include <RcppEigen.h>
#include "snp_filler.h"

#ifndef GWAS_LMM_SCORE
#define GWAS_LMM_SCORE

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_approx_pql {
  public:
  int n;
  VECTOR<scalar_t> Py;
  MATRIX<scalar_t> PP;
  VECTOR<scalar_t> SNP;
  snp_filler<scalar_t> & S;

  gwas_approx_pql(NumericVector PY, NumericMatrix P, snp_filler<scalar_t> & S_) 
  : n(PY.size()), Py(n), PP(n,n), SNP(n), S(S_) {
    if(P.nrow() != n || P.ncol() != n) 
      stop("Dimensions mismatch\n");

    // copy has to be done when scalar_t = float
    for(int i = 0; i < n; i++) Py(i) = PY[i];

    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) 
        PP(i,j) = P(i,j);
  }

  void run_tests() {
    std::vector<double> beta;
    std::vector<double> sd_beta;
    double t, v;  
    while( S.snp_fill( &SNP[0] ) ) {
      if( S.current_snp_monomorphic() ) {
        beta.push_back(NAN);
        sd_beta.push_back(NAN);
        continue;
      }
      v = (PP.template selfadjointView<Eigen::Lower>() * SNP).dot(SNP); // GPG
      t = SNP.dot(Py); // GPZ
      beta.push_back(t/v);
      sd_beta.push_back(1/sqrt(v));
    }
 
    S.L["beta"] = wrap(beta);
    S.L["sd"] = wrap(sd_beta);
  }
};

#endif
