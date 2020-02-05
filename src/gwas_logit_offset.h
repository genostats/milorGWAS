#include <Rcpp.h>
#include <RcppEigen.h>
#include "snp_filler.h"
#include "logit_model.h"

#ifndef GWAS_LMM_SCORE
#define GWAS_LMM_SCORE

template<typename scalar_t>
using MATRIX = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template<typename scalar_t>
using VECTOR = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

template<typename scalar_t>
class gwas_logit_offset {
  public:
  int n, r;
  VECTOR<scalar_t> y;
  VECTOR<scalar_t> offset;
  MATRIX<scalar_t> q;

  VECTOR<scalar_t> SNP;
  snp_filler<scalar_t> & S;

  scalar_t tol;
  int max_iter;  


  // Q matrice avec Q'Q = Id 
  // [Issue de la décomposition QR de la matrice de covariables]
  gwas_logit_offset(NumericVector Y, NumericVector Offset, NumericMatrix Q, scalar_t tol_, int max_iter_, snp_filler<scalar_t> & S_) 
  : n(Y.size()), r(Q.ncol()), tol(tol_), max_iter(max_iter_), y(n), offset(n), q(n,r), SNP(n), S(S_) {
    if(Q.nrow() != n || Offset.size() != n) 
      stop("Dimensions mismatch\n");

    // copy has to be done when scalar_t = float
    for(int i = 0; i < n; i++) {
      y(i) = (scalar_t) Y[i];
      offset(i) = (scalar_t) Offset[i];
    }
    for(int i = 0; i < n; i++)
      for(int j = 0; j < r; j++)
        q(i,j) = (scalar_t) Q(i,j);
  }

  void run_tests() {
    std::vector<scalar_t> BETA;
    std::vector<scalar_t> SDBETA;

    VECTOR<scalar_t> beta(1);
    MATRIX<scalar_t> varbeta(1,1);

    while( S.snp_fill( &SNP(0,0) ) ) {
      if( S.current_snp_monomorphic() ) {
        BETA.push_back(NAN);
        SDBETA.push_back(NAN);
        continue;
      }
      // prendre le résidu de SNP par la régression sur Q
      MATRIX<scalar_t> res_SNP( SNP - q * (q.transpose() * SNP) );

      // régression
      logistic_model_offset<scalar_t>(y, offset, res_SNP, beta, varbeta, tol, max_iter);
      BETA.push_back( beta(0) );
      // correction de variance (inutile)
      // SDBETA.push_back( sqrt( varbeta(0,0)*(n-1)/(n-r-1) ) );
      SDBETA.push_back( sqrt( varbeta(0,0) ) );  
    }
 
    S.L["beta"] = wrap(BETA);
    S.L["sd"] = wrap(SDBETA);
  }
};

#endif
