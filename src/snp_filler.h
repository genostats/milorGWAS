#include <Rcpp.h>
#ifndef GASTON_SNP_FILLER
#define GASTON_SNP_FILLER

template<typename scalar_t>
class snp_filler {
  public:
  bool monomorphic;
  Rcpp::List L; // pour éventuellement mettre des infos sur les SNPs lus

  snp_filler() : monomorphic(true) {};

  virtual bool snp_fill(scalar_t * SNP) = 0;
  virtual bool current_snp_monomorphic() {
    return monomorphic;
  }
};

#endif
