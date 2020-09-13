#include "snp_filter.h"
#include <string>

using namespace Rcpp;

// ------------------------------------------------------------------
bool snp_filter::operator()(const std::string & snp, int chr, int bp) {

  if(t == nofilter || t == range_cm) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp)
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
 

  // if(!H.m) return true; // empty hash

  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  if(H.htype == snpid) {
    int a = H.lookup(snp);
    return (a != NA_INTEGER);
  }

  if(H.htype == snpid_chr_pos) {
    int a = H.lookup(snp, chr, bp);
    return (a != NA_INTEGER);
  }

  stop("Wrong hash type !");
  return false; // ???
}

// ------------------------------------------------------------------
bool snp_filter::operator()(const std::string & snp, int chr, int bp, double cm) {

  if(t == nofilter) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp)
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
 
  if(t == range_cm) 
    return( chr == chr_ && low_cm <= cm && cm <= high_cm );
 

  // if(!H.m) return true; // empty hash

  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  if(H.htype == snpid) {
    int a = H.lookup(snp);
    return (a != NA_INTEGER);
  }

  if(H.htype == snpid_chr_pos) {
    int a = H.lookup(snp, chr, bp);
    return (a != NA_INTEGER);
  }

  stop("Wrong hash type !");
  return false; // ???
}

// ------------------------------------------------------------------
bool snp_filter::operator()(int chr, int bp) {
  if(t == nofilter) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp) 
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
  
  if(t != hash) return true;

  int a = H.lookup(chr, bp);
  return (a != NA_INTEGER);
}

// ------------------------------------------------------------------
// le cas le plus complet... (pour le filtrage des VCF)
// il faudra compléter pour les divers types de hash... on n'a que snpid, chr_pos et chr_pos_al pour le moment
bool snp_filter::operator()(const std::string & id, int chr, int bp, const std::string & A1, const std::string & A2, bool & flip, bool & swap) {

  swap = false; // le défaut !
  flip = false;

  if(t == nofilter) return true;

  if(t == chr_filter) 
    return (chr == chr_);

  if(t == range_bp) 
    return( chr == chr_ && low_bp <= bp && bp <= high_bp );
  
  if(t != hash) return true;

  if(H.htype == snpid) {
    int a = H.lookup(id);
    return (a != NA_INTEGER);
  }

  if(H.htype == chr_pos) {
    int a = H.lookup(chr, bp);
    return (a != NA_INTEGER);
  }
  
  if(H.htype == chr_pos_al) {
    int a = H.lookup(chr, bp, A1, A2, flip, swap);
    return (a != NA_INTEGER);
  }

  stop("Wrong hash type !");
  return false; // ???
}


// ------------------------------------------------------------------
bool snp_filter::operator()(const std::string & snp) {
 
  if(t != hash) return true;

  int a = H.lookup(snp);
  return (a != NA_INTEGER);
}

