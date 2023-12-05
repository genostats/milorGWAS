#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include "dosage_files.h"

using namespace Rcpp;

//[[Rcpp::export]]
List read_dose_file(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<std::string> CHR, SNP_ID, AL1, AL2;
  std::vector<int> POS;
  std::vector<double> dosage;
  int last_dos_le(0), nb_ind(-1);
  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    // check for right number of dosages
    int dos_le = dosage.size();
    if(nb_ind < 0) {
      nb_ind = dos_le - last_dos_le;
    } else if(nb_ind != dos_le - last_dos_le) {
      Rcerr << "While reading SNP #" << POS.size()+1 << " with id = " << snp_id << "\n";
      Rcerr << "Read " << dos_le - last_dos_le << " dosages, instead of " << nb_ind << " on previous line(s)\n";
      stop("File format error");
    }
    last_dos_le = dos_le;

    SNP_ID.push_back(snp_id);
    POS.push_back(snp_pos);
    CHR.push_back(chr);
    AL1.push_back(A1);
    AL2.push_back(A2);
  }
  List L;
  L["snp.id"] = wrap(SNP_ID);
  L["chr"] = wrap(CHR);
  L["pos"] = wrap(POS);
  L["A1"] = wrap(AL1);
  L["A2"] = wrap(AL2);
  if(in.samples.size()> 0) L["samples"] = wrap(in.samples); // VCF ou PES
  NumericVector dos = wrap(dosage);
  dos.attr("dim") = Dimension( dosage.size()/POS.size(), POS.size() );
  L["dosages"] = dos;
  return L;
}

//[[Rcpp::export]]
NumericVector dose_file_dim(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  in.read_line(dosage, snp_id, snp_pos, chr, A1, A2);
  int nb_inds = dosage.size();
  int nb_snps = 1;
  dosage.clear();
  while( in.read_line(dosage, snp_id, snp_pos, chr, A1, A2) ) {
    nb_snps++;
    if(nb_inds != dosage.size()) {
      Rcerr << "While reading SNP #" << nb_snps << " with id = " << snp_id << "\n";
      Rcerr << "Read " << dosage.size() << " dosages, but there are " << nb_inds << " individuals\n";
      stop("File format error");
    }
    dosage.clear();
  }
  return NumericVector::create(nb_inds, nb_snps);
}

//[[Rcpp::export]]
int nb_inds_dose_file(CharacterVector filename) {
  dosages in( filename );
  std::string snp_id, chr, A1, A2;
  int snp_pos;
  std::vector<double> dosage;
  in.read_line(dosage, snp_id, snp_pos, chr, A1, A2);
  return dosage.size();
}

//[[Rcpp::export]]
CharacterVector samples_dose_file(CharacterVector filename) {
  dosages in( filename );
  return wrap(in.samples);
}

