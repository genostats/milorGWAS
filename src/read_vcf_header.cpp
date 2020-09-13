#include "read_vcf_header.h"

void read_vcf_samples(std::string line, std::vector<std::string> & samples) {
  std::istringstream li(line);
  std::string G;
  for(int i = 0; i < 9; i++) {
    if(!(li >> G))
      stop("VCF file format error");
  }

  while(li >> G) {
    samples.push_back( G );
  }
}


void read_vcf_header(igzstream & in, std::vector<std::string> & samples, std::vector<std::string> & format_ids, std::vector<std::string> & info_ids) {
  // skip header and read samples
  std::string line, id;
  if(!std::getline(in, line))
    stop("File is empty");

  while(std::getline(in, line)) {
    if(line.substr(0,1) != "#") stop("Bad VCF format");
    if(line.substr(0,2) == "##") {
      if(line.substr(0,11) == "##INFO=<ID=") {
        std::istringstream li(line);
        std::getline(li, id, ',');
        info_ids.push_back(id.substr(11));
      }
      if(line.substr(0,13) == "##FORMAT=<ID=") {
        std::istringstream li(line);
        std::getline(li, id, ',');
        format_ids.push_back(id.substr(13));
      }
    } else {
      read_vcf_samples(line, samples);
      break; // fin
    }
  }
}


//[[Rcpp::export]]
List read_vcf_head(std::string filename) {
  std::vector<std::string> SAMPLES, FORMAT_IDS, INFO_IDS;

  // open file
  igzstream in( filename.c_str() );
  if(!in.good()) stop("Can't open file "+filename);

  // read header
  read_vcf_header(in, SAMPLES, FORMAT_IDS, INFO_IDS);
  
  List L;
  L["samples"] = wrap(SAMPLES);
  L["info_ids"] = wrap(INFO_IDS);
  L["format_ids"] = wrap(FORMAT_IDS);

  return L;
}

RcppExport SEXP gg_read_vcf_head(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(read_vcf_head(filename));
    return rcpp_result_gen;
END_RCPP
}

