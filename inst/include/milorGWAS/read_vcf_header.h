#include <iostream>
#include <string>
#include <sstream>
#include "gaston/matrix4.h"
#include <fstream>
#include "gaston/gzstream.h"

#ifndef READ_VCF_HEADER
#define READ_VCF_HEADER


void read_vcf_samples(std::string line, std::vector<std::string> & samples);

void read_vcf_header(igzstream & in, std::vector<std::string> & samples, std::vector<std::string> & format_ids, std::vector<std::string> & info_ids);

#endif
