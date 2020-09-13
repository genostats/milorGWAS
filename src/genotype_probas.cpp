#include <iostream>
#include <string>
#include "stringstream_lite.h"
#include "genotype_probas.h"

double genotype_probas_to_dosage(std::string & S) {
  stringstream_lite li(S, ',');
  double p0(0), p1(0), p2(0);
  if(li >> p0 >> p1 >> p2) {
    return p1 + 2.0*p2;
  } else {
    return -1.0;
  }
}

