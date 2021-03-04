#include <iostream>
#include "gtest/gtest.h"
#include <omp.h>
#include "BgenParser.h"

TEST(BgenTest, read_variants_openmp) {
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  genfile::bgen::BgenParser bgen_parser("/Users/ml18692/projects/jlst_cpp/lib/bgen/example/example.v11.bgen");

  // parse BGEN file
  while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    bgen_parser.read_probs(&probs);
    std::cout << chromosome << " " << position << " " << rsid << " " << alleles[0] << " " << alleles[1] << std::endl;
  }

}