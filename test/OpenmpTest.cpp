#include <iostream>
#include <utility>
#include "gtest/gtest.h"
#include <omp.h>
#include <Result.h>
#include "BgenParser.h"

vargwas::Result func(std::string chromosome,
                  uint32_t position,
                  std::string rsid,
                  std::vector<std::string> alleles,
                  std::vector<std::vector<double>> probs) {

  vargwas::Result res;

  res.chromosome = std::move(chromosome);
  res.position = position;
  res.rsid = std::move(rsid);
  res.effect_allele = alleles[0];
  res.other_allele = alleles[1];
  res.beta = probs[0][0];
  res.se = probs[0][0];
  res.pval = probs[0][0];
  res.n = probs[0][0];
  res.eaf = probs[0][0];

  return (res);
}

TEST(OpenMPTest, process_variants) {
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  genfile::bgen::BgenParser bgen_parser("/Users/ml18692/projects/vargwas/lib/bgen/example/example.v11.bgen");
  std::stringstream buf;
  unsigned i = 0;

  // adapted from https://stackoverflow.com/questions/10678605/idiomatic-way-to-parallelize-function-across-file-lines-in-c
#pragma omp parallel
  {
#pragma omp master
    while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
      bgen_parser.read_probs(&probs);
#pragma omp task
      {
        vargwas::Result res = func(chromosome, position, rsid, alleles, probs);
#pragma omp critical
        {
          buf << res.chromosome << "\t" << res.position << "\t" << res.rsid << "\t" << res.other_allele << "\t"
              << res.effect_allele << "\t" << res.beta << "\t" << res.se << "\t" << res.pval << "\t" << res.n
              << "\t"
              << res.eaf << "\n";
          buf.flush();
        }
      }
      i++;
#pragma omp taskwait
    }
  }

  std::string line;
  int n_lines = 0;
  while (std::getline(buf, line)) {
    n_lines++;
  }

  ASSERT_EQ(i, n_lines);
}
