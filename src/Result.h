//
// Created by Matt Lyon on 10/02/2021.
//

#ifndef VARGWAS_CPP_SRC_RESULT_H_
#define VARGWAS_CPP_SRC_RESULT_H_

namespace vargwas {
struct Result {
  std::string chromosome{};
  uint32_t position{};
  std::string rsid{};
  std::string effect_allele{};
  std::string other_allele{};
  double beta{}; // mean effect of SNP on outcome
  double se{}; // SE for mean effect of SNP on outcome
  double t{}; // t-stat for mean effect of SNP on outcome
  double pval{}; // p-value for mean effect of SNP on outcome
  double theta{}; // median effect of SNP on outcome
  double phi_x1{}; // var(Y|G==1) - var(Y|G==0)
  double se_x1{}; // SE of phi_x1
  double phi_x2{}; // var(Y|G==2) - var(Y|G==0)
  double se_x2{}; // SE of phi_x2
  double phi_f{}; // F-stat for difference in variance
  double phi_pval{}; // p-value for difference in variance
  int n{}; // sample size
  double eaf{}; // effect allele freq
};
}

#endif //VARGWAS_CPP_SRC_RESULT_H_
