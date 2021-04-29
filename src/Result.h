//
// Created by Matt Lyon on 10/02/2021.
//

#ifndef JLST_CPP_SRC_RESULT_H_
#define JLST_CPP_SRC_RESULT_H_

namespace jlst {
struct Result {
  std::string chromosome{};
  uint32_t position{};
  std::string rsid{};
  std::string effect_allele{};
  std::string other_allele{};
  double beta{};
  double se{};
  double t{};
  double pval{};
  double phi_x{};
  double se_x{};
  double phi_xsq{};
  double se_xsq{};
  double phi_f{};
  double phi_pval{};
  int n{};
  double eaf{};
};
}

#endif //JLST_CPP_SRC_RESULT_H_
