//
// Created by Matt Lyon on 10/02/2021.
//
#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "PhenotypeFile.h"
#include "BgenParser.h"
#include "Result.h"

#ifndef JLST_CPP_SRC_ASSOC_H_
#define JLST_CPP_SRC_ASSOC_H_

namespace jlst {
class Model {
 public:
  static void run(jlst::PhenotypeFile &phenotype_file, genfile::bgen::BgenParser &bgen_parser);
  static Result fit(std::string chromosome,
                  uint32_t position,
                  std::string rsid,
                  std::vector<std::string> alleles,
                  std::vector<double> dosages,
                  Eigen::MatrixXd X,
                  const Eigen::VectorXd &y);
 private:
};
}

#endif //JLST_CPP_SRC_ASSOC_H_