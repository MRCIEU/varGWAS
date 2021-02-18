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
#include "SynchronizedFile.h"

#ifndef JLST_CPP_SRC_ASSOC_H_
#define JLST_CPP_SRC_ASSOC_H_

namespace jlst {
class Model {
 public:
  static void run(jlst::PhenotypeFile &phenotype_file,
           genfile::bgen::BgenParser &bgen_parser,
           jlst::SynchronizedFile &output_file,
           int threads);
  static void fit(Result &result, std::vector<double> dosages, Eigen::MatrixXd X, Eigen::VectorXd y);
  static std::vector<double> get_p(Eigen::VectorXd &tstat, int n, int p);
  static void remove_row_mat(Eigen::MatrixXd &matrix, unsigned int rowToRemove);
  static void remove_row_vec(Eigen::VectorXd &vec, unsigned int rowToRemove);
 private:
};
}

#endif //JLST_CPP_SRC_ASSOC_H_
