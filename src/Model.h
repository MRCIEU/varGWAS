//
// Created by Matt Lyon on 10/02/2021.
//
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <set>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include "PhenotypeFile.h"
#include "BgenParser.h"
#include "PlinkParser.h"
#include "Result.h"

#ifndef JLST_CPP_SRC_MODEL_H_
#define JLST_CPP_SRC_MODEL_H_

namespace jlst {
class Model {
 public:
  Model(jlst::PhenotypeFile &phenotype_file,
        std::set<unsigned> &non_null_idx,
        std::string &output_file,
        int threads)
      : _phenotype_file(phenotype_file),
        _non_null_idx(non_null_idx),
        _output_file(output_file),
        _threads(threads) {
  }
  void parse_bgen(genfile::bgen::BgenParser &bgen_parser);
  void parse_plink(jlst::PlinkParser &plink_parser);
  static Result fit(std::string &chromosome,
                    uint32_t position,
                    std::string &rsid,
                    std::string &effect_allele,
                    std::string &other_allele,
                    std::vector<double> dosages,
                    std::set<unsigned> non_null_idx,
                    Eigen::MatrixXd X,
                    Eigen::VectorXd y);
 private:
  jlst::PhenotypeFile &_phenotype_file;
  std::set<unsigned> &_non_null_idx;
  int _threads;
  std::string &_output_file;
};
}

#endif //JLST_CPP_SRC_MODEL_H_