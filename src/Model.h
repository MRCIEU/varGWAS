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
#include <ThreadPool.h>
#include "PhenotypeFile.h"
#include "BgenParser.h"
#include "Result.h"
#include "SynchronizedFile.h"

#ifndef JLST_CPP_SRC_MODEL_H_
#define JLST_CPP_SRC_MODEL_H_

namespace jlst {
class Model {
 public:
  Model(jlst::PhenotypeFile &phenotype_file,
        genfile::bgen::BgenParser &bgen_parser,
        std::set<unsigned> &non_null,
        std::shared_ptr<SynchronizedFile> sf,
        int threads)
      : _phenotype_file(phenotype_file),
        _bgen_parser(bgen_parser),
        _non_null_idx(non_null_idx),
        _sf(std::move(sf)),
        _threads(threads) {}
  void run();
  static Result fit(std::string &chromosome,
                    uint32_t position,
                    std::string &rsid,
                    std::string &effect_allele,
                    std::string &other_allele,
                    std::vector<double> dosages,
                    std::set<unsigned> non_null_idx,
                    Eigen::MatrixXd X,
                    Eigen::VectorXd y);
  static std::vector<double> get_p(Eigen::VectorXd &tstat, int n, int p);
 private:
  jlst::PhenotypeFile &_phenotype_file;
  genfile::bgen::BgenParser &_bgen_parser;
  std::set<unsigned> &_non_null_idx;
  int _threads;
  std::shared_ptr<SynchronizedFile> _sf;
};
}

#endif //JLST_CPP_SRC_MODEL_H_
