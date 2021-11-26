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
#include "Result.h"

#ifndef VARGWAS_CPP_SRC_MODEL_H_
#define VARGWAS_CPP_SRC_MODEL_H_

namespace vargwas {
class Model {
 public:
  Model(vargwas::PhenotypeFile &phenotype_file,
        genfile::bgen::BgenParser &bgen_parser,
        std::set<unsigned> &non_null_idx,
        std::string &output_file,
        int threads,
        double maf_threshold,
        bool flip)
      : _phenotype_file(phenotype_file),
        _bgen_parser(bgen_parser),
        _non_null_idx(non_null_idx),
        _output_file(output_file),
        _threads(threads),
        _maf_threshold(maf_threshold),
        _flip(flip)
        {}
  void run();
  static Result fit(std::string &chromosome,
                    uint32_t position,
                    std::string &rsid,
                    std::string &effect_allele,
                    std::string &other_allele,
                    std::vector<double> dosages,
                    std::set<unsigned> non_null_idx,
                    Eigen::MatrixXd X1,
                    Eigen::MatrixXd X2,
                    Eigen::VectorXd y,
                    double maf_threshold);
 private:
  vargwas::PhenotypeFile &_phenotype_file;
  genfile::bgen::BgenParser &_bgen_parser;
  std::set<unsigned> &_non_null_idx;
  int _threads;
  double _maf_threshold;
  bool _flip;
  std::string &_output_file;
  static void first_stage_ols(
      const Eigen::MatrixXd &X_complete1,
      const Eigen::VectorXd &y_complete,
      Eigen::VectorXd &fs_se, Eigen::VectorXd &fs_fit,
      double &pval, double &t_stat
  );
  static void first_stage_lad(
      const Eigen::MatrixXd &X_complete1,
      const Eigen::VectorXd &y_complete,
      const Eigen::VectorXd &fs_fit,
      Eigen::VectorXd &fs_resid2,
      double &theta
  );
  static void second_stage(
      const Eigen::MatrixXd &X_complete2,
      const Eigen::VectorXd &fs_resid2,
      Eigen::VectorXd &ss_resid,
      Eigen::VectorXd &ss_fit
  );
  static Eigen::MatrixXd white_vcov(const Eigen::MatrixXd &X_complete2, const Eigen::VectorXd &ss_resid);
  static void delta_method(const Eigen::VectorXd &ss_fit,
                           const Eigen::MatrixXd &hc_vcov,
                           double &s1_dummy,
                           double &s2_dummy);
  static void null_fit(const Eigen::MatrixXd &X_complete2,
                       const Eigen::VectorXd &fs_resid2,
                       double &phi_pval,
                       double &f,
                       int n,
                       const Eigen::VectorXd &ss_resid);
};
}

#endif //VARGWAS_CPP_SRC_MODEL_H_
