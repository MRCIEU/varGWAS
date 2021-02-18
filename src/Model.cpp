//
// Created by Matt Lyon on 10/02/2021.
//

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <utility>
#include <vector>
#include <iostream>
#include <boost/math/distributions/students_t.hpp>
#include "Model.h"
#include "PhenotypeFile.h"
#include "Result.h"

/*
 * Class to perform association testing
 * */
namespace jlst {

void Model::run() {
  std::vector<Result> results;
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;

  // Create Eigen matrix of phenotypes wo dosage
  Eigen::MatrixXd X = Eigen::MatrixXd(_phenotype_file.GetNSamples(), _phenotype_file.GetCovariateColumn().size() + 2);
  Eigen::VectorXd y = Eigen::VectorXd(_phenotype_file.GetNSamples());

  // Populate Eigen matrices
  for (unsigned i = 0; i < _phenotype_file.GetNSamples(); i++) {
    X(i, 0) = 1; // intercept
    X(i, 1) = 0; // dosage values are initially set to zero
    for (unsigned j = 0; j < _phenotype_file.GetCovariateColumn().size(); j++) {
      X(i, j + 2) = _phenotype_file.GetCovariateColumn()[j][i]; // covariates
    }
    y(i, 0) = _phenotype_file.GetOutcomeColumn()[i]; // outcome
  }

  // Read variant-by-variant
  while (_bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    LOG_EVERY_N(INFO, 1000) << "Processed " << google::COUNTER << "th variant";

    // only support bi-allelic variants
    if (alleles.size() != 2) {
      LOG(WARNING) << "Skipping multi-allelic variant: " << rsid;
      continue;
    }

    // convert probabilities to dosage values
    _bgen_parser.read_probs(&probs);

    dosages.clear();
    for (auto &prob : probs) {
      // only support bi-allelic variants [0, 1, 2 copies of alt]
      assert(prob.size() == 3);

      // convert genotype probabilities to copies of alt
      if (prob[0] == -1 && prob[1] == -1 && prob[2] == -1) {
        dosages.push_back(-1);
      } else {
        dosages.push_back(prob[1] + (2 * prob[2]));
      }
    }

    // check no missing values between sample list and dosage
    assert(dosages.size() == _phenotype_file.GetNSamples());

    // Create result struct
    Result res;
    res.chromosome = chromosome;
    res.position = position;
    res.rsid = rsid;
    res.effect_allele = alleles[1];
    res.other_allele = alleles[0];

    // TODO multi-thread function

    // store result
    results.push_back(res);

    // pass result-by-reference to allow population in the fit function
    Model::fit(results.back(), dosages, X, y);

    // write to file
    _sf->write(results.back());
  }
}

void Model::fit(Result &result, std::vector<double> dosages, Eigen::MatrixXd X, Eigen::VectorXd y) {
  int n = X.rows();
  int p = X.cols();
  std::vector<unsigned> nulls;

  // set dosage values
  // X is passed without reference to allow for modification on each thread
  assert(dosages.size() == n);
  for (unsigned i = 0; i < dosages.size(); i++) {
    if (dosages[i] == -1) {
      nulls.push_back(i);
    } else {
      X(i, 1) = dosages[i];
    }
  }

  // subset data with missing values
  for (unsigned i = 0; i < nulls.size(); ++i) {
    remove_row_mat(X, i);
    remove_row_vec(y, i);
  }

  // first stage model
  Eigen::BDCSVD<Eigen::MatrixXd> solver(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd fit1 = solver.solve(y);
  Eigen::VectorXd y_hat = X * fit1;
  Eigen::VectorXd y_delta = y - y_hat;
  Eigen::VectorXd y_deltasq = y_delta.array().square();

  // second stage model
  Eigen::MatrixXd fit2 = solver.solve(y_deltasq);

  // epsilon variance
  double e_var = (y_deltasq - X * fit2).squaredNorm() / (n - p);

  // coeff se
  Eigen::MatrixXd ViD = solver.matrixV() * solver.singularValues().asDiagonal().inverse();
  Eigen::MatrixXd varBetasSvd = e_var * ViD * ViD.transpose();
  Eigen::VectorXd se = varBetasSvd.diagonal().array().sqrt();

  // t-stat
  Eigen::VectorXd tstat = fit2.array() / se.array();

  // P
  std::vector<double> pvalues = get_p(tstat, n, p);

  // set results
  result.beta = fit2(1, 0);
  result.se = se(1, 0);
  result.pval = pvalues[1];
  result.n = X.rows();
  result.eaf = X.col(1).mean() * 0.5;
}

std::vector<double> Model::get_p(Eigen::VectorXd &tstat, int n, int p) {
  boost::math::students_t dist(n - (p + 1)); // use student's t-distribution to compute p-value
  std::vector<double> pvalues;
  for (int i = 0; i < tstat.size(); i++) {
    double pval = 2.0 * cdf(complement(dist, tstat[i] > 0 ? tstat[i] : (0 - tstat[i])));
    pvalues.push_back(pval);
  }
  return pvalues;
}

void Model::remove_row_mat(Eigen::MatrixXd &matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();

  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

  matrix.conservativeResize(numRows, numCols);
}

void Model::remove_row_vec(Eigen::VectorXd &vec, unsigned int rowToRemove) {
  unsigned int numRows = vec.rows() - 1;
  unsigned int numCols = vec.cols();

  if (rowToRemove < numRows)
    vec.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        vec.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

  vec.conservativeResize(numRows, numCols);
}

}