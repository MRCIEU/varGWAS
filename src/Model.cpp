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
#include "BgenParser.h"
#include "Result.h"

/*
 * Class to perform association testing
 * */
namespace jlst {
std::vector<Result> Model::run(jlst::PhenotypeFile &phenotype_file,
                               genfile::bgen::BgenParser &bgen_parser,
                               int threads) {

  std::vector<Result> results;
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;

  // Create Eigen matrix of phenotypes wo dosage
  Eigen::MatrixXd X = Eigen::MatrixXd(phenotype_file.GetNSamples(), phenotype_file.GetCovariateColumn().size() + 2);
  Eigen::VectorXd y = Eigen::VectorXd(phenotype_file.GetNSamples());

  // Populate Eigen matrices
  for (unsigned i = 0; i < phenotype_file.GetNSamples(); i++) {
    X(i, 0) = 1; // intercept
    X(i, 1) = 0; // dosage values are initially set to zero
    for (unsigned j = 0; j < phenotype_file.GetCovariateColumn().size(); j++) {
      X(i, j + 2) = phenotype_file.GetCovariateColumn()[j][i]; // covariates
    }
    y(i, 0) = phenotype_file.GetOutcomeColumn()[i]; // outcome
  }

  // Read variant-by-variant
  while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    LOG_EVERY_N(INFO, 1000) << "Processed " << google::COUNTER << "th variant";

    std::cout << "chr " << chromosome << std::endl;

    // only support bi-allelic variants
    if (alleles.size() != 2) {
      LOG(WARNING) << "Skipping multi-allelic variant: " << rsid;
      continue;
    }

    // convert probabilities to dosage values
    bgen_parser.read_probs(&probs);

    dosages.clear();
    for (auto &prob : probs) {
      // only support bi-allelic variants [0, 1, 2 copies of alt]
      assert(prob.size() == 3);

      // convert genotype probabilities to copies of alt
      // TODO check for null values
      dosages.push_back(prob[1] + (2 * prob[2]));
    }

    // check no missing values between sample list and dosage
    assert(dosages.size() == phenotype_file.GetNSamples());

    // Create result struct
    Result res;
    res.chromosome = chromosome;
    res.position = position;
    res.rsid = rsid;
    res.effect_allele = alleles[1]; // TODO check allele order
    res.other_allele = alleles[0];

    // store result
    results.push_back(res);

    // pass result-by-reference to allow population in the fit function
    // TODO multi-thread function
    Model::fit(results.back(), dosages, X, y);
  }

  // TODO write to file to avoid memory issues
  return results;
}

void Model::fit(Result &result, std::vector<double> dosages, Eigen::MatrixXd X, const Eigen::VectorXd &y) {
  int n = X.rows();
  int p = X.cols();

  // set dosage values
  // X is passed without reference to allow for modification on each thread
  for (unsigned i = 0; i < dosages.size(); i++) {
    X(i, 1) = dosages[i];
  }

  // first stage model
  Eigen::BDCSVD<Eigen::MatrixXd> solver(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd fit1 = solver.solve(y);
  Eigen::VectorXd y_hat = X * fit1;
  Eigen::VectorXd y_delta = y - y_hat;
  Eigen::VectorXd y_deltasq = y_delta.cwiseProduct(y_delta);

  // second stage model
  Eigen::MatrixXd fit2 = solver.solve(y_deltasq);

  // epsilon variance
  double e_var = (y - X * fit2).squaredNorm() / (n - p);

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

}