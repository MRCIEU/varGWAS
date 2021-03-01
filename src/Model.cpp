//
// Created by Matt Lyon on 10/02/2021.
//

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <boost/math/distributions/students_t.hpp>
#include "Model.h"
#include "PhenotypeFile.h"
#include "Result.h"
#include "spdlog/spdlog.h"


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
  unsigned n = 0;

  spdlog::info("Starting {} thread(s)", _threads);
  ThreadPool pool(_threads);

  // Create Eigen matrix of phenotypes wo dosage
  // p+2 for dosage and intercept
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

  spdlog::info("Estimating model with {} samples and {} parameters", _non_null_idx.size(), X.cols());

  // Read variant-by-variant
  while (_bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    n++;
    spdlog::debug("Testing {}th variant: {}", n, rsid);

    // only support bi-allelic variants
    if (alleles.size() != 2) {
      spdlog::warn("Skipping multi-allelic variant: {}", rsid);
      continue;
    }

    // convert probabilities to dosage values
    _bgen_parser.read_probs(&probs);

    spdlog::debug("Converting probabilities to dosage values");
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

    // enqueue and store future
    spdlog::debug("Submitting job to queue");
    auto result = pool.enqueue(fit, chromosome, position, rsid, alleles[1], alleles[0], dosages, _non_null_idx, X, y);

    // write to file
    _sf->write(result.get());
  }
}

Result Model::fit(std::string &chromosome,
                  uint32_t position,
                  std::string &rsid,
                  std::string &effect_allele,
                  std::string &other_allele,
                  std::vector<double> dosages,
                  std::set<unsigned> non_null_idx,
                  Eigen::MatrixXd X,
                  Eigen::VectorXd y) {
  Result res;
  res.chromosome = chromosome;
  res.position = position;
  res.rsid = rsid;
  res.effect_allele = effect_allele;
  res.other_allele = other_allele;

  // set dosage values
  // X is passed without reference to allow for modification on each thread
  spdlog::debug("Checking for null dosage values");
  assert(dosages.size() == X.rows());
  for (unsigned i = 0; i < dosages.size(); i++) {
    if (dosages[i] == -1) {
      non_null_idx.erase(i);
    } else {
      X(i, 1) = dosages[i];
    }
  }

  // subset data with non-null values
  spdlog::debug("Selecting non-null values for model");
  std::vector<unsigned> non_null_idx_vec(non_null_idx.begin(), non_null_idx.end());
  Eigen::MatrixXd X_complete = X(non_null_idx_vec, Eigen::all).eval();
  Eigen::VectorXd y_complete = y(non_null_idx_vec, Eigen::all).eval();

  // model
  spdlog::debug("Checking for rank deficiency");
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_complete);
  if (qr.rank() < X_complete.cols()) {
    throw std::runtime_error("rank-deficient matrix");
  }

  // first stage model
  spdlog::debug("Estimating first stage model");
  Eigen::MatrixXd fs_fit = qr.solve(y_complete);
  Eigen::VectorXd fs_fitted = X_complete * fs_fit;
  Eigen::VectorXd fs_resid = y_complete - fs_fitted;
  Eigen::VectorXd fs_resid2 = fs_resid.array().square();

  // second stage model
  spdlog::debug("Estimating second stage model");
  Eigen::MatrixXd fit2 = qr.solve(fs_resid2);
  Eigen::VectorXd ss_fitted = X_complete * fit2;
  Eigen::VectorXd ss_resid = fs_resid2 - ss_fitted;

  // se
  spdlog::debug("Estimating SE and P value");
  // TODO check df is correct as using multiple models - do we include the second-stage intercept and slope?
  double sig2 = ss_resid.squaredNorm();
  long df = X_complete.rows() - X_complete.cols();
  Eigen::MatrixXd XtX = (X_complete.transpose() * X_complete);
  Eigen::MatrixXd vcov = XtX.inverse();
  Eigen::VectorXd se = (vcov * (sig2 / df)).diagonal().cwiseSqrt();

  // t-stat
  Eigen::VectorXd tstat = fit2.array() / se.array();

  // P
  std::vector<double> pvalues = get_p(tstat, X_complete.rows(), X_complete.cols());

  // set results
  res.beta = fit2(1, 0);
  res.se = se(1, 0);
  res.pval = pvalues[1];
  res.n = X_complete.rows();
  res.eaf = X_complete.col(1).mean() * 0.5;

  return res;
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