//
// Created by Matt Lyon on 10/02/2021.
//

#include <omp.h>
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
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;

  spdlog::info("Using {} thread(s)", _threads);
  omp_set_num_threads(_threads);

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

  // open output file
  std::ofstream file(_output_file);
  if (file.is_open()) {
    file << "CHR\tPOS\tRSID\tOA\tEA\tBETA_x\tSE_x\tBETA_xsq\tSE_xsq\tP\tN\tEAF" << std::endl;
    file.flush();
#pragma omp parallel default(none) shared(file, X, y) private(chromosome, position, rsid, alleles, probs, dosages)
    {
#pragma omp master
      // Read variant-by-variant
      while (_bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {

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
          if (prob.size() != 3) {
            throw std::runtime_error("Found " + std::to_string(prob.size()) + " genotypes but we expect three");
          }

          // convert genotype probabilities to copies of alt
          if (prob[0] == -1 && prob[1] == -1 && prob[2] == -1) {
            dosages.push_back(-1);
          } else {
            dosages.push_back(prob[1] + (2 * prob[2]));
          }
        }

        // check no missing values between sample list and dosage
        if (dosages.size() != _phenotype_file.GetNSamples()) {
          throw std::runtime_error("Number of dosage values does not match number of participants");
        }

        // run test and write to file
#pragma omp task
        {
          spdlog::debug("rsid = {}, thread = {}", rsid, omp_get_thread_num());
          Result res = fit(chromosome, position, rsid, alleles[1], alleles[0], dosages, _non_null_idx, X, y);
#pragma omp critical
          {
            // TODO implement file buffer to improve performance
            file << res.chromosome << "\t" << res.position << "\t" << res.rsid << "\t" << res.other_allele << "\t"
                 << res.effect_allele << "\t" << res.beta_x << "\t" << res.se_x << "\t" << res.beta_xsq << "\t"
                 << res.se_xsq << "\t" << res.pval << "\t" << res.n << "\t"
                 << res.eaf << "\n";
            file.flush();
          }
        }
      }
#pragma omp taskwait
    }

    file.close();
  } else {
    throw std::runtime_error("Could not open file: " + _output_file);
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
  res.beta_x = -1;
  res.se_x = -1;
  res.beta_xsq = -1;
  res.se_xsq = -1;
  res.pval = -1;
  res.n = -1;
  res.eaf = -1;

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

  // first stage model
  spdlog::debug("Checking for rank deficiency for first fit");
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_complete);
  if (qr.rank() < X_complete.cols()) {
    return res;
  }
  spdlog::debug("Estimating first stage model");
  Eigen::MatrixXd fs_fit = qr.solve(y_complete);
  Eigen::VectorXd fs_fitted = X_complete * fs_fit;
  Eigen::VectorXd fs_resid = y_complete - fs_fitted;
  Eigen::VectorXd fs_resid2 = fs_resid.array().square();

  // second stage model
  /*Eigen::MatrixXd X_complete2 = Eigen::MatrixXd(X_complete.rows(), 3);
  X_complete2[0] = 1;
  X_complete2[1] = X_complete.col(1);
  X_complete2[2] = X_complete.col(1).square();

  spdlog::debug("Estimating second stage model");
  Eigen::MatrixXd fit2 = qr.solve(fs_resid2);
  Eigen::VectorXd ss_fitted = X_complete * fit2;
  Eigen::VectorXd ss_resid = fs_resid2 - ss_fitted;

  // se
  spdlog::debug("Estimating SE");
  double sig2 = ss_resid.squaredNorm();
  long df = X_complete.rows() - X_complete.cols();
  Eigen::MatrixXd XtX = (X_complete.transpose() * X_complete);
  Eigen::MatrixXd vcov = XtX.inverse();
  Eigen::VectorXd se = (vcov * (sig2 / df)).diagonal().cwiseSqrt();

  // F-test


  // set results
  //res.beta_x = fit2(1, 0);
  //res.beta_xsq = fit2(2, 0);
  res.se_x = se(1, 0);
  res.se_xsq = se(2, 0);
  //res.pval = pvalues[1];
  res.n = X_complete.rows();
  res.eaf = X_complete.col(1).mean() * 0.5;*/

  return res;
}

}