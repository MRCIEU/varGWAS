//
// Created by Matt Lyon on 10/02/2021.
//

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <utility>
#include <vector>
#include <iostream>
#include "Model.h"
#include "PhenotypeFile.h"
#include "BgenParser.h"
#include "Result.h"

/*
 * Class to perform association testing
 * */
namespace jlst {
void Model::run(jlst::PhenotypeFile &phenotype_file, genfile::bgen::BgenParser &bgen_parser, int threads,
                std::string &out_file) {
  LOG(INFO) << "Running model";

  // Create Eigen matrix of phenotypes wo dosage
  Eigen::MatrixXd X = Eigen::MatrixXd(phenotype_file.GetNSamples(), phenotype_file.GetCovariateColumn().size() + 2);
  Eigen::VectorXd y = Eigen::VectorXd(phenotype_file.GetNSamples());

  // Populate matrix
  for (unsigned i = 0; i < phenotype_file.GetNSamples(); i++) {
    X(i, 0) = 1; // intercept
    X(i, 1) = 0; // dosage set to zero
    for (unsigned j = 0; j < phenotype_file.GetCovariateColumn().size(); j++) {
      X(i, j + 2) = phenotype_file.GetCovariateColumn()[j][i]; // covariates
    }
    y(i, 0) = phenotype_file.GetOutcomeColumn()[i]; // outcome
  }

  // Read variant-by-variant
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;
  while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    LOG_EVERY_N(INFO, 1000) << "Processed " << google::COUNTER << "th variant";

    // only support bi-allelic variants
    if (alleles.size() != 2) {
      LOG(WARNING) << "Skipping non bi-allelic variant: " << rsid;
      continue;
    }

    // convert probabilities to dosage values
    bgen_parser.read_probs(&probs);

    for (auto &prob : probs) {
      // only support bi-allelic variants [0, 1, 2 copies of alt]
      assert(prob.size() == 3);

      // convert genotype probabilities to copies of alt
      // TODO check for null values
      dosages.push_back(prob[1] + (2 * prob[2]));
    }

    // check no missing values between sample list and dosage
    assert(dosages.size() == phenotype_file.GetNSamples());

    // TODO multi-thread
    LOG(INFO) << "Running with " << threads << " threads";
    Result result = Model::fit(chromosome, position, rsid, alleles, dosages, X, y);
  }
}

Result Model::fit(std::string chromosome,
                  uint32_t position,
                  std::string rsid,
                  std::vector<std::string> alleles,
                  std::vector<double> dosages,
                  Eigen::MatrixXd X,
                  const Eigen::VectorXd &y) {

  // set dosage values
  for (unsigned i = 0; i < dosages.size(); i++) {
    X(i, 1) = dosages[i];
  }

  // linear regression using SVD
  // adapted from: https://genome.sph.umich.edu/w/images/2/2c/Biostat615-lecture14-presentation.pdf
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd betasSvd = svd.solve(y);

  // calculate VDˆ{-1}
  Eigen::MatrixXd ViD = svd.matrixV() * svd.singularValues().asDiagonal().inverse();
  double sigmaSvd = (y - X * betasSvd).squaredNorm() / (X.rows() - X.cols()); // compute \sigmaˆ2
  Eigen::MatrixXd varBetasSvd = sigmaSvd * ViD * ViD.transpose(); // Cov(\hat{beta})

  // Build results struct
  jlst::Result res;
  res.chromosome = std::move(chromosome);
  res.position = position;
  res.rsid = std::move(rsid);
  res.effect_allele = alleles[1];
  res.other_allele = alleles[0];
  res.beta = betasSvd(1, 0);
  res.se = varBetasSvd(1, 0);
  // TODO
  res.pval = 0;

  // second stage model
  // TODO

  return res;
}
}