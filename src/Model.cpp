//
// Created by Matt Lyon on 10/02/2021.
//

#include <glog/logging.h>
#include "Model.h"
#include "PhenotypeFile.h"
#include "BgenParser.h"
#include "Result.h"
#include <Eigen/Core>
#include <Eigen/SVD>

/*
 * Class to perform association testing
 * */
namespace jlst {
void Model::run(jlst::PhenotypeFile &phenotype_file, genfile::bgen::BgenParser &bgen_parser) {

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

  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;

  // Read variant-by-variant
  while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {
    LOG_EVERY_N(INFO, 1000) << "Read the " << google::COUNTER << "th variant";

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

    // fit model
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

  // Build results struct
  jlst::Result res;
  res.chromosome = chromosome;
  res.position = position;
  res.rsid = rsid;
  res.effect_allele = alleles[1];
  res.other_allele = alleles[0];
  // Perform model
  // TODO add dosage values to Eigen matrix
  return res;
}
}