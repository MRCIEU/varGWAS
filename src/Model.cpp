//
// Created by Matt Lyon on 10/02/2021.
//
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include "Model.h"
#include "PhenotypeFile.h"
#include "Result.h"
#include "spdlog/spdlog.h"
extern "C" {
#include <plinkio/plinkio.h>
}

/*
 * Class to perform association testing
 * */
namespace jlst {

void Model::parse_bgen(genfile::bgen::BgenParser &bgen_parser) {
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<std::vector<double>> probs;
  std::vector<double> dosages;

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
    file << "chr\tpos\trsid\toa\tea\tbeta\tse\tt\tp\tphi_x\tse_x\tphi_xsq\tse_xsq\tphi_f\tphi_p\tn\teaf" << std::endl;
    file.flush();
#pragma omp parallel default(none) shared(file, bgen_parser, X, y) private(chromosome, position, rsid, alleles, probs, dosages)
    {
#pragma omp master
      // Read variant-by-variant
      while (bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {

        // only support bi-allelic variants
        if (alleles.size() != 2) {
          spdlog::warn("Skipping multi-allelic variant: {}", rsid);
          continue;
        }

        // convert probabilities to dosage values
        bgen_parser.read_probs(&probs);

        spdlog::debug("Converting probabilities to dosage values");
        dosages.clear();
        for (auto &prob : probs) {
          // only support bi-allelic variants [0, 1, 2 copies of alt]
          if (prob.size() != 3) {
            throw std::runtime_error("Found " + std::to_string(prob.size()) + " genotypes but we expect three");
          }

          // check dosage values are within expected range
          if (prob[0] < 0 || prob[0] > 1 || prob[1] < 0 || prob[1] > 1 || prob[2] < 0 || prob[2] > 1) {
            throw std::runtime_error(
                "Dosage value outside expected range: " + std::to_string(prob[0]) + " " + std::to_string(prob[1]) + " "
                    + std::to_string(prob[2])
            );
          }

          // convert genotype probabilities to copies of alt
          // TODO how are null values recorded in BGEN?
          if ((prob[0] == -1 && prob[1] == -1 && prob[2] == -1) || (prob[0] == 0 && prob[1] == 0 && prob[2] == 0)) {
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
                 << res.effect_allele << "\t" << res.beta << "\t" << res.se << "\t" << res.t << "\t"
                 << res.pval << "\t" << res.phi_x << "\t" << res.se_x << "\t" << res.phi_xsq << "\t"
                 << res.se_xsq << "\t" << res.phi_f << "\t" << res.phi_pval << "\t" << res.n << "\t" << res.eaf << "\n";
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

void Model::parse_plink(std::string &file_path) {
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<double> dosages;
  struct pio_file_t plink_file;
  snp_t *snp_buffer;
  int sample_id;
  int locus_id;

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
    file << "chr\tpos\trsid\toa\tea\tbeta\tse\tt\tp\tphi_x\tse_x\tphi_xsq\tse_xsq\tphi_f\tphi_p\tn\teaf" << std::endl;
    file.flush();
#pragma omp parallel default(none) shared(file, file_path, plink_file, X, y) private(snp_buffer, sample_id, locus_id, chromosome, position, rsid, alleles, dosages)
    {
#pragma omp master
      // Read variant-by-variant
      if (pio_open(&plink_file, file_path.c_str()) != PIO_OK) {
        throw std::runtime_error("Could not open file: " + file_path);
      }

      if (!pio_one_locus_per_row(&plink_file)) {
        throw std::runtime_error("This script requires that snps are rows and samples columns");
      }

      locus_id = 0;
      snp_buffer = (snp_t *) malloc(pio_row_size(&plink_file));

      // loop over variants
      while (pio_next_row(&plink_file, snp_buffer) == PIO_OK) {
        struct pio_locus_t *locus = pio_get_locus(&plink_file, locus_id);
        dosages.clear();

        // loop over samples
        for (sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++) {
          struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
          int dosage = snp_buffer[sample_id];
          if (dosage = 3) {
            dosages.push_back(-1);
          } else {
            dosages.push_back(dosage);
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
                 << res.effect_allele << "\t" << res.beta << "\t" << res.se << "\t" << res.t << "\t"
                 << res.pval << "\t" << res.phi_x << "\t" << res.se_x << "\t" << res.phi_xsq << "\t"
                 << res.se_xsq << "\t" << res.phi_f << "\t" << res.phi_pval << "\t" << res.n << "\t" << res.eaf
                 << "\n";
            file.flush();
          }
        }

        locus_id++;
      }

      free(snp_buffer);
      pio_close(&plink_file);
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
  res.eaf = -1;
  res.n = -1;
  res.beta = -1;
  res.se = -1;
  res.t = -1;
  res.pval = -1;
  res.phi_x = -1;
  res.se_x = -1;
  res.phi_xsq = -1;
  res.se_xsq = -1;
  res.phi_f = -1;
  res.phi_pval = -1;

  // set dosage values
  // X is passed without reference to allow for modification on each thread
  spdlog::debug("Checking for null dosage values");
  if (dosages.size() != X.rows()) {
    throw std::runtime_error("Number of dosage values does not match number of participants");
  }
  int j = 0;
  for (unsigned i = 0; i < dosages.size(); i++) {
    if (dosages[i] == -1) {
      j++;
      non_null_idx.erase(i);
    } else {
      X(i, 1) = dosages[i];
    }
  }
  spdlog::debug("Found {} null values", j);

  // subset data with non-null values
  spdlog::debug("Selecting non-null values for model");
  std::vector<unsigned> non_null_idx_vec(non_null_idx.begin(), non_null_idx.end());
  Eigen::MatrixXd X_complete = X(non_null_idx_vec, Eigen::all).eval();
  Eigen::VectorXd y_complete = y(non_null_idx_vec, Eigen::all).eval();
  int n = X_complete.rows();

  // first stage model
  spdlog::debug("Checking for rank deficiency for first fit");
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr1(X_complete);
  if (qr1.rank() < X_complete.cols()) {
    return res;
  }
  spdlog::debug("Estimating first stage model");
  Eigen::MatrixXd fs_fit = qr1.solve(y_complete);
  Eigen::VectorXd fs_fitted = X_complete * fs_fit;
  Eigen::VectorXd fs_resid = y_complete - fs_fitted;
  Eigen::VectorXd fs_resid2 = fs_resid.array().square();

  // se
  spdlog::debug("Estimating SE");
  double fs_sig2 = fs_resid.squaredNorm();
  long fs_df = X_complete.rows() - X_complete.cols();
  Eigen::MatrixXd fs_XtX = (X_complete.transpose() * X_complete);
  Eigen::MatrixXd fs_vcov = fs_XtX.inverse();
  Eigen::VectorXd fs_se = (fs_vcov * (fs_sig2 / fs_df)).diagonal().cwiseSqrt();

  // pval
  spdlog::debug("Estimating P");
  double t_stat = fs_fit(1, 0) / fs_se(1, 0);
  boost::math::students_t t_dist(fs_df);
  double pval = 2.0 * boost::math::cdf(boost::math::complement(t_dist, fabs(t_stat)));

  // second stage model
  spdlog::debug("Estimating second stage model");
  Eigen::MatrixXd X_complete2 = Eigen::MatrixXd(n, 3);
  X_complete2.col(0).setOnes();
  X_complete2.col(1) = X_complete.col(1).array();
  X_complete2.col(2) = X_complete.col(1).array().square();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr2(X_complete2);
  if (qr2.rank() < X_complete2.cols()) {
    return res;
  }
  Eigen::MatrixXd ss_fit = qr2.solve(fs_resid2);
  Eigen::VectorXd ss_fitted = X_complete2 * ss_fit;
  Eigen::VectorXd ss_resid = fs_resid2 - ss_fitted;
  Eigen::VectorXd ss_resid2 = ss_resid.array().square();

  // se
  spdlog::debug("Estimating SE");
  double ss_sig2 = ss_resid.squaredNorm();
  long ss_df = X_complete2.rows() - X_complete2.cols();
  Eigen::MatrixXd ss_XtX = (X_complete2.transpose() * X_complete2);
  Eigen::MatrixXd ss_vcov = ss_XtX.inverse();
  Eigen::VectorXd ss_se = (ss_vcov * (ss_sig2 / ss_df)).diagonal().cwiseSqrt();

  // intercept only model
  Eigen::MatrixXd X_complete3 = Eigen::MatrixXd(n, 1);
  X_complete3.col(0).setOnes();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr3(X_complete3);
  if (qr3.rank() < X_complete3.cols()) {
    return res;
  }
  Eigen::MatrixXd null_fit = qr3.solve(fs_resid2);
  Eigen::VectorXd null_fitted = X_complete3 * null_fit;
  Eigen::VectorXd null_resid = fs_resid2 - null_fitted;
  Eigen::VectorXd null_resid2 = null_resid.array().square();

  // F-test
  // adapted from http://people.reed.edu/~jones/Courses/P24.pdf
  int df_f = n - 3;
  int df_r = n - 1;
  int df_n = df_r - df_f;
  double rss_f = ss_resid.squaredNorm();
  double rss_r = null_resid.squaredNorm();
  double f = ((rss_r - rss_f) / (df_r - df_f)) / (rss_f / df_f);
  if (f < 0) {
    return res;
  }
  boost::math::fisher_f f_dist(df_n, df_f);
  double phi_pval = boost::math::cdf(boost::math::complement(f_dist, f));

  // set results
  res.beta = fs_fit(1, 0);
  res.se = fs_se(1, 0);
  res.t = t_stat;
  res.pval = pval;
  res.phi_x = ss_fit(1, 0);
  res.phi_xsq = ss_fit(2, 0);
  res.se_x = ss_se(1, 0);
  res.se_xsq = ss_se(2, 0);
  res.phi_f = f;
  res.phi_pval = phi_pval;
  res.n = n;
  res.eaf = X_complete.col(1).mean() * 0.5;

  return res;
}

}