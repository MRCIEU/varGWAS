//
// Created by Matt Lyon on 10/02/2021.
//

#include <omp.h>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <math.h>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include "Model.h"
#include "PhenotypeFile.h"
#include "Result.h"
#include "QRMM.h"
#include "spdlog/spdlog.h"

/*
 * Class to perform association testing
 * */
namespace vargwas {

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
  // create place for dosage and intercept
  Eigen::MatrixXd
      X1 = Eigen::MatrixXd(_phenotype_file.GetNSamples(), 1 + 1 + _phenotype_file.GetCovariateColumn().size());
  Eigen::MatrixXd
      X2 = Eigen::MatrixXd(_phenotype_file.GetNSamples(), 1 + 2 + _phenotype_file.GetCovariateColumn().size());
  Eigen::VectorXd y = Eigen::VectorXd(_phenotype_file.GetNSamples());

  // Populate Eigen matrices
  for (unsigned i = 0; i < _phenotype_file.GetNSamples(); i++) {
    X1(i, 0) = 1; // intercept
    X2(i, 0) = 1; // intercept
    X1(i, 1) = 0; // dosage values are initially set to zero
    X2(i, 1) = 0; // dosage values are initially set to zero
    X2(i, 2) = 0; // dosage values are initially set to zero
    for (unsigned j = 0; j < _phenotype_file.GetCovariateColumn().size(); j++) {
      X1(i, j + 2) = _phenotype_file.GetCovariateColumn()[j][i]; // covariates
      X2(i, j + 3) = _phenotype_file.GetCovariateColumn()[j][i]; // covariates
    }
    y(i, 0) = _phenotype_file.GetOutcomeColumn()[i]; // outcome
  }

  spdlog::info("Estimating model with {} samples and {} parameters", _non_null_idx.size(), X1.cols());

  // open output file
  std::ofstream file(_output_file);
  if (file.is_open()) {
    file << "chr\tpos\trsid\toa\tea\tn\teaf\tbeta\tse\tt\tp\ttheta\tphi_x1\tse_x1\tphi_x2\tse_x2\tphi_f\tphi_p"
         << std::endl;
    file.flush();
#pragma omp parallel default(none) shared(file, X1, X2, y) private(chromosome, position, rsid, alleles, probs, dosages)
    {
#pragma omp master
      // Read variant-by-variant
      while (_bgen_parser.read_variant(&chromosome, &position, &rsid, &alleles)) {

        // only support bi-allelic variants
        if (alleles.size() != 2) {
          spdlog::warn("Skipping multi-allelic variant: {}, thread = {}", rsid, omp_get_thread_num());
          continue;
        }

        if (flip) {
          std::string a0 = alleles[1];
          std::string a1 = alleles[0];
        } else {
          std::string a0 = alleles[0];
          std::string a1 = alleles[1];
        }

        // convert probabilities to dosage values
        _bgen_parser.read_probs(&probs);

        spdlog::debug("Converting probabilities to dosage values, thread = {}", omp_get_thread_num());
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

          // convert genotype probabilities to copies of A1
          // TODO how are null values recorded in BGEN?
          if ((prob[0] == -1 && prob[1] == -1 && prob[2] == -1) || (prob[0] == 0 && prob[1] == 0 && prob[2] == 0)) {
            dosages.push_back(-1);
          } else {
            if (flip) {
              dosages.push_back(prob[1] + (2 * prob[2]));
            } else {
              dosages.push_back(prob[1] + (2 * prob[0]));
            }
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
          Result
              res = fit(chromosome,
                        position,
                        rsid,
                        a0,
                        a1,
                        dosages,
                        _non_null_idx,
                        X1,
                        X2,
                        y,
                        _maf_threshold);
          if (!isnan(res.phi_pval)) {
#pragma omp critical
            {
              // TODO implement file buffer to improve performance
              file << res.chromosome << "\t" << res.position << "\t" << res.rsid << "\t" << res.other_allele << "\t"
                   << res.effect_allele << "\t" << res.n << "\t" << res.eaf << "\t" << res.beta << "\t" << res.se
                   << "\t" << res.t << "\t" << res.pval << "\t" << res.theta << "\t" << res.phi_x1 << "\t" << res.se_x1
                   << "\t" << res.phi_x2 << "\t" << res.se_x2 << "\t" << res.phi_f << "\t" << res.phi_pval << "\n";
              file.flush();
            }
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

void Model::first_stage_ols(
    const Eigen::MatrixXd &X_complete1,
    const Eigen::VectorXd &y_complete,
    Eigen::VectorXd &fs_se, Eigen::VectorXd &fs_fit,
    double &pval, double &t_stat
) {

  spdlog::debug("Checking for rank deficiency for first fit, thread = {}", omp_get_thread_num());
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr1(X_complete1);
  if (qr1.rank() < X_complete1.cols()) {
    throw std::runtime_error("Rank deficiency");
  }
  spdlog::debug("Estimating first stage model using OLS regression, thread = {}", omp_get_thread_num());
  fs_fit = qr1.solve(y_complete);
  Eigen::VectorXd fs_fitted = X_complete1 * fs_fit;
  Eigen::VectorXd fs_resid = y_complete - fs_fitted;

  // se
  spdlog::debug("Estimating SE, thread = {}", omp_get_thread_num());
  double fs_sig2 = fs_resid.squaredNorm();
  long fs_df = X_complete1.rows() - X_complete1.cols();
  Eigen::MatrixXd fs_vcov = (X_complete1.transpose() * X_complete1).inverse();
  fs_se = (fs_vcov * (fs_sig2 / fs_df)).diagonal().cwiseSqrt();

  // pval
  spdlog::debug("Estimating P, thread = {}", omp_get_thread_num());
  t_stat = fs_fit(1, 0) / fs_se(1, 0);
  boost::math::students_t t_dist(fs_df);
  pval = 2.0 * boost::math::cdf(boost::math::complement(t_dist, fabs(t_stat)));
}
void Model::first_stage_lad(
    const Eigen::MatrixXd &X_complete1,
    const Eigen::VectorXd &y_complete,
    const Eigen::VectorXd &fs_fit,
    Eigen::VectorXd &fs_resid2,
    double &theta
) {
  spdlog::debug("Fitting quantile model, thread = {}", omp_get_thread_num());
  Eigen::VectorXd fs_fit_lad = CqrReg::QRMM::fit(X_complete1, y_complete, fs_fit, 0.001, 200, 0.5);
  Eigen::VectorXd fs_fitted = X_complete1 * fs_fit_lad;
  Eigen::VectorXd fs_resid = y_complete - fs_fitted;
  // estimate absolute residuals
  fs_resid2 = fs_resid.array().abs();
  theta = fs_fit_lad(1, 0);
}
void Model::second_stage(
    const Eigen::MatrixXd &X_complete2,
    const Eigen::VectorXd &fs_resid2,
    Eigen::VectorXd &ss_resid,
    Eigen::VectorXd &ss_fit
) {
  spdlog::debug("Estimating second stage model, thread = {}", omp_get_thread_num());
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr2(X_complete2);
  if (qr2.rank() < X_complete2.cols()) {
    throw std::runtime_error("Rank deficiency");
  }
  ss_fit = qr2.solve(fs_resid2);
  Eigen::VectorXd ss_fitted = X_complete2 * ss_fit;
  ss_resid = fs_resid2 - ss_fitted;
  Eigen::VectorXd ss_resid2 = ss_resid.array().square();
}
Eigen::MatrixXd Model::white_vcov(const Eigen::MatrixXd &X_complete2, const Eigen::VectorXd &ss_resid) {
  spdlog::debug("Estimating robust SE, thread = {}", omp_get_thread_num());
  Eigen::MatrixXd hc_vcov = (X_complete2.transpose() * X_complete2).inverse() * X_complete2.transpose()
      * ss_resid.cwiseProduct(ss_resid).asDiagonal() * X_complete2 * (X_complete2.transpose() * X_complete2).inverse();
  return (hc_vcov);
}
void Model::delta_method(const Eigen::VectorXd &ss_fit,
                         const Eigen::MatrixXd &hc_vcov,
                         double &s1_dummy,
                         double &s2_dummy) {
  spdlog::debug("Preparing partial derivatives for deltamethod, thread = {}", omp_get_thread_num());
  const double pi = boost::math::constants::pi<double>();

  Eigen::MatrixXd grad1 = Eigen::MatrixXd(ss_fit.rows(), 1);
  grad1(0, 0) = 0; // intercept
  grad1(1, 0) = (2 * ss_fit(0, 0) + 2 * ss_fit(1, 0)) / (2 / pi); // SNP=1
  grad1(2, 0) = 0; // SNP=2

  Eigen::MatrixXd grad2 = Eigen::MatrixXd(ss_fit.rows(), 1);
  grad2(0, 0) = 0; // intercept
  grad2(1, 0) = 0; // SNP=1
  grad2(2, 0) = (2 * ss_fit(0, 0) + 2 * ss_fit(2, 0)) / (2 / pi); // SNP=2

  // fill in covariates as zeros
  for (unsigned j = 3; j < ss_fit.rows(); j++) {
    grad1(j, 0) = 0;
    grad2(j, 0) = 0;
  }

  if (grad1.rows() != ss_fit.rows()) {
    throw std::runtime_error("Number of observations for G1 and X differ");
  }
  if (grad2.rows() != ss_fit.rows()) {
    throw std::runtime_error("Number of observations for G2 and X differ");
  }

  spdlog::debug("Estimating HC SE_1 using the deltamethod, thread = {}", omp_get_thread_num());
  Eigen::MatrixXd s1_hc0 = grad1.transpose() * hc_vcov * grad1;
  spdlog::debug("Estimating s1_dummy, thread = {}", omp_get_thread_num());
  s1_dummy = sqrt(s1_hc0(0, 0));
  spdlog::debug("Estimating HC SE_2 using the deltamethod, thread = {}", omp_get_thread_num());
  Eigen::MatrixXd s2_hc0 = grad2.transpose() * hc_vcov * grad2;
  spdlog::debug("Estimating s2_dummy, thread = {}", omp_get_thread_num());
  s2_dummy = sqrt(s2_hc0(0, 0));
}
void Model::null_fit(const Eigen::MatrixXd &X_complete2,
                     const Eigen::VectorXd &fs_resid2,
                     double &phi_pval,
                     double &f,
                     int n,
                     const Eigen::VectorXd &ss_resid) {
  spdlog::debug("Preparing null model, thread = {}", omp_get_thread_num());
  std::vector<unsigned> covariates_to_keep;
  for (unsigned i = 0; i < X_complete2.cols(); ++i) {
    // drop x==1 and x==2 but keep intercept & all covariates
    if (i == 1 || i == 2) {
      continue;
    } else {
      covariates_to_keep.push_back(i);
    }
  }
  Eigen::MatrixXd X_complete3 = X_complete2(Eigen::all, covariates_to_keep).eval();
  spdlog::debug("Fitting null model, thread = {}", omp_get_thread_num());
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr3(X_complete3);
  if (qr3.rank() < X_complete3.cols()) {
    throw std::runtime_error("Rank deficiency");
  }
  Eigen::VectorXd null_fit = qr3.solve(fs_resid2);
  spdlog::debug("Estimating variance of the null model, thread = {}", omp_get_thread_num());
  Eigen::VectorXd null_fitted = X_complete3 * null_fit;
  Eigen::VectorXd null_resid = fs_resid2 - null_fitted;
  Eigen::VectorXd null_resid2 = null_resid.array().square();
  // F-test
  // adapted from http://people.reed.edu/~jones/Courses/P24.pdf
  spdlog::debug("F-test, thread = {}", omp_get_thread_num());
  int df_f = n - X_complete2.cols();
  spdlog::debug("DF full model " + std::to_string(df_f) + ", thread = {}", omp_get_thread_num());
  int df_r = n - X_complete3.cols();
  spdlog::debug("DF restricted model " + std::to_string(df_r) + ", thread = {}", omp_get_thread_num());
  int df_n = df_r - df_f;
  double rss_f = ss_resid.squaredNorm();
  double rss_r = null_resid.squaredNorm();
  f = ((rss_r - rss_f) / (df_r - df_f)) / (rss_f / df_f);
  if (f < 0) {
    throw std::runtime_error("Invalid F-stat");
  }
  boost::math::fisher_f f_dist(df_n, df_f);
  phi_pval = boost::math::cdf(boost::math::complement(f_dist, f));
}

Result Model::fit(std::string &chromosome,
                  uint32_t position,
                  std::string &rsid,
                  std::string &effect_allele,
                  std::string &other_allele,
                  std::vector<double> dosages,
                  std::set<unsigned> non_null_idx,
                  Eigen::MatrixXd X1,
                  Eigen::MatrixXd X2,
                  Eigen::VectorXd y,
                  double maf_threshold) {
  Result res;
  res.chromosome = chromosome;
  res.position = position;
  res.rsid = rsid;
  res.effect_allele = effect_allele;
  res.other_allele = other_allele;
  res.eaf = NAN;
  res.n = -1;
  res.beta = NAN;
  res.theta = NAN;
  res.se = NAN;
  res.t = NAN;
  res.pval = NAN;
  res.phi_x1 = NAN;
  res.se_x1 = NAN;
  res.phi_x2 = NAN;
  res.se_x2 = NAN;
  res.phi_f = NAN;
  res.phi_pval = NAN;
  const double pi = boost::math::constants::pi<double>();

  // set dosage values
  // X is passed without reference to allow for modification on each thread
  spdlog::debug("Checking X has same number of observations, thread = {}", omp_get_thread_num());
  if (X1.rows() != X2.rows()) {
    throw std::runtime_error("Number of observations for X1 and X2 differ");
  }
  spdlog::debug("Checking for null dosage values, thread = {}", omp_get_thread_num());
  if (dosages.size() != X1.rows()) {
    throw std::runtime_error("Number of dosage values does not match number of participants");
  }
  spdlog::debug("Preparing matrices, thread = {}", omp_get_thread_num());
  int j = 0;
  for (unsigned i = 0; i < dosages.size(); i++) {
    if (dosages[i] == -1) {
      j++;
      non_null_idx.erase(i);
    } else {
      int g = int(round(dosages[i]));
      X1(i, 1) = dosages[i]; // dosage as continuous value
      if (g == 0) {
        X2(i, 1) = 0; // Dummy G==1
        X2(i, 2) = 0; // Dummy G==2
      } else if (g == 1) {
        X2(i, 1) = 1; // Dummy G==1
        X2(i, 2) = 0; // Dummy G==2
      } else if (g == 2) {
        X2(i, 1) = 0; // Dummy G==1
        X2(i, 2) = 1; // Dummy G==2
      } else {
        throw std::runtime_error("Rounded dosage value: " + std::to_string(g));
      }

    }
  }
  spdlog::debug("Found {} null values, thread = {}", j, omp_get_thread_num());

  // subset data with non-null values
  spdlog::debug("Selecting non-null values for model, thread = {}", omp_get_thread_num());
  std::vector<unsigned> non_null_idx_vec(non_null_idx.begin(), non_null_idx.end());
  Eigen::MatrixXd X_complete1 = X1(non_null_idx_vec, Eigen::all).eval();
  Eigen::MatrixXd X_complete2 = X2(non_null_idx_vec, Eigen::all).eval();
  Eigen::VectorXd y_complete = y(non_null_idx_vec, Eigen::all).eval();
  res.n = X_complete1.rows();

  // filter out low MAF
  res.eaf = X_complete1.col(1).mean() * 0.5;
  if (res.eaf < maf_threshold || (1 - res.eaf) < maf_threshold) {
    spdlog::warn("Skipping variant " + res.rsid + " with EAF of " + std::to_string(res.eaf) + ", thread = {}",
                 omp_get_thread_num());
    return res;
  }

  // first stage model ols
  Eigen::VectorXd fs_se;
  Eigen::VectorXd fs_fit;
  try {
    Model::first_stage_ols(X_complete1, y_complete, fs_se, fs_fit, res.pval, res.t);
    res.beta = fs_fit(1, 0);
    res.se = fs_se(1, 0);
  }
  catch (...) {
    return res;
  }

  // LAD regression
  Eigen::VectorXd fs_resid2;
  Model::first_stage_lad(X_complete1, y_complete, fs_fit, fs_resid2, res.theta);

  // second stage model
  Eigen::VectorXd ss_resid;
  Eigen::VectorXd ss_fit;
  try {
    Model::second_stage(X_complete2, fs_resid2, ss_resid, ss_fit);
  }
  catch (...) {
    return res;
  }

  // HC0 White SE for second-stage estimates
  Eigen::MatrixXd hc_vcov = Model::white_vcov(X_complete2, ss_resid);

  // transform MAD estimates to variance effects (assuming normality)
  spdlog::debug("Transforming MAD values to variance, thread = {}", omp_get_thread_num());
  res.phi_x1 = (2 * ss_fit(0, 0) * ss_fit(1, 0) + pow(ss_fit(1, 0), 2.0)) / (2 / pi);
  res.phi_x2 = (2 * ss_fit(0, 0) * ss_fit(2, 0) + pow(ss_fit(2, 0), 2.0)) / (2 / pi);

  // deltamethod to obtain the variance SE
  Model::delta_method(ss_fit, hc_vcov, res.se_x1, res.se_x2);

  // null model (intercept & covariates)
  try {
    Model::null_fit(X_complete2, fs_resid2, res.phi_pval, res.phi_f, res.n, ss_resid);
  }
  catch (...) {
    return res;
  }

  return res;
}

}
