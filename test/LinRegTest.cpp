#include "csv.h"
#include "gtest/gtest.h"
#include "iostream"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Dense>
#include <boost/math/distributions/students_t.hpp>

/*
 * Test for performing linear regression model
 * */

// Helpful notes https://web.stanford.edu/~mrosenfe/soc_meth_proj3/matrix_OLS_NYU_notes.pdf
// TODO random effect for pop strat https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf

TEST(LinRegTest, svd) {
  const double intercept = 1.0;
  double x_f;
  double c1_f;
  double c2_f;
  double y_f;
  int n = 1000;
  int p = 3;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p + 1);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/linreg.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  int t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    X(t, 2) = c1_f;
    X(t, 3) = c2_f;
    y(t, 0) = y_f;
    t++;
  }

  // linear regression using SVD
  // adapted from: https://genome.sph.umich.edu/w/images/2/2c/Biostat615-lecture14-presentation.pdf
  Eigen::BDCSVD<Eigen::MatrixXd> solver(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // beta
  Eigen::MatrixXd betas = solver.solve(y); // betas
  assert(betas.size() == p + 1);
  ASSERT_NEAR(betas(0, 0), 4, 4 * .2);
  ASSERT_NEAR(betas(1, 0), 0.6, 0.6 * .2);
  ASSERT_NEAR(betas(2, 0), 2, 2 * 0.2);
  ASSERT_NEAR(betas(3, 0), 0.3, .3 * .2);

  // predicted Y
  Eigen::VectorXd y_hat = X * betas;
  assert(y_hat.size() == n);

  // residuals
  Eigen::VectorXd y_delta = y - y_hat;
  assert(y_delta.size() == n);

  // squared residuals
  Eigen::VectorXd y_deltasq = y_delta.cwiseProduct(y_delta);
  assert(y_deltasq.size() == n);
  assert(y_deltasq[0] == y_delta[0] * y_delta[0]);
  assert(y_deltasq[100] == y_delta[100] * y_delta[100]);

  // unbiased variance of error term
  double e_var = (y - X * betas).squaredNorm() / (n - p);

  // se
  Eigen::MatrixXd ViD = solver.matrixV() * solver.singularValues().asDiagonal().inverse();
  Eigen::MatrixXd varBetasSvd = e_var * ViD * ViD.transpose();
  Eigen::VectorXd se = varBetasSvd.diagonal().array().sqrt();
  assert(se.size() == p + 1);
  ASSERT_NEAR(se(0, 0), 0.0557, 0.0557 * .2);
  ASSERT_NEAR(se(1, 0), 0.0314, 0.0314 * .2);
  ASSERT_NEAR(se(2, 0), 0.0451, 0.0451 * .2);
  ASSERT_NEAR(se(3, 0), 0.0303, 0.0303 * .2);

  // t-stat
  Eigen::VectorXd tstat = betas.array() / se.array();

  // pval
  boost::math::students_t dist(n - (p + 1)); // use student's t-distribution to compute p-value
  std::vector<double> pvalues;
  for (int i = 0; i < tstat.size(); i++) {
    double pval = 2.0 * cdf(complement(dist, tstat[i] > 0 ? tstat[i] : (0 - tstat[i])));
    pvalues.push_back(pval);
  }
  assert(pvalues.size() == p + 1);
  ASSERT_NEAR(pvalues[0], 0, 5e-5);
  ASSERT_NEAR(pvalues[1], 7.81e-68, 5e-5);
  ASSERT_NEAR(pvalues[2], 4.58e-233, 5e-5);
  ASSERT_NEAR(pvalues[3], 2.39e-17, 5e-5);
}

TEST(LinRegTest, normal_eq) {
  const double intercept = 1.0;
  double x_f;
  double c1_f;
  double c2_f;
  double y_f;
  int n = 1000;
  int p = 3;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p + 1);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/data.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  int t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    X(t, 2) = c1_f;
    X(t, 3) = c2_f;
    y(t, 0) = y_f;
    t++;
  }

  // linear regression using normal equations
  Eigen::MatrixXd betas = (X.transpose() * X).ldlt().solve(X.transpose() * y);

  // beta
  assert(betas.size() == p + 1);
  ASSERT_NEAR(betas(0, 0), 4, 4 * .2);
  ASSERT_NEAR(betas(1, 0), 0.6, 0.6 * .2);
  ASSERT_NEAR(betas(2, 0), 2, 2 * 0.2);
  ASSERT_NEAR(betas(3, 0), 0.3, .3 * .2);

  // predicted Y
  Eigen::VectorXd y_hat = X * betas;
  assert(y_hat.size() == n);

  // residuals
  Eigen::VectorXd y_delta = y - y_hat;
  assert(y_delta.size() == n);

  // unbiased variance of error term
  double y_var = (y - X * betas).squaredNorm() / (n - p);

  // https://stats.stackexchange.com/questions/236437/how-to-compute-the-standard-error-of-a-predictor-variable
  Eigen::VectorXd ssq = (X.rowwise() - (X.colwise().sum() / X.rows()).eval()).colwise().squaredNorm();
  std::vector<double> se;
  for (int i = 0; i < p + 1; i++) {
    double std_error = sqrt(y_var / ssq[i]);
    se.push_back(std_error);
    //std::cout << std_error << std::endl;
  }

  assert(se.size() == p + 1);
  ASSERT_NEAR(se[1], 0.03136, 0.03136 * .2);
  ASSERT_NEAR(se[2], 0.04509, 0.04509 * 0.2);
  ASSERT_NEAR(se[3], 0.03030, 0.03030 * .2);
}

TEST(LinRegTest, ColPivHouseholderQR) {
  const double intercept = 1.0;
  double x_f;
  double c1_f;
  double c2_f;
  double y_f;
  int n = 1000;
  int p = 4;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/data.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  int t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    X(t, 2) = c1_f;
    X(t, 3) = c2_f;
    y(t, 0) = y_f;
    t++;
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
  Eigen::MatrixXd betahat = qr.solve(y);

  if (qr.rank() < X.cols()) {
    throw std::runtime_error("rank-deficient matrix");
  }

  Eigen::VectorXd fitted = X * betahat;
  Eigen::VectorXd resid = y - fitted;

  // unbiased variance of error term
  double e = resid.squaredNorm() / (n - p);

  // TODO ask why cov of intercept is not 0
  // https://stats.stackexchange.com/questions/236437/how-to-compute-the-standard-error-of-a-predictor-variable
  Eigen::VectorXd ssq = (X.rowwise() - (X.colwise().sum() / X.rows()).eval()).colwise().squaredNorm();
  std::vector<double> se;
  for (int i = 0; i < p; i++) {
    double std_error = sqrt(e / ssq[i]);
    se.push_back(std_error);
  }

  // beta
  assert(betahat.size() == p);
  ASSERT_NEAR(betahat(0, 0), 4, 4 * .2);
  ASSERT_NEAR(betahat(1, 0), 0.6, 0.6 * .2);
  ASSERT_NEAR(betahat(2, 0), 2, 2 * 0.2);
  ASSERT_NEAR(betahat(3, 0), 0.3, .3 * .2);

  // se
  assert(se.size() == p);
  ASSERT_NEAR(se[1], 0.03136, 0.03136 * .2);
  ASSERT_NEAR(se[2], 0.04509, 0.04509 * 0.2);
  ASSERT_NEAR(se[3], 0.03030, 0.03030 * .2);
}

TEST(LinRegTest, vcov_se) {
  const double intercept = 1.0;
  double x_f;
  double c1_f;
  double c2_f;
  double y_f;
  int n = 1000;
  int p = 4;
  int df = n - p;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/data.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  int t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    X(t, 2) = c1_f;
    X(t, 3) = c2_f;
    y(t, 0) = y_f;
    t++;
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
  Eigen::MatrixXd betahat = qr.solve(y);

  if (qr.rank() < X.cols()) {
    throw std::runtime_error("rank-deficient matrix");
  }

  Eigen::VectorXd fitted = X * betahat;
  Eigen::VectorXd resid = y - fitted;

  // variance of error term
  double sig2 = resid.squaredNorm();

  // taken from https://stats.stackexchange.com/a/65444
  Eigen::MatrixXd vcov = (X.transpose() * X).inverse();
  Eigen::VectorXd se = (vcov * (sig2 / df)).diagonal().cwiseSqrt();

  // beta
  assert(betahat.size() == p);
  ASSERT_NEAR(betahat(0, 0), 4, 4 * .15);
  ASSERT_NEAR(betahat(1, 0), 0.6, 0.6 * .15);
  ASSERT_NEAR(betahat(2, 0), 2, 2 * .15);
  ASSERT_NEAR(betahat(3, 0), 0.3, .3 * .15);

  // se
  assert(se.size() == p);
  ASSERT_NEAR(se(0, 0), 0.05570, 0.05570 * .1);
  ASSERT_NEAR(se(1, 0), 0.03136, 0.03136 * .1);
  ASSERT_NEAR(se(2, 0), 0.04509, 0.04509 * .1);
  ASSERT_NEAR(se(3, 0), 0.03030, 0.03030 * .1);
}

TEST(LinRegTest, check_se_with_r) {
  const double intercept = 1.0;
  double x_f;
  double c1_f;
  double c2_f;
  double y_f;
  int n = 1000;
  int p = 4;
  int df = n - p;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/data.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  int t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    X(t, 2) = c1_f;
    X(t, 3) = c2_f;
    y(t, 0) = y_f;
    t++;
  }

  Eigen::MatrixXd betahat = Eigen::MatrixXd(4, 1);
  betahat(0, 0) = 4.06302;
  betahat(1, 0) = 0.59024;
  betahat(2, 0) = 1.96525;
  betahat(3, 0) = 0.26150;

  Eigen::VectorXd fitted = X * betahat;
  Eigen::VectorXd resid = y - fitted;

  // variance of error term
  double sig2 = resid.squaredNorm();

  // taken from https://stats.stackexchange.com/a/65444
  Eigen::MatrixXd vcov = (X.transpose() * X).inverse();
  Eigen::VectorXd se = (vcov * (sig2 / df)).diagonal().cwiseSqrt();

  // se
  assert(se.size() == p);
  ASSERT_NEAR(se(0, 0), 0.05570, 0.05570 * .001);
  ASSERT_NEAR(se(1, 0), 0.03136, 0.03136 * .001);
  ASSERT_NEAR(se(2, 0), 0.04509, 0.04509 * .001);
  ASSERT_NEAR(se(3, 0), 0.03030, 0.03030 * .001);
}