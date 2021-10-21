#include "csv.h"
#include <Eigen/Core>
#include <Eigen/QR>
#include "gtest/gtest.h"
#include "cqrReg.h"

/*
 * Test for performing quantile regression
 * */

TEST(QuantRegTest, slope_residual) {
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
  Eigen::VectorXd betahat = qr.solve(y);

  // model
  double b1 = betahat(0, 0); // check that qrmm does not modify the vector which happens when pass by reference
  Eigen::VectorXd b = cqrReg::cqrReg::qrmm(X, y, betahat, 0.001, 200, 0.5);

  ASSERT_NEAR(b(0, 0), 4, 0.1);
  ASSERT_NEAR(b(1, 0), 0.6, 0.05);
  ASSERT_NEAR(b(2, 0), 2, 0.05);
  ASSERT_NEAR(b(3, 0), 0.3, 0.05);
  ASSERT_EQ(b1, betahat(0, 0));
}

TEST(QuantRegTest, slope_residual_outlier) {
  const double intercept = 1.0;
  double x_f;
  double y_f;
  int n = 100;
  int p = 2;

  Eigen::MatrixXd X = Eigen::MatrixXd(n, p);
  Eigen::VectorXd y = Eigen::VectorXd(n);

  // get data (see data/data.R)
  io::CSVReader<2> in("data-outlier.csv");
  in.read_header(io::ignore_extra_column, "x", "y");
  int t = 0;
  while (in.read_row(x_f, y_f)) {
    X(t, 0) = intercept;
    X(t, 1) = x_f;
    y(t, 0) = y_f;
    t++;
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
  Eigen::VectorXd betahat = qr.solve(y);

  // model
  Eigen::VectorXd b = cqrReg::cqrReg::qrmm(X, y, betahat, 0.001, 200, 0.5);

  ASSERT_NEAR(b(0, 0), 0, 0.1);
  ASSERT_NEAR(b(1, 0), 1, 0.1);
}