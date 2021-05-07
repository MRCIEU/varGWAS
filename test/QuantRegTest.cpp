#include "csv.h"
#include "libscl.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <Eigen/Dense>
#include "gtest/gtest.h"

/*
 * Test for performing quantile regression
 * */

TEST(QuantRegTest, slope_residual) {
  REAL x_f;
  REAL c1_f;
  REAL c2_f;
  REAL y_f;

  INTEGER n = 1000;
  INTEGER K = 4;

  scl::realmat y(n, 1);
  scl::realmat X(n, K);
  scl::realmat b(K, 1);
  scl::realmat d(n, 1);

  // get data (see data/data.R)
  io::CSVReader<4> in("data.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y");
  INTEGER t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f)) {
    t++;
    // intercept
    X(t, 1) = 1.0;

    // slope
    X(t, 2) = x_f;
    X(t, 3) = c1_f;
    X(t, 4) = c2_f;
    y[t] = y_f;
  }

  // model
  REAL p = 0.5;
  b = scl::quantreg(y, X, p);

  ASSERT_NEAR(b[1], 4, 4 * .2);
  ASSERT_NEAR(b[2], 0.6, 0.6 * .2);
  ASSERT_NEAR(b[3], 2, 2 * 0.2);
  ASSERT_NEAR(b[4], 0.3, .3 * .2);
}

TEST(QuantRegTest, eigen2realmat) {
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

  scl::realmat scl_X(X.rows(), X.cols());
  scl::realmat scl_y(X.rows(), 1);
  for (int i = 0; i < X.rows(); i++) {
    scl_y[i + 1] = y(i, 0);
    for (int j = 0; j < X.cols(); j++) {
      scl_X(i + 1, j + 1) = X(i, j);
    }
  }

  ASSERT_EQ(X(0, 0), scl_X(1, 1));
  ASSERT_EQ(X(0, 1), scl_X(1, 2));
  ASSERT_EQ(X(0, 2), scl_X(1, 3));
  ASSERT_EQ(X(0, 3), scl_X(1, 4));
  ASSERT_EQ(y(0, 0), scl_y(1, 1));

  ASSERT_EQ(X(X.rows() - 1, 0), scl_X(scl_X.get_rows(), 1));
  ASSERT_EQ(X(X.rows() - 1, 1), scl_X(scl_X.get_rows(), 2));
  ASSERT_EQ(X(X.rows() - 1, 2), scl_X(scl_X.get_rows(), 3));
  ASSERT_EQ(X(X.rows() - 1, 3), scl_X(scl_X.get_rows(), 4));
  ASSERT_EQ(y(X.rows() - 1, 0), scl_y(scl_y.get_rows(), 1));

  // model
  scl::realmat b(4, 1);
  REAL tau = 0.5;
  b = scl::quantreg(scl_y, scl_X, tau);

  ASSERT_NEAR(b[1], 4, 4 * .2);
  ASSERT_NEAR(b[2], 0.6, 0.6 * .2);
  ASSERT_NEAR(b[3], 2, 2 * 0.2);
  ASSERT_NEAR(b[4], 0.3, .3 * .2);
}