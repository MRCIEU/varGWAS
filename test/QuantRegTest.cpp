#include "csv.h"
#include "libscl.h"
#include "gtest/gtest.h"
#include <cmath>

/*
 * Test for performing quantile regression
 * */

TEST(QuantRegTest, slope_residual) {
  INTEGER x_f;
  INTEGER c1_f;
  REAL c2_f;
  REAL y_f;
  REAL d_f;

  INTEGER n = 50000;
  INTEGER K = 4;

  scl::realmat y(n, 1);
  scl::realmat X(n, K);
  scl::realmat b(K, 1);
  scl::realmat d(n, 1);

  // get data (see data/regression.R)
  io::CSVReader<5> in("data/regression.csv");
  in.read_header(io::ignore_extra_column, "x", "c1", "c2", "y", "d");
  INTEGER t = 0;
  while (in.read_row(x_f, c1_f, c2_f, y_f, d_f)) {
    t++;
    // intercept
    X(t, 1) = 1.0;

    // slope
    X(t, 2) = x_f;
    X(t, 3) = c1_f;
    X(t, 4) = c2_f;
    y[t] = y_f;

    // absolute residuals
    d[t] = d_f;
  }

  // model
  REAL p = 0.5;
  b = scl::quantreg(y, X, p);

  ASSERT_NEAR(b[1], 25, 0.1);
  ASSERT_NEAR(b[2], 0.6, 0.1);
  ASSERT_NEAR(b[3], 2, 0.1);
  ASSERT_NEAR(b[4], 0.05, 0.002);

  // absolute residual assertions
  scl::realmat e = y - X * b;
  for (t = 1; t <= n; ++t) {
    ASSERT_NEAR(std::abs(e[t]), d[t], 0.005);
  }
}