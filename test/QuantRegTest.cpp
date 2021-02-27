#include "csv.h"
#include "libscl.h"
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