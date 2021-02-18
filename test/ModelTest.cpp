#include "gtest/gtest.h"
#include "iostream"
#include "Model.h"
#include "Result.h"

/*
 * Test for performing Breusch-Pagan model
 * */

TEST(ModelTest, fit) {
  jlst::Result result;
  std::vector<double> dosages = {

  };
  std::vector<double> pheno = {

  };
  Eigen::MatrixXd X;
  Eigen::VectorXd y;

  // initialise empty matrix
  assert(dosages.size() == pheno.size());
  for (unsigned n = 0; n < dosages.size(); ++n) {
    X(n, 0) = 1;
    X(n, 1) = 0;
    y(n, 0) = pheno[n];
  }

  // fit B-P model
  jlst::Model::fit(result, dosages, X, y);

  // print summary stats
  std::cout << result.beta << std::endl;
  std::cout << result.se << std::endl;
  std::cout << result.pval << std::endl;

  // check output
  ASSERT_NEAR(result.beta, 0, 0.01);
  ASSERT_NEAR(result.se, 0, 0.01);
  ASSERT_NEAR(result.pval, 0, 0.01);
}