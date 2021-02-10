#include <iostream>
#include "Result.h"
#include "Model.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include "gtest/gtest.h"

/*
 * Test for performing variance model
 * */

TEST(ModelTest, fit_should_function) {
  Eigen::MatrixXd X = Eigen::MatrixXd(3, 2);
  Eigen::VectorXd y = Eigen::VectorXd(3);

  // intercept
  X(0, 0) = 1;
  X(1, 0) = 1;
  X(2, 0) = 1;

  // dosage (X)
  X(0, 1) = 0;
  X(1, 1) = 0;
  X(2, 1) = 0;

  // Y
  y(0, 0) = 0;
  y(1, 0) = 5;
  y(2, 0) = 10;

  jlst::Result result = jlst::Model::fit("22",
                             1000000,
                             "rs1234",
                             std::vector<std::string>{"A", "G"},
                             std::vector<double>{0, 1, 2},
                             X,
                             y);

  ASSERT_EQ(result.chromosome, "22");
  ASSERT_EQ(result.position, 1000000);
  ASSERT_EQ(result.rsid, "rs1234");
  ASSERT_EQ(result.effect_allele, "G");
  ASSERT_EQ(result.other_allele, "A");
  ASSERT_EQ(result.beta, 5);
}