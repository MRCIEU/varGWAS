#include <set>
#include <vector>
#include "gtest/gtest.h"
#include <Eigen/Core>

// Requires Eigen >=3.4
TEST(EigenTest, mask_matrix) {
  Eigen::MatrixXd X = Eigen::MatrixXd(5, 2);
  X << 1, 1, 2, 2, 3, 3, 4, 4, 5, 5;
  std::set<unsigned> idx_set = {0, 1, 2, 3};
  std::vector<unsigned> idx(idx_set.begin(), idx_set.end());
  Eigen::MatrixXd subset = X(idx, Eigen::all);
  ASSERT_EQ(subset.rows(), 4);
  ASSERT_EQ(subset.cols(), 2);
  ASSERT_EQ(subset(0, 0), 1);
  ASSERT_EQ(subset(1, 0), 2);
  ASSERT_EQ(subset(2, 0), 3);
}