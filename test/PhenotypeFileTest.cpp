#include "PhenotypeFile.h"
#include "gtest/gtest.h"

TEST(PhenotypeFileTest, load) {
  std::string filename = "data/regression.csv";
  std::string outcomeColumnHeader = "y";
  std::vector<std::string> covariateColumnHeaders = {"c1", "c2"};
  char sep = ',';

  // read file
  jlst::PhenotypeFile p(filename, covariateColumnHeaders, outcomeColumnHeader, sep);
  p.load();

  // assertions
  ASSERT_EQ(p.GetFileHeader()[0], "x");
  ASSERT_EQ(p.GetFileHeader()[1], "c1");
  ASSERT_EQ(p.GetFileHeader()[2], "c2");
  ASSERT_EQ(p.GetFileHeader()[3], "y");
  ASSERT_EQ(p.GetFileHeader()[4], "d");
  ASSERT_NEAR(p.GetFileBody()[0][0], (double) 1, 0.001);
  ASSERT_NEAR(p.GetFileBody()[0][1], (double) 0, 0.001);
  ASSERT_NEAR(p.GetFileBody()[0][2], (double) 32.435, 0.001);
  ASSERT_NEAR(p.GetFileBody()[0][3], (double) 17.733, 0.001);
  ASSERT_NEAR(p.GetFileBody()[0][4], (double) 9.527, 0.001);
}