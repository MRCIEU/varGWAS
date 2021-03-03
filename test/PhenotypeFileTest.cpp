#include <fstream>
#include <vector>
#include <stdexcept>
#include <set>
#include <algorithm>
#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "gtest/gtest.h"

/*
 * Tests for reading phenotype file
 * */

TEST(PhenotypeFileTest, parse_should_function) {
  static std::string filePath = "data.csv";
  static std::set<std::string> covariateColumnHeaders = {"c1", "c2"};
  static std::string outcomeColumnHeader = "y";
  static std::string idColumnHeader = "id";
  static char sep = ',';

  jlst::PhenotypeFile phenotypeFile(filePath,
                                    covariateColumnHeaders,
                                    outcomeColumnHeader,
                                    idColumnHeader,
                                    sep);
  phenotypeFile.parse();

  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[0], "S1");
  ASSERT_EQ(phenotypeFile.GetCovariateColumn().size(), 2);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][0], 0, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][0], 1.79407452864789, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[0], 2.44774533568743, 1e-10);
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[999], "S1000");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][999], 1, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][999], 0.692115714615254, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[999], 5.15096671268665, 1e-10);
}

TEST(PhenotypeFileTest, parse_missing_exposure_field) {
  static std::string filePath = "data.csv";
  static std::set<std::string> covariateColumnHeaders = {"c1", "c2"};
  static std::string outcomeColumnHeader = "NA";
  static std::string idColumnHeader = "id";
  static char sep = ',';
  try {
    jlst::PhenotypeFile phenotypeFile(filePath,
                                      covariateColumnHeaders,
                                      outcomeColumnHeader,
                                      idColumnHeader,
                                      sep);
    phenotypeFile.parse();
    FAIL() << "Expected PhenotypeFileException";
  } catch (jlst::PhenotypeFileException const &err) {
    EXPECT_EQ(err.what(), std::string("Field missing from phenotype file: NA"));
  } catch (...) {
    FAIL() << "Expected PhenotypeFileException";
  }
}

TEST(PhenotypeFileTest, subset_samples_should_function) {
  static std::string filePath = "data.csv";
  static std::set<std::string> covariateColumnHeaders = {"c1", "c2"};
  static std::string outcomeColumnHeader = "y";
  static std::string idColumnHeader = "id";
  static char sep = ',';

  jlst::PhenotypeFile phenotypeFile(filePath,
                                    covariateColumnHeaders,
                                    outcomeColumnHeader,
                                    idColumnHeader,
                                    sep);
  phenotypeFile.parse();
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn().size(), 1000);
  std::set<unsigned> non_null_idx = phenotypeFile.join(std::vector<std::string>{"S1", "S100", "S4"});
  ASSERT_EQ(non_null_idx.size(), 3);
  ASSERT_EQ(non_null_idx.count(0), 1);
  ASSERT_EQ(non_null_idx.count(1), 1);
  ASSERT_EQ(non_null_idx.count(2), 1);
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn().size(), 3);
  ASSERT_EQ(phenotypeFile.GetCovariateColumn()[0].size(), 3);
  ASSERT_EQ(phenotypeFile.GetCovariateColumn()[1].size(), 3);
  ASSERT_EQ(phenotypeFile.GetNSamples(), 3);

  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[1], "S100");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][1], 1, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][1], -0.921834963877281, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[1], 8.25428258247076, 1e-10);
}