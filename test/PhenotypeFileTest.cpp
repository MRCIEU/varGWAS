#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "gtest/gtest.h"

/*
 * Tests for reading phenotype file
 * */

TEST(PhenotypeFileTest, parse_should_function) {
  static std::string filePath = "data.csv";
  static std::vector<std::string> covariateColumnHeaders = {"c1", "c2"};
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
  static std::vector<std::string> covariateColumnHeaders = {"c1", "c2"};
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
  static std::vector<std::string> covariateColumnHeaders = {"c1", "c2"};
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
  phenotypeFile.subset_samples(std::vector<std::string>{"S1", "S4", "S100"});
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn().size(), 3);

  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[2], "S100");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][2], 1, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][2], -0.921834963877281, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[2], 8.25428258247076, 1e-10);
}