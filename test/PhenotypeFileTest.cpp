#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "gtest/gtest.h"

/*
 * Tests for reading phenotype file
 * */

// TODO update

TEST(PhenotypeFileTest, parse_should_function) {
  static std::string filePath = "data/regression.csv";
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
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][0], 32.4359168112278, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[0], 17.7335116415805, 1e-10);
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[49999], "S50000");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][49999], 0, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][49999], 52.4099412281066, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[49999], 35.0527037914899, 1e-10);
}

TEST(PhenotypeFileTest, parse_missing_exposure_field) {
  static std::string filePath = "data/regression.csv";
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
  static std::string filePath = "data/regression.csv";
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
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn().size(), 50000);
  phenotypeFile.subset_samples(std::vector<std::string>{"S1", "S4", "S100"});
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn().size(), 3);

  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[2], "S100");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][2], 0, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][2], 39.3593977019191, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[2], 32.1264679946589, 1e-10);
}

TEST(PhenotypeFileTest, subset_samples_missing_sample) {
  static std::string filePath = "data/regression.csv";
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
  try {
    phenotypeFile.subset_samples(std::vector<std::string>{"abc123"});
    FAIL() << "Expected std::runtime_error";
  } catch (std::runtime_error &err) {
    EXPECT_EQ(err.what(), std::string("Missing sample from phenotype file: abc123"));
  } catch (...) {
    FAIL() << "Expected std::runtime_error";
  }
}