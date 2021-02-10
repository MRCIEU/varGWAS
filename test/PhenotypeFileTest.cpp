#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "gtest/gtest.h"

/*
 * Tests for reading phenotype file
 * */

TEST(PhenotypeFileTest, parse) {
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

TEST(PhenotypeFileTest, missing_exposure_field) {
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
  } catch(jlst::PhenotypeFileException const &err) {
    EXPECT_EQ(err.what(),std::string("Field missing from phenotype file: NA"));
  } catch(...){
    FAIL() << "Expected PhenotypeFileException";
  }
}