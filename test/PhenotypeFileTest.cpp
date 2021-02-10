#include "PhenotypeFile.h"
#include "gtest/gtest.h"

/*
 * Test for reading phenotype file
 * */

TEST(PhenotypeFileTest, load) {
  std::string filePath = "data/regression.csv";
  std::string exposureColumnHeader = "x";
  std::vector<std::string> covariateColumnHeaders = {"c1", "c2"};
  std::string outcomeColumnHeader = "y";
  std::string idColumnHeader = "id";
  char sep = ',';

  // parse phenotype file
  jlst::PhenotypeFile phenotypeFile(filePath,
                        covariateColumnHeaders,
                        outcomeColumnHeader,
                        idColumnHeader,
                        sep);
  phenotypeFile.parse();

  // check values
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[0], "S1");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][0], 0, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][0], 32.4359168112278, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[0], 17.7335116415805, 1e-10);
  ASSERT_EQ(phenotypeFile.GetSampleIdentifierColumn()[49999], "S50000");
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[0][49999], 0, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetCovariateColumn()[1][49999], 52.4099412281066, 1e-10);
  ASSERT_NEAR(phenotypeFile.GetOutcomeColumn()[49999], 35.0527037914899, 1e-10);
}