#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

#ifndef JLST_CPP_PHENOTYPEFILE_H
#define JLST_CPP_PHENOTYPEFILE_H

namespace jlst {
class PhenotypeFile {

 public:
  PhenotypeFile(
      const std::string &phenoFilePath,
      const std::string &exposureColumnHeader,
      const std::vector<std::string> &covariateColumnHeaders,
      const std::string &outcomeColumnHeader,
      const std::string &idColumnHeader,
      const char &sep
  );
  void parse();
  const std::vector<std::string> &GetSampleIdentifierColumn() const;
  const std::vector<double> &GetExposureColumn() const;
  const std::vector<double> &GetOutcomeColumn() const;
  const std::vector<std::vector<double>> &GetCovariateColumn() const;

 private:
  std::string phenoFilePath;
  std::vector<std::string> covariateColumnHeaders;
  std::string exposureColumnHeader;
  std::string outcomeColumnHeader;
  std::string idColumnHeader;
  std::vector<std::string> sampleIdentifierColumn;
  std::vector<double> exposureColumn;
  std::vector<double> outcomeColumn;
  std::vector<std::vector<double>> covariateColumn;
  char sep;

};
}

#endif //JLST_CPP_PHENOTYPEFILE_H
