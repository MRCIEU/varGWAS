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
      const std::vector<std::string> &covariateColumnHeaders,
      const std::string &outcomeColumnHeader,
      const std::string &idColumnHeader,
      const char &sep
  );
  void parse();
  std::vector<unsigned> subset_samples(const std::vector<std::string> &samples);
  const std::vector<std::string> &GetSampleIdentifierColumn() const;
  const std::vector<double> &GetOutcomeColumn() const;
  const std::vector<std::vector<double>> &GetCovariateColumn() const;
  int GetNSamples() const;

 private:
  std::string phenoFilePath;
  std::vector<std::string> covariateColumnHeaders;
  std::string outcomeColumnHeader;
  std::string idColumnHeader;
  std::vector<std::string> sampleIdentifierColumn;
  std::vector<double> outcomeColumn;
  std::vector<std::vector<double>> covariateColumn;
  char sep;
  int n_samples;
};
}

#endif //JLST_CPP_PHENOTYPEFILE_H
