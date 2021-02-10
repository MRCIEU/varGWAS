#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

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
  void GetMatrix(std::vector<std::string>samples);
  const std::vector<std::string> &GetCovariateColumnHeaders() const;
  const std::string &GetOutcomeColumnHeader() const;
  const std::vector<std::string> &GetFileHeader() const;
  const std::vector<std::vector<double>> &GetFileBody() const;

 private:
  std::string phenoFilePath;
  std::vector<std::string> covariateColumnHeaders;
  std::string outcomeColumnHeader;
  std::string idColumnHeader;
  char sep;
  std::vector<std::string> fileHeader;
  std::vector<std::vector<double>> fileBody;

};
}

#endif //JLST_CPP_PHENOTYPEFILE_H
