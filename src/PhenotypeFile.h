#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <exception>

#ifndef JLST_CPP_PHENOTYPEFILE_H
#define JLST_CPP_PHENOTYPEFILE_H

namespace jlst {
class PhenotypeFile {

 public:
  PhenotypeFile(
      const std::string &pheno_file_path,
      const std::vector<std::string> &covariate_column_headers,
      const std::string &outcome_column_header,
      const std::string &id_column_header,
      const char &sep
  );
  void parse();
  std::set<unsigned> join(const std::vector<std::string> &samples);
  const std::vector<std::string> &GetSampleIdentifierColumn() const;
  const std::vector<double> &GetOutcomeColumn() const;
  const std::vector<std::vector<double>> &GetCovariateColumn() const;
  int GetNSamples() const;

 private:
  std::string pheno_file_path;
  std::vector<std::string> covariate_column_headers;
  std::string outcome_column_header;
  std::string id_column_header;
  std::vector<std::string> sample_identifier_column;
  std::vector<double> outcome_column;
  std::vector<std::vector<double>> covariate_column;
  char sep;
  int n_samples;
};
}

#endif //JLST_CPP_PHENOTYPEFILE_H
