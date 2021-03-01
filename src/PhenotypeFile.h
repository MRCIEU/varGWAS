#include <string>
#include <vector>
#include <set>

#ifndef JLST_CPP_PHENOTYPEFILE_H
#define JLST_CPP_PHENOTYPEFILE_H

namespace jlst {
class PhenotypeFile {

 public:
  PhenotypeFile(
      std::string &pheno_file_path,
      std::vector<std::string> &covariate_column_headers,
      std::string &outcome_column_header,
      std::string &id_column_header,
      char sep)
      : _pheno_file_path(pheno_file_path),
        _covariate_column_headers(covariate_column_headers),
        _outcome_column_header(outcome_column_header),
        _id_column_header(id_column_header),
        _sep(sep),
        _n_samples(-1) {}
  void parse();
  std::set<unsigned> join(const std::vector<std::string> &samples);
  const std::vector<std::string> &GetSampleIdentifierColumn() const;
  const std::vector<double> &GetOutcomeColumn() const;
  const std::vector<std::vector<double>> &GetCovariateColumn() const;
  int GetNSamples() const;

 private:
  std::string &_pheno_file_path;
  std::vector<std::string> &_covariate_column_headers;
  std::string &_outcome_column_header;
  std::string &_id_column_header;
  std::vector<std::string> _sample_identifier_column;
  std::vector<double> _outcome_column;
  std::vector<std::vector<double>> _covariate_column;
  char _sep;
  int _n_samples;
};
}

#endif //JLST_CPP_PHENOTYPEFILE_H
