#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <unordered_map>
#include <algorithm>
#include "spdlog/spdlog.h"
#include "PhenotypeFileException.h"
#include "PhenotypeFile.h"

/*
 * Class to read in outcome, covariate(s) and sample identifier into memory
 * */

namespace jlst {
/*
 * Function to parse file
 * */
void PhenotypeFile::parse() {
  spdlog::info("Parsing phenotype from: {}", _pheno_file_path);
  for (auto &c: _covariate_column_headers) {
    spdlog::info("Including covariate: {}", c);
  }
  spdlog::info("Outcome variable: {}", _outcome_column_header);
  std::ifstream file(_pheno_file_path.c_str());
  int out_idx = -1;
  int sid_idx = -1;
  std::vector<int> cov_idx;
  int i = -1;
  int cov_offset = -1;

  if (file.is_open()) {
    bool passed_first_line = false;
    std::string line;
    std::string token;

    // read file line-by-line
    while (getline(file, line)) {
      spdlog::trace(line);
      std::istringstream token_stream(line);

      if (passed_first_line) { // read file body
        i = 0;

        // TODO implement using boost to allow for quotes in the file
        while (std::getline(token_stream, token, _sep)) {
          spdlog::trace("token={}, value={}", i, token);

          // is the current token from the outcome column
          if (i == out_idx) {
            try {
              spdlog::trace("outcome value={}", token);
              long double val = std::stold(token);
              _outcome_column.push_back(val);
            } catch (...) {
              throw std::runtime_error("Could not cast outcome value to numeric: " + token);
            }
          }

          // is the current token from the sample column
          if (i == sid_idx) {
            spdlog::trace("sample ID value={}", token);
            _sample_identifier_column.push_back(token);
          }

          // is the current token from a covariate column
          if (cov_offset > -1 && std::find(cov_idx.begin(), cov_idx.end(), i) != cov_idx.end()) {
            int idx = i - cov_offset;
            assert(_covariate_column.size() > idx);
            try {
              spdlog::trace("covariate n={}, value={}", idx, token);
              long double val = std::stold(token);
              _covariate_column[idx].push_back(val);
            } catch (...) {
              throw std::runtime_error("Could not cast covariate value to numeric: " + token);
            }
          }

          i++;
        }

      } else { // read file header
        i = 0;

        while (std::getline(token_stream, token, _sep)) {
          spdlog::debug("Phenotype file header n={}: {}", i, token);

          // record file column numbers of model variables
          if (std::find(_covariate_column_headers.begin(), _covariate_column_headers.end(), token)
              != _covariate_column_headers.end()) {
            cov_idx.push_back(i);
            _covariate_column.emplace_back(); // instantiate v of v
            spdlog::debug("Found covariate index: {}", i);
          } else if (token == _outcome_column_header) {
            out_idx = i;
            spdlog::debug("Found outcome index: {}", out_idx);
          } else if (token == _id_column_header) {
            sid_idx = i;
            spdlog::debug("Found sample id index: {}", sid_idx);
          }

          i++;
        }

        if (out_idx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + _outcome_column_header);
        }
        if (sid_idx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + _id_column_header);
        }
        if (_covariate_column.size() != _covariate_column_headers.size()
            || cov_idx.size() != _covariate_column.size()) {
          throw jlst::PhenotypeFileException(
              "Expected n=" + std::to_string(_covariate_column_headers.size()) + " covariate(s) but n="
                  + std::to_string(_covariate_column.size()) + " were found in the file");
        }
        passed_first_line = true;

        // get minimum covariate index in vector
        for (int idx : cov_idx) {
          if (cov_offset == -1 || idx < cov_offset) {
            cov_offset = idx;
          }
        }
        if (!_covariate_column_headers.empty()) {
          assert(cov_offset != -1);
        }
      }

    }

    file.close();
  } else {
    throw std::runtime_error("Could not open file: " + _pheno_file_path);
  }
  _n_samples = _sample_identifier_column.size();
  spdlog::info("Found {} samples in phenotype file", _n_samples);
}
/*
 * Function to subset data using provided list of sample identifiers
 * @return vector of sample indexes with missing phenotypes to mask in the model
 * */
std::set<unsigned> PhenotypeFile::join(const std::vector<std::string> &samples) {
  std::set<unsigned> non_null_idx;
  spdlog::info("Joining phenotypes with {} samples", samples.size());

  // create sample identifier-to-index mapping
  std::unordered_map<std::string, unsigned> mapping;
  for (unsigned i = 0; i < _sample_identifier_column.size(); ++i) {
    if (mapping.count(_sample_identifier_column[i])) {
      throw std::runtime_error("Duplicate identifiers were detected for: " + _sample_identifier_column[i]);
    }
    // store array index against column identifier
    mapping[_sample_identifier_column[i]] = i;
  }

  // create tmp data stores
  std::vector<std::string> sample_identifier_column_tmp;
  std::vector<long double> outcome_column_tmp;
  std::vector<std::vector<long double>> covariate_column_tmp;
  for (unsigned i = 0; i < _covariate_column.size(); ++i) {
    covariate_column_tmp.emplace_back();
  }

  // subset data using new ordering
  for (unsigned i = 0; i < samples.size(); ++i) {

    // get index for sample
    if (mapping.count(samples[i]) == 0) {
      spdlog::info("Sample missing from phenotype file: {}", samples[i]);
      // add missing data with null values
      sample_identifier_column_tmp.push_back(samples[i]);
      outcome_column_tmp.push_back(-1);
      for (unsigned j = 0; j < _covariate_column.size(); ++j) {
        covariate_column_tmp[j].push_back(-1);
      }
    } else {
      // subset data
      unsigned idx = mapping[samples[i]];
      sample_identifier_column_tmp.push_back(_sample_identifier_column[idx]);
      outcome_column_tmp.push_back(_outcome_column[idx]);
      for (unsigned j = 0; j < _covariate_column.size(); ++j) {
        covariate_column_tmp[j].push_back(_covariate_column[j][idx]);
      }
      non_null_idx.insert(i);
    }

  }

  // check variables are the same length
  assert(sample_identifier_column_tmp.size() == samples.size());
  assert(outcome_column_tmp.size() == samples.size());
  for (unsigned i = 0; i < _covariate_column.size(); ++i) {
    assert(covariate_column_tmp[i].size() == samples.size());
  }

  // set variables
  _sample_identifier_column = sample_identifier_column_tmp;
  _outcome_column = outcome_column_tmp;
  _covariate_column = covariate_column_tmp;
  _n_samples = _sample_identifier_column.size();

  spdlog::info("Remaining samples after join: {}", _n_samples);
  spdlog::info("Samples with non-null phenotypes: {}", non_null_idx.size());

  return (non_null_idx);
}
const std::vector<std::string> &PhenotypeFile::GetSampleIdentifierColumn() const {
  return _sample_identifier_column;
}
const std::vector<long double> &PhenotypeFile::GetOutcomeColumn() const {
  return _outcome_column;
}
const std::vector<std::vector<long double>> &PhenotypeFile::GetCovariateColumn() const {
  return _covariate_column;
}
int PhenotypeFile::GetNSamples() const {
  if (_n_samples == -1) {
    throw std::runtime_error("Value not set");
  }
  return _n_samples;
}
}