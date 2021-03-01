#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <unordered_map>
#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "spdlog/spdlog.h"


/*
 * Class to read in outcome, covariate(s) and sample identifier into memory
 * */

namespace jlst {
/*
 * Class for working with phenotype files
 * */
PhenotypeFile::PhenotypeFile(const std::string &pheno_file_path,
                             const std::vector<std::string> &covariate_column_headers,
                             const std::string &outcome_column_header,
                             const std::string &id_column_header,
                             const char &sep) {
  this->pheno_file_path = pheno_file_path;
  this->covariate_column_headers = covariate_column_headers;
  this->outcome_column_header = outcome_column_header;
  this->id_column_header = id_column_header;
  this->sep = sep;
  this->n_samples = -1;
};

/*
 * Function to parse file
 * */
void PhenotypeFile::parse() {
  spdlog::info("Parsing phenotype from: {}", pheno_file_path);
  for (auto &c: covariate_column_headers) {
    spdlog::info("Including covariate: {}", c);
  }
  spdlog::info("Outcome variable: {}", outcome_column_header);
  std::ifstream file(pheno_file_path.c_str());
  int out_idx = -1;
  int sid_idx = -1;
  std::vector<int> cov_idx;
  int i;
  unsigned min_cov_idx = UINT_MAX;

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
        while (std::getline(token_stream, token, sep)) {
          spdlog::trace("token={}, value={}", i, token);
          if (i == out_idx) {
            try {
              spdlog::trace("outcome value={}", token);
              outcome_column.push_back(std::stod(token));
            } catch (...) {
              throw std::runtime_error("Could not cast outcome value to numeric: " + token);
            }
          }
          if (i == sid_idx) {
            spdlog::trace("sample ID value={}", token);
            sample_identifier_column.push_back(token);
          }
          if (std::find(cov_idx.begin(), cov_idx.end(), i) != cov_idx.end()) {
            try {
              spdlog::trace("covariate n={}, value={}", i - min_cov_idx, token);
              covariate_column[i - min_cov_idx].push_back(std::stod(token));
            } catch (...) {
              throw std::runtime_error("Could not cast covariate value to numeric: " + token);
            }
          }
          i++;
        }

      } else { // read file header
        i = 0;

        while (std::getline(token_stream, token, sep)) {
          spdlog::debug("Phenotype file header n={}: {}", i, token);

          // record file column numbers of model variables
          if (std::find(covariate_column_headers.begin(), covariate_column_headers.end(), token)
              != covariate_column_headers.end()) {
            cov_idx.push_back(i);
            covariate_column.emplace_back(); // instantiate v of v
            spdlog::debug("Found sample covariate index: {}", i);
          } else if (token == outcome_column_header) {
            out_idx = i;
            spdlog::debug("Found outcome index: {}", out_idx);
          } else if (token == id_column_header) {
            sid_idx = i;
            spdlog::debug("Found sample id index: {}", sid_idx);
          }

          i++;
        }

        if (out_idx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + outcome_column_header);
        }
        if (sid_idx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + id_column_header);
        }
        passed_first_line = true;

        // get minimum covariate index in vector
        for (int idx : cov_idx) {
          if (idx < min_cov_idx) {
            min_cov_idx = idx;
          }
        }

      }

    }

    file.close();
  } else {
    throw std::runtime_error("Could not open file: " + pheno_file_path);
  }
  n_samples = sample_identifier_column.size();
  spdlog::info("Found {} samples in phenotype file", n_samples);
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
  for (unsigned i = 0; i < sample_identifier_column.size(); ++i) {
    if (mapping.count(sample_identifier_column[i])) {
      throw std::runtime_error("Duplicate identifiers were detected for: " + sample_identifier_column[i]);
    }
    // store array index against column identifier
    mapping[sample_identifier_column[i]] = i;
  }

  // create tmp data stores
  std::vector<std::string> sample_identifier_column_tmp;
  std::vector<double> outcome_column_tmp;
  std::vector<std::vector<double>> covariate_column_tmp;
  for (unsigned i = 0; i < covariate_column.size(); ++i) {
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
      for (unsigned j = 0; j < covariate_column.size(); ++j) {
        covariate_column_tmp[j].push_back(-1);
      }
    } else {
      // subset data
      unsigned idx = mapping[samples[i]];
      sample_identifier_column_tmp.push_back(sample_identifier_column[idx]);
      outcome_column_tmp.push_back(outcome_column[idx]);
      for (unsigned j = 0; j < covariate_column.size(); ++j) {
        covariate_column_tmp[j].push_back(covariate_column[j][idx]);
      }
      non_null_idx.insert(i);
    }

  }

  // check variables are the same length
  assert(sample_identifier_column_tmp.size() == samples.size());
  assert(outcome_column_tmp.size() == samples.size());
  for (unsigned i = 0; i < covariate_column.size(); ++i) {
    assert(covariate_column_tmp[i].size() == samples.size());
  }

  // set variables
  sample_identifier_column = sample_identifier_column_tmp;
  outcome_column = outcome_column_tmp;
  covariate_column = covariate_column_tmp;
  n_samples = sample_identifier_column.size();

  spdlog::info("Remaining samples after join: {}", n_samples);
  spdlog::info("Samples with non-null phenotypes: {}", non_null_idx.size());

  return (non_null_idx);
}
const std::vector<std::string> &PhenotypeFile::GetSampleIdentifierColumn() const {
  return sample_identifier_column;
}
const std::vector<double> &PhenotypeFile::GetOutcomeColumn() const {
  return outcome_column;
}
const std::vector<std::vector<double>> &PhenotypeFile::GetCovariateColumn() const {
  return covariate_column;
}
int PhenotypeFile::GetNSamples() const {
  if (n_samples == -1) {
    throw std::runtime_error("Value not set");
  }
  return n_samples;
}
}