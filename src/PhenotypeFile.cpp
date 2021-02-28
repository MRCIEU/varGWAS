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
PhenotypeFile::PhenotypeFile(const std::string &phenoFilePath,
                             const std::vector<std::string> &covariateColumnHeaders,
                             const std::string &outcomeColumnHeader,
                             const std::string &idColumnHeader,
                             const char &sep) {
  this->phenoFilePath = phenoFilePath;
  this->covariateColumnHeaders = covariateColumnHeaders;
  this->outcomeColumnHeader = outcomeColumnHeader;
  this->idColumnHeader = idColumnHeader;
  this->sep = sep;
  this->n_samples = -1;
};

/*
 * Function to parse file
 * */
void PhenotypeFile::parse() {
  spdlog::info("Parsing phenotype from: {}", phenoFilePath);
  std::ifstream file(phenoFilePath.c_str());
  int outIdx = -1;
  int sidIdx = -1;
  std::vector<int> covIdx;
  int i;
  unsigned minCovIdx = UINT_MAX;

  if (file.is_open()) {
    bool passedFirstLine = false;
    std::string line;
    std::string token;

    // read file line-by-line
    while (getline(file, line)) {
      std::istringstream tokenStream(line);

      if (passedFirstLine) { // read file body
        i = 0;

        // TODO implement using boost to allow for quotes in the file
        while (std::getline(tokenStream, token, sep)) {
          if (i == outIdx) {
            outcomeColumn.push_back(std::stod(token));
          }
          if (i == sidIdx) {
            sampleIdentifierColumn.push_back(token);
          }
          if (std::find(covIdx.begin(), covIdx.end(), i) != covIdx.end()) {
            covariateColumn[i - minCovIdx].push_back(std::stod(token));
          }
          i++;
        }

      } else { // read file header
        i = 0;

        while (std::getline(tokenStream, token, sep)) {

          // record file column numbers of model variables
          if (std::find(covariateColumnHeaders.begin(), covariateColumnHeaders.end(), token)
              != covariateColumnHeaders.end()) {
            covIdx.push_back(i);
            covariateColumn.emplace_back(); // instantiate v of v
          } else if (token == outcomeColumnHeader) {
            outIdx = i;
          } else if (token == idColumnHeader) {
            sidIdx = i;
          }

          i++;
        }

        if (outIdx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + outcomeColumnHeader);
        }
        if (sidIdx == -1) {
          throw jlst::PhenotypeFileException("Field missing from phenotype file: " + idColumnHeader);
        }
        passedFirstLine = true;

        // get minimum covariate index in vector
        for (auto &idx : covIdx) {
          if (idx < minCovIdx) {
            minCovIdx = idx;
          }
        }

      }

    }

    file.close();
  } else {
    throw std::runtime_error("Could not open file: " + phenoFilePath);
  }
  n_samples = sampleIdentifierColumn.size();
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
  for (unsigned i = 0; i < sampleIdentifierColumn.size(); ++i) {
    if (mapping.count(sampleIdentifierColumn[i])) {
      throw std::runtime_error("Duplicate identifiers were detected for: " + sampleIdentifierColumn[i]);
    }
    // store array index against column identifier
    mapping[sampleIdentifierColumn[i]] = i;
  }

  // create tmp data stores
  std::vector<std::string> sampleIdentifierColumnTmp;
  std::vector<double> outcomeColumnTmp;
  std::vector<std::vector<double>> covariateColumnTmp;
  for (unsigned i = 0; i < covariateColumn.size(); ++i) {
    covariateColumnTmp.emplace_back();
  }

  // subset data using new ordering
  for (unsigned i = 0; i < samples.size(); ++i) {

    // get index for sample
    if (mapping.count(samples[i]) == 0) {
      spdlog::info("Sample missing from phenotype file: {}", samples[i]);
      // add missing data with null values
      sampleIdentifierColumnTmp.push_back(samples[i]);
      outcomeColumnTmp.push_back(-1);
      for (unsigned j = 0; j < covariateColumn.size(); ++j) {
        covariateColumnTmp[j].push_back(-1);
      }
    } else {
      // subset data
      unsigned idx = mapping[samples[i]];
      sampleIdentifierColumnTmp.push_back(sampleIdentifierColumn[idx]);
      outcomeColumnTmp.push_back(outcomeColumn[idx]);
      for (unsigned j = 0; j < covariateColumn.size(); ++j) {
        covariateColumnTmp[j].push_back(covariateColumn[j][idx]);
      }
      non_null_idx.insert(i);
    }

  }

  // check variables are the same length
  assert(sampleIdentifierColumnTmp.size() == samples.size());
  assert(outcomeColumnTmp.size() == samples.size());
  for (unsigned i = 0; i < covariateColumn.size(); ++i) {
    assert(covariateColumnTmp[i].size() == samples.size());
  }

  // set variables
  sampleIdentifierColumn = sampleIdentifierColumnTmp;
  outcomeColumn = outcomeColumnTmp;
  covariateColumn = covariateColumnTmp;
  n_samples = sampleIdentifierColumn.size();

  spdlog::info("Remaining samples after subset: {}", n_samples);
  spdlog::info("Samples with non-null phenotypes: {}", non_null_idx.size());

  return (non_null_idx);
}
const std::vector<std::string> &PhenotypeFile::GetSampleIdentifierColumn() const {
  return sampleIdentifierColumn;
}
const std::vector<double> &PhenotypeFile::GetOutcomeColumn() const {
  return outcomeColumn;
}
const std::vector<std::vector<double>> &PhenotypeFile::GetCovariateColumn() const {
  return covariateColumn;
}
int PhenotypeFile::GetNSamples() const {
  if (n_samples == -1) {
    throw std::runtime_error("Value not set");
  }
  return n_samples;
}
}