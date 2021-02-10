#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <glog/logging.h>
#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"

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
};

/*
 * Function to parse file
 * */
void PhenotypeFile::parse() {
  LOG(INFO) << "Parsing phenotype from: " << phenoFilePath;
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

}
/*
 * Function to subset data using provided list of sample identifiers
 * */
void PhenotypeFile::subset_samples(const std::vector<std::string> &samples) {
  // create sample identifier-to-index mapping
  std::unordered_map<std::string, unsigned> mapping;
  for (unsigned i = 0; i < sampleIdentifierColumn.size(); ++i) {
    if (mapping.count(sampleIdentifierColumn[i])){
      throw std::runtime_error("Duplicate identifiers were detected for: " + sampleIdentifierColumn[i]);
    }
    // store array index against column identifier
    mapping[sampleIdentifierColumn[i]] = i;
  }

  // create tmp data stores
  std::vector<std::string> sampleIdentifierColumnTmp;
  std::vector<double> outcomeColumnTmp;
  std::vector<std::vector<double>> covariateColumnTmp;
  for (unsigned i = 0; i < covariateColumn.size(); ++i){
    covariateColumnTmp.emplace_back();
  }

  // subset data using new ordering
  for (auto &sample : samples) {
    // get index for sample
    if (mapping.count(sample) == 0){
      throw std::runtime_error("Missing sample from phenotype file: " + sample);
    }

    // subset data
    unsigned idx = mapping[sample];
    sampleIdentifierColumnTmp.push_back(sampleIdentifierColumn[idx]);
    outcomeColumnTmp.push_back(outcomeColumn[idx]);
    for (unsigned i = 0; i < covariateColumn.size(); ++i){
      covariateColumnTmp[i].push_back(covariateColumn[i][idx]);
    }
  }

  // check variables are the same length
  assert(sampleIdentifierColumnTmp.size() == samples.size());
  assert(outcomeColumnTmp.size() == samples.size());
  for (unsigned i = 0; i < covariateColumn.size(); ++i){
    assert(covariateColumnTmp[i].size() == samples.size());
  }

  // set variables
  sampleIdentifierColumn = sampleIdentifierColumnTmp;
  outcomeColumn = outcomeColumnTmp;
  covariateColumn = covariateColumnTmp;
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
}