#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>
#include <glog/logging.h>
#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"

/*
 * Class to read in outcome, covariate(s) and sample identifier into memory
 * */

namespace jlst {
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

// TODO implement using boost to allow for quotes in the file
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
          throw PhenotypeFileException("Field missing from phenotype file: " + outcomeColumnHeader);
        }
        if (sidIdx == -1) {
          throw PhenotypeFileException("Field missing from phenotype file: " + idColumnHeader);
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