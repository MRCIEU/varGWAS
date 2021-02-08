#include <stdexcept>
#include <vector>
#include <iostream>
#include <glog/logging.h>
#include "PhenotypeFile.h"

/*
 * Class to read in phenotype and covariate(s) into memory
 * */

namespace jlst {
PhenotypeFile::PhenotypeFile(const std::string &phenoFilePath,
                             const std::vector<std::string> &covariateColumnHeaders,
                             const std::string &outcomeColumnHeader,
                             const char &sep) {
    this->phenoFilePath = phenoFilePath;
    this->covariateColumnHeaders = covariateColumnHeaders;
    this->outcomeColumnHeader = outcomeColumnHeader;
    this->sep = sep;
};

void PhenotypeFile::load() {
    // TODO implement using boost to allow for quotes in the file
    // TODO improve performance
    LOG(INFO) << "Parsing phenotype from: " << phenoFilePath;
    static std::ifstream file(phenoFilePath.c_str());

    if (file.is_open()) {
        static bool passedFirstLine = false;
        static std::string line;
        static std::string token;

        while (getline(file, line)) {
            std::istringstream tokenStream(line);

            if (passedFirstLine) {
                static std::vector<double> fields;
                while (std::getline(tokenStream, token, sep)) {
                    fields.push_back(std::stod(token));
                }
                fileBody.push_back(fields);
            } else {
                while (std::getline(tokenStream, token, sep)) {
                    fileHeader.push_back(token);
                }
                passedFirstLine = true;
            }

        }

        file.close();
    } else {
        throw std::runtime_error("Could not open file: " + phenoFilePath);
    }

}
const std::vector<std::string> &PhenotypeFile::GetCovariateColumnHeaders() const {
    return covariateColumnHeaders;
}
const std::string &PhenotypeFile::GetOutcomeColumnHeader() const {
    return outcomeColumnHeader;
}
const std::vector<std::string> &PhenotypeFile::GetFileHeader() const {
    return fileHeader;
}
const std::vector<std::vector<double>> &PhenotypeFile::GetFileBody() const {
    return fileBody;
};
}