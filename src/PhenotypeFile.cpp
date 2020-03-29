#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>

/*
 * Class to read phenotype CSV
 *
 * */

namespace jlst {
    class PhenotypeFile {
    public:
        PhenotypeFile(const std::string &phenoFile, const std::vector<std::string> &covariateFields,
                      const std::string &outcomeField, const std::string &sep) {
            this->phenoFile = phenoFile;
            this->covariateFields = covariateFields;
            this->outcomeField = outcomeField;
            this->sep = sep;
        };

        void readfile() {
            std::ifstream file(phenoFile);
            if (file.is_open()) {
                std::string line;
                while (getline(file, line)) {

                }
                file.close();
            }
        };
    private:
        std::string phenoFile;
        std::vector<std::string> covariateFields;
        std::string outcomeField;
        std::string sep;
    };
}