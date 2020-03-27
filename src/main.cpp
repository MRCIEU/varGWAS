//
// Created by Matthew Lyon on 19/03/2020.
//

#include <iostream>
#include <cassert>
#include <stdexcept>
#include <cxxopts.hpp>
#include "genfile/bgen/bgen.hpp"
#include "csv.h"
#include <glog/logging.h>
#include "BgenParser.h"

bool file_exists(const std::string &name) {
    std::ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char **argv) {
    // Initialize Google's logging library.
    google::InitGoogleLogging(argv[0]);

    // Parse arguments
    cxxopts::Options options("JLST C++", "Program to perform vGWAS of trait against variants in the BGEN format");
    options.add_options()
            ("v,variable_file", "Path to phenotype file", cxxopts::value<std::string>())
            ("c,covariates", "List of covariates column names separated by a comma (whitespace is not permitted).",
             cxxopts::value<std::vector<std::string>>())
            ("o,output_file", "Path to output file", cxxopts::value<std::string>())
            ("b,bgen_file", "Path to BGEN file", cxxopts::value<std::string>())
            ("p,phenotype", "Column name for phenotype", cxxopts::value<std::string>())
            ("t,threads", "Number of threads", cxxopts::value<int>()->default_value("1"));
    auto result = options.parse(argc, argv);

    std::string variable_file = result["variable_file"].as<std::string>();
    std::vector<std::string> covariates = result["covariates"].as<std::vector<std::string>>();
    std::string output_file = result["output_file"].as<std::string>();
    std::string bgen_file = result["bgen_file"].as<std::string>();
    std::string phenotype = result["phenotype"].as<std::string>();
    int threads = result["threads"].as<int>();

    // check files exist
    if (!file_exists(variable_file)) {
        LOG(FATAL) << "File does not exist or is not readable: " << variable_file;
        return -1;
    }
    if (!file_exists(bgen_file)) {
        LOG(FATAL) << "File does not exist or is not readable: " << bgen_file;
        return -1;
    }

    // Read phenotype and covariates into memory
    LOG(INFO) << "Reading variables from: " << variable_file;
    io::LineReader in(variable_file);
    while (char *line = in.next_line()) {
        std::cout << line << std::endl;
    }

    // TODO initalize static JSLP obj with Y and covar

    // TODO init multiple threads
    LOG(INFO) << "Running with " << threads << " threads";

    // TODO foreach variant estimate SNP-trait effect adjusted for covar one on each thread

    // TODO write out SNP-trait assoc to CSV using mutex?

    // TODO check file exists before passing to BgenParser

    try {
        LOG(INFO) << "Reading variants from: " << bgen_file;
        genfile::bgen::BgenParser bgenParser(bgen_file);

        bgenParser.summarise(std::cerr);

        // Output header
        std::cout << "CHROM\tPOS\tID\tREF\tALT";
        bgenParser.get_sample_ids(
                [](std::string const &id) { std::cout << "\t" << id; }
        );
        std::cout << "\n";

        // Output variants
        std::string chromosome;
        uint32_t position;
        std::string rsid;
        std::vector<std::string> alleles;
        std::vector<std::vector<double> > probs;

        while (bgenParser.read_variant(&chromosome, &position, &rsid, &alleles)) {
            LOG_EVERY_N(INFO, 10) << "Got the " << google::COUNTER << "th cookie";

            // only support bi-allelic variants
            assert(alleles.size() == 2);

            // print variant
            std::cout << chromosome << '\t'
                      << position << '\t'
                      << rsid << '\t'
                      << alleles[0] << '\t'
                      << alleles[1] << '\t';

            bgenParser.read_probs(&probs);
            for (auto &prob : probs) {
                std::cout << '\t';

                // only support bi-allelic variants [0, 1, 2 copies of alt]
                assert(prob.size() == 3);

                // convert genotype probabilities to copies of alt
                std::cout << prob[1] + (2 * prob[2]);
            }

            std::cout << "\n";
        }

        return 0;
    } catch (genfile::bgen::BGenError const &e) {
        std::cerr << "Error parsing bgen file.\n";
        return -1;
    }

}