//
// Created by Matthew Lyon on 19/03/2020.
//

#include <iostream>
#include <stdexcept>
#include <thread>
#include <cxxopts.hpp>
#include "genfile/bgen/bgen.hpp"
#include <glog/logging.h>
#include "ThreadPool.h"
#include "BgenParser.h"
#include "PhenotypeFile.h"
#include "PhenotypeFileException.h"
#include "Model.h"

bool file_exists(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int main(int argc, char **argv) {
  // Initialize Google's logging library.
  // TODO fix logging
  google::InitGoogleLogging(argv[0]);

  // Configure arguments
  cxxopts::Options options("JLST C++", "Program to perform vGWAS of trait against variants in the BGEN format");
  options.add_options()
      ("v,variable_file", "Path to phenotype file", cxxopts::value<std::string>())
      ("s,sep", "File separator", cxxopts::value<char>())
      ("c,covariates",
       "List of covariates column names separated by a comma (whitespace and quotes are not permitted).",
       cxxopts::value<std::vector<std::string>>())
      ("o,output_file", "Path to output file", cxxopts::value<std::string>())
      ("b,bgen_file", "Path to BGEN file", cxxopts::value<std::string>())
      ("p,phenotype", "Column name for phenotype", cxxopts::value<std::string>())
      ("i,id", "Column name for genotype identifier", cxxopts::value<std::string>())
      ("t,threads",
       "Number of threads",
       cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())));
  auto result = options.parse(argc, argv);

  // Parse arguments
  std::string variable_file = result["variable_file"].as<std::string>();
  char sep = result["sep"].as<char>();
  std::vector<std::string> covariates;
  if (result.count("covariates")) {
    covariates = result["covariates"].as<std::vector<std::string>>();
  }
  std::string output_file = result["output_file"].as<std::string>();
  std::string bgen_file = result["bgen_file"].as<std::string>();
  std::string phenotype = result["phenotype"].as<std::string>();
  std::string id = result["id"].as<std::string>();
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

  try {
    // Read BGEN
    genfile::bgen::BgenParser bgen_parser(bgen_file);
    static std::vector<std::string> samples;
    bgen_parser.get_sample_ids(
        [](std::string const &id) { samples.push_back(id); }
    );

    // Read phenotypes
    jlst::PhenotypeFile phenotype_file(variable_file, covariates, phenotype, id, sep);
    phenotype_file.parse();
    phenotype_file.subset_samples(samples);

    // create thread pool with N worker threads
    LOG(INFO) << "Running with " << threads << " threads";
    //ThreadPool pool(threads);

    // Perform locus association tests
    LOG(INFO) << "Running model";
    jlst::Model::run(phenotype_file, bgen_parser);

    // write output to CSV
    // TODO

  } catch (jlst::PhenotypeFileException const &e) {
    LOG(FATAL) << "Error parsing phenotype file. " << e.what();
    return -1;
  } catch (genfile::bgen::BGenError const &e) {
    LOG(FATAL) << "Error parsing BGEN file: " << e.what();
    return -1;
  }

}