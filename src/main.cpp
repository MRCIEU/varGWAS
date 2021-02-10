//
// Created by Matthew Lyon on 19/03/2020.
//

#include <iostream>
#include <cassert>
#include <stdexcept>
#include <thread>
#include <cxxopts.hpp>
#include "genfile/bgen/bgen.hpp"
#include <glog/logging.h>
#include "ThreadPool.h"
#include "BgenParser.h"
#include "PhenotypeFile.h"

bool file_exists(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int main(int argc, char **argv) {
  // Initialize Google's logging library.
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
      ("t,threads", "Number of threads", cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())));
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

  // Read phenotype data
  try {
    jlst::PhenotypeFile phenotype_file(variable_file, covariates, phenotype, id, sep);
    phenotype_file.parse();
  } catch (std::runtime_error const &e) {
    LOG(FATAL) << "Error parsing phenotype file: " << e.what();
    return -1;
  }

  // Read sample list from BGEN
  try {
    genfile::bgen::BgenParser bgenParser(bgen_file);
    static std::vector<std::string> samples;
    bgenParser.get_sample_ids(
        [](std::string const &id) { samples.push_back(id); }
    );
  } catch (genfile::bgen::BGenError const &e) {
    LOG(FATAL) << "Error parsing BGEN file: " << e.what();
    return -1;
  }

  // Extract samples from phenotype data


    // create thread pool with N worker threads
  LOG(INFO) << "Running with " << threads << " threads";
  ThreadPool pool(threads);

  try {
    // Parse BGEN
    LOG(INFO) << "Reading variants from: " << bgen_file;
    genfile::bgen::BgenParser bgenParser(bgen_file);

    // Parse sample list
    static std::vector<std::string> samples;
    bgenParser.get_sample_ids(
        [](std::string const &id) { samples.push_back(id); }
    );

    // Read variants
    std::string chromosome;
    uint32_t position;
    std::string rsid;
    std::vector<std::string> alleles;
    std::vector<std::vector<double> > probs;
    static std::vector<double> dosages;

    while (bgenParser.read_variant(&chromosome, &position, &rsid, &alleles)) {
      LOG_EVERY_N(INFO, 1000) << "Read the " << google::COUNTER << "th variant";

      // only support bi-allelic variants
      if (alleles.size() != 2){
        LOG(WARNING) << "Skipping non bi-allelic variant: " << rsid;
        continue;
      }

      // print variant
      /*std::cout << chromosome << '\t'
                << position << '\t'
                << rsid << '\t'
                << alleles[0] << '\t'
                << alleles[1] << '\t';*/

      // convert probabilities to dosage values
      bgenParser.read_probs(&probs);
      dosages.clear();
      for (auto &prob : probs) {
        // only support bi-allelic variants [0, 1, 2 copies of alt]
        assert(prob.size() == 3);

        // convert genotype probabilities to copies of alt
        // TODO check for null values
        dosages.push_back(prob[1] + (2 * prob[2]));
      }

      // check no missing values between sample list and dosage
      assert(dosages.size() == samples.size());

      // TODO add dosage values to Eigen matrix

      // enqueue and store future
      // auto assoc = pool.enqueue([](int answer) { return answer; }, 42);
    }

    // get result from future
    // std::cout << assoc.get() << std::endl;

    return 0;
  } catch (genfile::bgen::BGenError const &e) {
    LOG(FATAL) << "Error parsing BGEN file: " << e.what();
    return -1;
  }

  // TODO write output to CSV

}