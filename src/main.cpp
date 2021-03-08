//
// Created by Matthew Lyon on 19/03/2020.
//

#include <vector>
#include <string>
#include <set>
#include <memory>
#include <cxxopts.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/cfg/env.h"
#include "PhenotypeFile.h"
#include "Model.h"
#include "PhenotypeFileException.h"

bool file_exists(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int main(int argc, char **argv) {
  static std::string VERSION = "v1.0.0";
  static std::string PROGRAM_NAME = "JLST C++";
  spdlog::cfg::load_env_levels();
  static bool no_args = (argc == 1);

  // Configure arguments
  cxxopts::Options
      options(PROGRAM_NAME + " " + VERSION, "Program to perform vGWAS of trait against variants in the BGEN format");
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
      ("h,help", "Print usage")
      ("t,threads",
       "Number of threads",
       cxxopts::value<int>()->default_value(std::to_string(std::thread::hardware_concurrency())));
  auto result = options.parse(argc, argv);

  if (no_args || result.count("help") == 1) {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  if (result.count("variable_file") == 0) {
    spdlog::error("Phenotype file not provided");
    exit(1);
  }

  if (result.count("sep") == 0) {
    spdlog::error("Phenotype file separator not provided");
    exit(1);
  }

  if (result.count("output_file") == 0) {
    spdlog::error("Output file not provided");
    exit(1);
  }

  if (result.count("bgen_file") == 0) {
    spdlog::error("BGEN file not provided");
    exit(1);
  }

  if (result.count("phenotype") == 0) {
    spdlog::error("Column name of the phenotype not provided");
    exit(1);
  }

  if (result.count("id") == 0) {
    spdlog::error("Column name of the sample identifier not provided");
    exit(1);
  }

  // Parse arguments
  spdlog::info(PROGRAM_NAME + " " + VERSION);
  spdlog::debug("Log level: {}", spdlog::get_level());
  std::string variable_file = result["variable_file"].as<std::string>();
  spdlog::debug("Variable file: {}", variable_file);
  char sep = result["sep"].as<char>();
  spdlog::debug("Variable file sep: {}", sep);
  std::set<std::string> covariates;
  if (result.count("covariates") > 0) {
    std::vector<std::string> tmp = result["covariates"].as<std::vector<std::string>>();
    for (auto &e : tmp) {
      covariates.insert(e);
    }
  }
  for (auto &c:covariates) {
    spdlog::debug("Covariate(s): {}", c);
  }
  std::string output_file = result["output_file"].as<std::string>();
  spdlog::debug("Output file: {}", output_file);
  std::string bgen_file = result["bgen_file"].as<std::string>();
  spdlog::debug("BGEN file: {}", bgen_file);
  std::string phenotype = result["phenotype"].as<std::string>();
  spdlog::debug("Phenotype column: {}", phenotype);
  std::string id = result["id"].as<std::string>();
  spdlog::debug("ID column: {}", id);
  int threads = result["threads"].as<int>();
  spdlog::debug("Threads n={}", threads);

  // check files exist
  if (!file_exists(variable_file)) {
    spdlog::error("File does not exist or is not readable: {}", variable_file);
    return -1;
  }
  if (!file_exists(bgen_file)) {
    spdlog::error("File does not exist or is not readable: {}", bgen_file);
    return -1;
  }
  if (threads < 1) {
    spdlog::error("Invalid value for threads: {}", threads);
    return -1;
  }

  try {
    // Open BGEN and read sample list
    spdlog::info("Reading samples from BGEN: {}", bgen_file);
    genfile::bgen::BgenParser bgen_parser(bgen_file);
    static std::vector<std::string> samples;
    bgen_parser.get_sample_ids(
        [](std::string const &id) { samples.push_back(id); }
    );

    // Read phenotypes and subset using sample list
    jlst::PhenotypeFile phenotype_file(variable_file, covariates, phenotype, id, sep);
    phenotype_file.parse();
    std::set<unsigned> non_null_idx = phenotype_file.join(samples);

    // Perform locus association tests & write to file
    jlst::Model model(phenotype_file, bgen_parser, non_null_idx, output_file, threads);
    model.run();

    spdlog::info("Analysis complete!");

  } catch (jlst::PhenotypeFileException const &e) {
    spdlog::error("Error parsing phenotype file: {}", e.what());
    return -1;
  } catch (genfile::bgen::BGenError const &e) {
    spdlog::error("Error parsing BGEN file: {}", e.what());
    return -1;
  }

  return 0;
}