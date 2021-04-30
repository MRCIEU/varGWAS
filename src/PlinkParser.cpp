//
// Created by Matt Lyon on 30/04/2021.
//

#include <vector>
#include <string>
#include <plinkio/plinkio.h>
#include <stdexcept>
#include "PlinkParser.h"

std::vector<std::string> jlst::PlinkParser::get_samples() {
  struct pio_file_t plink_file{};
  int sample_id;
  std::vector<std::string> samples;

  if (pio_open(&plink_file, _file_path.c_str()) != PIO_OK) {
    throw std::runtime_error("Could not open file: " + _file_path);
  }

  if (!pio_one_locus_per_row(&plink_file)) {
    throw std::runtime_error("This script requires that snps are rows and samples columns");
  }

  for (sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++) {
    struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
    samples.emplace_back(sample->iid);
  }

  pio_close(&plink_file);

  return(samples);
}