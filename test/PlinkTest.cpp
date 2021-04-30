#include "gtest/gtest.h"
#include "PlinkParser.h"
extern "C" {
#include <plinkio/plinkio.h>
}

TEST(PlinkTest, read_bed) {
  std::string file_name = "/Users/ml18692/projects/jlst_cpp/lib/libplinkio-0.9.8/tests/data/wgas";
  jlst::PlinkParser plink_parser(file_name);
  for (auto &s : plink_parser.get_samples()) {
    ASSERT_EQ(s, "NA18526");
    break;
  }
}

TEST(PlinkTest, model){
  std::string file_path = "/Users/ml18692/projects/jlst_cpp/lib/libplinkio-0.9.8/tests/data/wgas";
  std::string chromosome;
  uint32_t position;
  std::string rsid;
  std::vector<std::string> alleles;
  std::vector<double> dosages;
  struct pio_file_t plink_file;
  snp_t *snp_buffer;
  int sample_id;
  int locus_id;

  // Read variant-by-variant
  if (pio_open(&plink_file, file_path.c_str()) != PIO_OK) {
    throw std::runtime_error("Could not open file: " + file_path);
  }

  if (!pio_one_locus_per_row(&plink_file)) {
    throw std::runtime_error("This script requires that snps are rows and samples columns");
  }

  locus_id = 0;
  snp_buffer = (snp_t *) malloc(pio_row_size(&plink_file));

  // loop over variants
  while (pio_next_row(&plink_file, snp_buffer) == PIO_OK) {
    struct pio_locus_t *locus = pio_get_locus(&plink_file, locus_id);
    int chr = locus->chromosome;
    chromosome = std::to_string(chr);
    position = locus->bp_position;
    rsid = locus->name;
    alleles.emplace_back(locus->allele1);
    alleles.emplace_back(locus->allele2);

    // loop over samples
    for (sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++) {
      struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
      int dosage = snp_buffer[sample_id];
      if (dosage == 3) {
        dosages.push_back(-1.0);
      } else {
        dosages.push_back((double) dosage);
      }
    }
    locus_id++;
    dosages.clear();
    alleles.clear();
  }

  free(snp_buffer);
  pio_close(&plink_file);
}