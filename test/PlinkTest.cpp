#include "gtest/gtest.h"
#include <plinkio/plinkio.h>

TEST(PlinkTest, read_bed) {
  std::string file_name = "/Users/ml18692/projects/jlst_cpp/lib/libplinkio-0.9.8/tests/data/wgas";
  struct pio_file_t plink_file{};
  snp_t *snp_buffer;
  int sample_id;
  int locus_id;

  if (pio_open(&plink_file, file_name.c_str()) != PIO_OK) {
    throw std::runtime_error("Could not open file: " + file_name);
  }

  if (!pio_one_locus_per_row(&plink_file)) {
    throw std::runtime_error("This script requires that snps are rows and samples columns");
  }

  for (sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++) {
    struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
    printf("%s\n", sample->iid);
  }

  locus_id = 0;
  snp_buffer = (snp_t *) malloc(pio_row_size(&plink_file));
  while (pio_next_row(&plink_file, snp_buffer) == PIO_OK) {
    for (sample_id = 0; sample_id < pio_num_samples(&plink_file); sample_id++) {
      struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
      struct pio_locus_t *locus = pio_get_locus(&plink_file, locus_id);
      printf("Individual %s has genotype %d for snp %s.\n", sample->iid, snp_buffer[sample_id], locus->name);
    }

    locus_id++;
  }

  free(snp_buffer);
  pio_close(&plink_file);
}