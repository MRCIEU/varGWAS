#include "gtest/gtest.h"
#include "PlinkParser.h"

TEST(PlinkTest, read_bed) {
  std::string file_name = "/Users/ml18692/projects/jlst_cpp/lib/libplinkio-0.9.8/tests/data/wgas";
  jlst::PlinkParser plink_parser(file_name);
  for (auto &s : plink_parser.get_samples()) {
    ASSERT_EQ(s, "NA18526");
    break;
  }
}